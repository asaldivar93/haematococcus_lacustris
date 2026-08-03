"""
Utility for checking whether InChIKeys resolve to valid compounds in UniChem.

UniChem Compound Search API docs:
https://chembl.gitbook.io/unichem/api/compound-search

Entry points:

- check_inchikeys_unichem(...)          - simple synchronous version (requests),
                                          fine for a handful to a few hundred keys.
- check_inchikeys_unichem_bulk(...)     - async (aiohttp) with bounded
                                          concurrency, token-bucket rate limiting,
                                          progress logging and crash-safe
                                          checkpointing. Use for thousands of keys.
- check_inchikeys_unichem_bulk_sync(...) - blocking wrapper around the above,
                                          for calling from ordinary code.
- bulk_fetch_json(specs, ...)           - the same engine, generalised: any list
                                          of independent HTTP calls.

ASYNC IS THE DEFAULT FOR SLOW APIS IN THIS SKILL.
Every one of these calls is pure network wait, so concurrency is the only
thing that matters and asyncio does it with the least machinery: one event
loop, one connection pool, a semaphore for in-flight bounding and a token
bucket for the polite rate. Do not reach for threads or multiprocessing here
— a thread pool pays a stack and a GIL handoff per request to sit on a
socket, and every HTTP client session then needs its own pool because it is
not thread-safe. Route new bulk API work through `bulk_fetch_json` rather
than writing another parallelisation.
"""

import asyncio
import json
import os
import time
from typing import Any, Dict, List, Optional

import requests


def check_inchikeys_unichem(
    inchikeys: List[str],
    timeout: int = 15,
    max_retries: int = 3,
    retry_backoff: float = 2.0,
    request_delay: float = 0.0,
) -> Dict[str, bool]:
    """
    Query the UniChem Compound Search API (v1) for a list of InChIKeys and
    report whether each one resolves to a valid, known compound.

    UniChem returns HTTP 200 with an empty "compounds" list (and
    "response": "Not found") for InChIKeys that are malformed, empty, or
    simply not present in any of its source databases. So "found" is
    determined by whether "compounds" is non-empty, not by HTTP status.

    Parameters
    ----------
    inchikeys : list of str
        InChIKeys to look up (e.g. "BSYNRYMUTXBXSQ-UHFFFAOYSA-N").
    timeout : int
        Per-request timeout in seconds.
    max_retries : int
        Number of retry attempts for transient network/server errors
        (timeouts, connection errors, 5xx responses).
    retry_backoff : float
        Base seconds for exponential backoff between retries
        (attempt 1 waits retry_backoff, attempt 2 waits retry_backoff*2, ...).
    request_delay : float
        Optional pause (seconds) between successive requests, to be polite
        to the API when checking long lists.

    Returns
    -------
    dict
        Mapping {inchikey: bool}, True if UniChem returned at least one
        matching compound, False if not found (or if the lookup ultimately
        failed after retries -- in which case a warning is printed).
    """
    url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
    headers = {"Content-Type": "application/json"}
    results: Dict[str, bool] = {}

    for key in inchikeys:
        payload = {"type": "inchikey", "compound": key}
        found: Optional[bool] = None

        for attempt in range(max_retries):
            try:
                resp = requests.post(url, json=payload, headers=headers, timeout=timeout)
                if resp.status_code == 200:
                    data = resp.json()
                    found = bool(data.get("compounds"))
                    break
                elif 500 <= resp.status_code < 600:
                    # transient server error -> retry
                    raise requests.exceptions.RequestException(
                        f"Server error {resp.status_code}"
                    )
                else:
                    # non-retryable client error (e.g. bad request format)
                    print(f"[unichem] {key}: unexpected status {resp.status_code} - {resp.text[:200]}")
                    found = False
                    break
            except requests.exceptions.RequestException as e:
                if attempt < max_retries - 1:
                    time.sleep(retry_backoff * (2 ** attempt))
                    continue
                print(f"[unichem] {key}: failed after {max_retries} attempts ({e}); marking as False")
                found = False

        results[key] = bool(found)

        if request_delay:
            time.sleep(request_delay)

    return results


class _RateLimiter:
    """
    Token-bucket style rate limiter for asyncio: `acquire()` blocks just long
    enough to keep the long-run call rate at or below `rate` per second,
    regardless of how many coroutines are trying to acquire concurrently.
    """

    def __init__(self, rate: float):
        self.interval = 1.0 / rate if rate and rate > 0 else 0.0
        self._lock = asyncio.Lock()
        self._next_time = 0.0

    async def acquire(self):
        if self.interval == 0:
            return
        async with self._lock:
            loop = asyncio.get_event_loop()
            now = loop.time()
            wait = self._next_time - now
            if wait > 0:
                await asyncio.sleep(wait)
                now = loop.time()
            self._next_time = max(now, self._next_time) + self.interval


async def _fetch_one_async(session, sem, limiter, key, url, timeout, max_retries, retry_backoff):
    payload = {"type": "inchikey", "compound": key}
    data = await _request_async(session, sem, limiter, "POST", url,
                                json_body=payload, params=None, timeout=timeout,
                                max_retries=max_retries,
                                retry_backoff=retry_backoff, label=key)
    return key, bool((data or {}).get("compounds"))


# ==========================================================================
# Generic async engine
# ==========================================================================
#
# The machinery below is what makes the bulk InChIKey check fast, and none of
# it is specific to InChIKeys: bounded concurrency, a shared token bucket,
# retry with backoff, crash-safe checkpointing and progress logging apply to
# any list of independent HTTP calls. It is factored out here so every slow
# API in this skill goes through one implementation instead of each one
# growing its own.
#
# Async rather than threads: these calls are pure network wait. A thread pool
# pays a stack and a GIL handoff per request to do nothing but block on a
# socket, and requests.Session is not thread-safe so each thread needs its
# own connection pool. One event loop with a semaphore holds thousands of
# in-flight requests on a single shared pool, and the rate limiter is a plain
# await with no lock contention.


async def _request_async(session, sem, limiter, method, url, *, json_body,
                         params, timeout, max_retries, retry_backoff, label):
    """One rate-limited request with retries. Returns parsed JSON or None.

    None means "no usable answer": a 404, a non-retryable client error, or
    exhausted retries. Callers that must distinguish "absent" from "failed"
    should inspect their own payload semantics — UniChem, for instance,
    answers 200 with an empty list for an unknown key.
    """
    import aiohttp

    async with sem:
        for attempt in range(max_retries):
            await limiter.acquire()
            try:
                async with session.request(
                    method, url, json=json_body, params=params,
                    timeout=aiohttp.ClientTimeout(total=timeout),
                ) as resp:
                    if resp.status == 200:
                        try:
                            return await resp.json(content_type=None)
                        except (ValueError, aiohttp.ClientError):
                            return None
                    if resp.status == 404:
                        return None
                    if resp.status == 429:
                        # Honour Retry-After when the server sends it.
                        retry_after = resp.headers.get("Retry-After")
                        wait = (float(retry_after) if retry_after
                                else retry_backoff * (2 ** attempt))
                        await asyncio.sleep(wait)
                        continue
                    if 500 <= resp.status < 600:
                        await asyncio.sleep(retry_backoff * (2 ** attempt))
                        continue
                    return None                       # 4xx: will not improve
            except (asyncio.TimeoutError, OSError, aiohttp.ClientError) as e:
                if attempt < max_retries - 1:
                    await asyncio.sleep(retry_backoff * (2 ** attempt))
                    continue
                # aiohttp raises several exceptions with empty str(); the
                # class name is the only diagnostic in those cases.
                detail = str(e) or type(e).__name__
                print(f"[async] {label}: failed after {max_retries} "
                      f"attempts ({detail})")
                return None
        return None


async def bulk_fetch_json(
    specs: List[dict],
    *,
    concurrency: int = 10,
    rate_per_sec: float = 10.0,
    timeout: int = 15,
    max_retries: int = 3,
    retry_backoff: float = 2.0,
    progress_every: int = 0,
    label: str = "async",
) -> Dict[str, Any]:
    """Run many independent HTTP calls concurrently; return {key: json}.

    Parameters
    ----------
    specs : list of dict
        One request each: ``{"key": str, "url": str, "method": "GET"|"POST",
        "json": ..., "params": ...}``. ``key`` is how the caller finds the
        result and must be unique; only ``key`` and ``url`` are required.
    concurrency : int
        Maximum simultaneous in-flight requests.
    rate_per_sec : float
        Sustained ceiling on the request rate across all workers. This is the
        courtesy knob for a shared public service — raising `concurrency`
        without raising this just queues requests, which is the intended
        behaviour.
    progress_every : int
        Print a progress line every N completions; 0 disables.

    Returns
    -------
    dict
        ``{key: parsed_json_or_None}``, one entry per spec.

    Notes
    -----
    Different hosts are *not* rate-limited separately here. Group specs by
    host and call once per host when the limits differ materially.
    """
    import aiohttp

    if not specs:
        return {}

    sem = asyncio.Semaphore(concurrency)
    limiter = _RateLimiter(rate_per_sec)
    connector = aiohttp.TCPConnector(limit=concurrency)
    out: Dict[str, Any] = {}
    total = len(specs)

    async with aiohttp.ClientSession(
        connector=connector,
        headers={"Accept": "application/json"},
        trust_env=True,                      # honour HTTP(S)_PROXY env vars
    ) as session:
        async def one(spec):
            data = await _request_async(
                session, sem, limiter,
                spec.get("method", "GET"), spec["url"],
                json_body=spec.get("json"), params=spec.get("params"),
                timeout=timeout, max_retries=max_retries,
                retry_backoff=retry_backoff, label=spec["key"])
            return spec["key"], data

        tasks = [asyncio.create_task(one(s)) for s in specs]
        done = 0
        for coro in asyncio.as_completed(tasks):
            key, data = await coro
            out[key] = data
            done += 1
            if progress_every and done % progress_every == 0:
                print(f"[{label}] {done}/{total} done ({done/total:.1%})")
    return out


def run_async(coro):
    """Run a coroutine from ordinary synchronous code.

    Kernels and scripts calling into these helpers usually have no event loop
    running. When one *is* running (inside a notebook that has already
    started a loop, or another coroutine) this raises rather than silently
    deadlocking — await the coroutine directly in that case.
    """
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)
    raise RuntimeError(
        "run_async() called from inside a running event loop; "
        "await the coroutine directly instead.")


def _load_checkpoint(path: Optional[str]) -> Dict[str, bool]:
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def _save_checkpoint(path: Optional[str], results: Dict[str, bool]):
    if not path:
        return
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(results, f)
    os.replace(tmp, path)  # atomic on POSIX, avoids truncated files on crash


async def _check_inchikeys_unichem_bulk_async(
    inchikeys: List[str],
    concurrency: int,
    rate_per_sec: float,
    timeout: int,
    max_retries: int,
    retry_backoff: float,
    checkpoint_path: Optional[str],
    checkpoint_every: int,
    progress_every: int,
) -> Dict[str, bool]:
    import aiohttp  # imported lazily so the sync function has no hard dependency on it

    url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"

    results: Dict[str, bool] = _load_checkpoint(checkpoint_path)
    todo = [k for k in inchikeys if k not in results]
    total = len(inchikeys)
    if results:
        print(f"[unichem] resuming from checkpoint: {len(results)}/{total} already done, {len(todo)} remaining")

    if not todo:
        return {k: results[k] for k in inchikeys}

    sem = asyncio.Semaphore(concurrency)
    limiter = _RateLimiter(rate_per_sec)
    connector = aiohttp.TCPConnector(limit=concurrency)

    completed_since_checkpoint = 0
    done_count = len(results)

    async with aiohttp.ClientSession(
        connector=connector,
        headers={"Content-Type": "application/json"},
        trust_env=True,  # honor HTTP(S)_PROXY env vars
    ) as session:
        tasks = [
            asyncio.create_task(
                _fetch_one_async(session, sem, limiter, k, url, timeout, max_retries, retry_backoff)
            )
            for k in todo
        ]
        for coro in asyncio.as_completed(tasks):
            key, found = await coro
            results[key] = found
            done_count += 1
            completed_since_checkpoint += 1

            if progress_every and done_count % progress_every == 0:
                print(f"[unichem] {done_count}/{total} done ({done_count/total:.1%})")

            if checkpoint_path and completed_since_checkpoint >= checkpoint_every:
                _save_checkpoint(checkpoint_path, results)
                completed_since_checkpoint = 0

    if checkpoint_path:
        _save_checkpoint(checkpoint_path, results)

    # {inchikey: bool}, matching the documented contract and the synchronous
    # sibling. (An earlier revision returned a list of records here, which
    # disagreed with both the annotation and the docstring.)
    return {k: results[k] for k in inchikeys}


async def check_inchikeys_unichem_bulk(
    inchikeys: list,
    concurrency: int = 10,
    rate_per_sec: float = 10.0,
    timeout: int = 15,
    max_retries: int = 3,
    retry_backoff: float = 2.0,
    checkpoint_path: str|None = None,
    checkpoint_every: int = 100,
    progress_every: int = 50,
) -> Dict[str, bool]:
    """Check thousands of InChIKeys against UniChem concurrently.

    With bounded parallelism, a token-bucket rate limit, progress logging,
    and optional crash-safe checkpointing.

    This is a synchronous function (safe to call from ordinary scripts / a
    Jupyter-style kernel) that runs an asyncio event loop internally.

    Parameters
    ----------
    inchikeys : list of str
        InChIKeys to look up. Duplicates are only queried once.
    concurrency : int
        Max number of simultaneous in-flight requests (bounds parallelism
        independently of the rate limit -- useful if the server is slow to
        respond but you still don't want to open hundreds of sockets).
    rate_per_sec : float
        Maximum sustained request rate, in requests/second, across all
        concurrent workers. This is the main "rate limiting courtesy" knob --
        keep it modest (UniChem has no published public rate limit, so
        5-10 req/s is a reasonable default for a public EBI service).
    timeout : int
        Per-request timeout in seconds.
    max_retries : int
        Retry attempts per key on timeout, connection error, HTTP 5xx, or
        HTTP 429 (429 additionally honors a `Retry-After` header if present).
    retry_backoff : float
        Base seconds for exponential backoff between retries.
    checkpoint_path : str, optional
        If given, results are periodically written to this JSON file
        (atomically, via a temp-file + os.replace). If the file already
        exists when the function is called, previously-checkpointed keys
        are loaded and skipped -- so an interrupted run of thousands of
        keys can simply be re-launched with the same checkpoint_path.
    checkpoint_every : int
        Write the checkpoint file after this many additional keys complete.
    progress_every : int
        Print a progress line after this many additional keys complete
        (set to 0 to disable).

    Returns
    -------
    dict
        Mapping {inchikey: bool}, in the same key set as the input
        (duplicates in the input map to the same bool).

    Example
    -------
    >>> keys = [...]  # thousands of InChIKeys
    >>> results = check_inchikeys_unichem_bulk(
    ...     keys,
    ...     concurrency=10,
    ...     rate_per_sec=10.0,
    ...     checkpoint_path="unichem_checkpoint.json",
    ... )
    """
    return await _check_inchikeys_unichem_bulk_async(
        inchikeys,
        concurrency=concurrency,
        rate_per_sec=rate_per_sec,
        timeout=timeout,
        max_retries=max_retries,
        retry_backoff=retry_backoff,
        checkpoint_path=checkpoint_path,
        checkpoint_every=checkpoint_every,
        progress_every=progress_every,
    )


def check_inchikeys_unichem_bulk_sync(inchikeys: list, **kwargs) -> Dict[str, bool]:
    """Blocking wrapper around :func:`check_inchikeys_unichem_bulk`.

    Call this from ordinary synchronous code (a script, a kernel cell, a
    check function). Accepts the same keyword arguments.
    """
    return run_async(check_inchikeys_unichem_bulk(inchikeys, **kwargs))


if __name__ == "__main__":
    example_keys = [
        "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",  # aspirin - valid
        "AAAAAAAAAAAAAA-AAAAAAAAAA-A",  # invalid
    ]
    print(check_inchikeys_unichem(example_keys))
    #print(await check_inchikeys_unichem_bulk(example_keys, concurrency=5, rate_per_sec=5.0))
