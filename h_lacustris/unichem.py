"""
Utility for checking whether InChIKeys resolve to valid compounds in UniChem.

UniChem Compound Search API docs:
https://chembl.gitbook.io/unichem/api/compound-search

Two entry points are provided:

- check_inchikeys_unichem(...)        - simple synchronous version (requests),
                                         fine for a handful to a few hundred keys.
- check_inchikeys_unichem_bulk(...)   - async version (aiohttp) with bounded
                                         concurrency, token-bucket rate limiting,
                                         progress logging, and crash-safe
                                         checkpointing to disk. Use this for
                                         lists of thousands of InChIKeys.
"""

import asyncio
import json
import os
import time
from typing import Dict, List, Optional

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
    async with sem:
        await limiter.acquire()
        payload = {"type": "inchikey", "compound": key}
        for attempt in range(max_retries):
            try:
                async with session.post(
                    url, json=payload, timeout=timeout
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        return key, bool(data.get("compounds"))
                    elif resp.status == 429:
                        # Rate-limited by the server: back off, honoring
                        # Retry-After if UniChem sends it.
                        retry_after = resp.headers.get("Retry-After")
                        wait = float(retry_after) if retry_after else retry_backoff * (2 ** attempt)
                        await asyncio.sleep(wait)
                        continue
                    elif 500 <= resp.status < 600:
                        await asyncio.sleep(retry_backoff * (2 ** attempt))
                        continue
                    else:
                        return key, False
            except (asyncio.TimeoutError, Exception) as e:
                # aiohttp.ClientError subclasses Exception; keep this dependency-light
                if attempt < max_retries - 1:
                    await asyncio.sleep(retry_backoff * (2 ** attempt))
                    continue
                print(f"[unichem] {key}: failed after {max_retries} attempts ({e}); marking as False")
                return key, False
        return key, False


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
        return [{"inchikey": k, "valid_inchikey": results[k]} for k in inchikeys]

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

    return [{"inchikey": k, "valid_inchikey": results[k]} for k in inchikeys]


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
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None

    if loop and loop.is_running():
        task = loop.create_task(
            _check_inchikeys_unichem_bulk_async(
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
        )
        return await task

    return asyncio.run(
        _check_inchikeys_unichem_bulk_async(
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
    )


if __name__ == "__main__":
    example_keys = [
        "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",  # aspirin - valid
        "AAAAAAAAAAAAAA-AAAAAAAAAA-A",  # invalid
    ]
    print(check_inchikeys_unichem(example_keys))
    #print(await check_inchikeys_unichem_bulk(example_keys, concurrency=5, rate_per_sec=5.0))
