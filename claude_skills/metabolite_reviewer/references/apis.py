"""Reference patterns for propagating metabolite annotations through public
chemistry APIs: UniChem, ChEBI and PubChem.

Why this module exists
----------------------
A model's name / formula / charge / xref block is usually *derived* from the
MetaNetX id it was assigned, so comparing the block against that same MNXM
cannot detect a wrong assignment. These three APIs are independent of that
assignment and of each other, which makes them the evidence of record for:

  * is this InChIKey a real compound?            -> `unichem_validate_inchikeys`
  * what else is this compound called elsewhere? -> `unichem_xrefs`
  * which charge states does this species have?  -> `chebi_charge_states`
  * what formula/charge does the database say?   -> `chebi_compound`, `pubchem_by_*`

Curation rules encoded here
---------------------------
1. Database annotation beats a calculated correction. If a mass/charge
   imbalance implies a charge that no database records, the imbalance is in
   the reaction or in another metabolite — not here.
2. A compound may legitimately have several charge states (the neutral acid
   and its conjugate bases). `chebi_charge_states` walks the ChEBI conjugate
   acid/base relations to enumerate them, so a proposed charge can be checked
   for membership rather than equality.
3. Propagate through identifiers before guessing from names. An id-based
   lookup is evidence; a name search is a hint.

Rate limiting
-------------
Every call goes through `_Throttle`, a process-wide token bucket per host.
Defaults follow each service's published guidance:

  * PubChem PUG-REST: no more than 5 requests/second (documented limit).
  * EBI (UniChem, ChEBI): no published limit; 8 req/s is used here as a
    courteous default for a shared public service.

Responses are memoised on disk (`cache_dir`) so a re-run of a review costs no
requests at all. Delete the cache directory to force a refresh.

Networks
--------
Requires outbound access to `www.ebi.ac.uk` and `pubchem.ncbi.nlm.nih.gov`.
Every function degrades to `None` / empty rather than raising, so a review can
still run offline — but `annotation_evidence` then reports `source="offline"`
and the caller must not treat absence of evidence as evidence of an error.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import requests

UNICHEM = "https://www.ebi.ac.uk/unichem/api/v1"
CHEBI = "https://www.ebi.ac.uk/chebi/backend/api/public"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

DEFAULT_CACHE = Path("api_cache")

# requests/second per host — see module docstring
RATE_LIMITS = {
    "www.ebi.ac.uk": 8.0,
    "pubchem.ncbi.nlm.nih.gov": 5.0,
}


class _Throttle:
    """Process-wide token bucket, one bucket per host."""

    _locks: dict[str, threading.Lock] = {}
    _next: dict[str, float] = {}
    _guard = threading.Lock()

    @classmethod
    def wait(cls, host: str) -> None:
        rate = RATE_LIMITS.get(host, 5.0)
        interval = 1.0 / rate
        with cls._guard:
            lock = cls._locks.setdefault(host, threading.Lock())
        with lock:
            now = time.monotonic()
            nxt = cls._next.get(host, 0.0)
            if nxt > now:
                time.sleep(nxt - now)
                now = time.monotonic()
            cls._next[host] = max(now, nxt) + interval


@dataclass
class ApiClient:
    """Cached, rate-limited HTTP client for the three chemistry services.

    Parameters
    ----------
    cache_dir : directory for the on-disk response cache; None disables it.
    timeout   : per-request timeout in seconds.
    retries   : attempts on timeout / 5xx / 429 (429 honours `Retry-After`).
    offline   : never touch the network; serve from cache only. Useful for a
                reproducible re-run of a review.
    """

    cache_dir: Path | None = DEFAULT_CACHE
    timeout: int = 15
    retries: int = 4
    backoff: float = 1.5
    offline: bool = False
    _local: threading.local = field(default_factory=threading.local, repr=False)
    stats: dict = field(default_factory=lambda: {"hit": 0, "miss": 0, "fail": 0},
                        repr=False)

    def __post_init__(self):
        if self.cache_dir is not None:
            self.cache_dir = Path(self.cache_dir)
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    # -- transport session ------------------------------------------------
    @property
    def session(self) -> requests.Session:
        """One Session per thread.

        Bulk work in this module is async, not threaded (see
        ``unichem.py``), so this is defensive rather than load-bearing: it
        costs nothing and means a caller who does wrap the synchronous client
        in a thread cannot corrupt the connection pool, which
        ``requests.Session`` does not guard against.
        """
        s = getattr(self._local, "session", None)
        if s is None:
            s = requests.Session()
            self._local.session = s
        return s

    def _reset_session(self) -> None:
        s = getattr(self._local, "session", None)
        if s is not None:
            s.close()
        self._local.session = requests.Session()

    # -- cache ------------------------------------------------------------
    def _cache_file(self, key: str) -> Path | None:
        if self.cache_dir is None:
            return None
        h = hashlib.sha1(key.encode()).hexdigest()[:24]
        return self.cache_dir / f"{h}.json"

    def _read_cache(self, key: str):
        p = self._cache_file(key)
        if p and p.exists():
            try:
                return json.loads(p.read_text())
            except (json.JSONDecodeError, OSError):
                return None
        return None

    def _write_cache(self, key: str, value) -> None:
        p = self._cache_file(key)
        if p is None:
            return
        tmp = p.with_suffix(".tmp")
        tmp.write_text(json.dumps(value))
        os.replace(tmp, p)          # atomic; never leaves a truncated file

    def cache_key(self, method: str, url: str, *, json_body=None,
                  params=None) -> str:
        """The cache key for a request.

        Public so the async bulk helpers can consult and populate the same
        cache the synchronous path uses — a review re-run must not re-issue
        requests just because the first run went out concurrently.
        """
        return (f"{method} {url} {json.dumps(json_body, sort_keys=True)} "
                f"{json.dumps(params, sort_keys=True)}")

    def cache_put(self, method: str, url: str, json_body, params, data) -> None:
        """Store a response fetched outside :meth:`request`."""
        self._write_cache(self.cache_key(method, url, json_body=json_body,
                                         params=params), {"data": data})

    # -- transport --------------------------------------------------------
    def request(self, method: str, url: str, *, json_body=None, params=None):
        key = self.cache_key(method, url, json_body=json_body, params=params)
        cached = self._read_cache(key)
        if cached is not None:
            self.stats["hit"] += 1
            return cached.get("data")
        if self.offline:
            return None

        host = url.split("/")[2]
        for attempt in range(self.retries):
            _Throttle.wait(host)
            try:
                resp = self.session.request(
                    method, url, json=json_body, params=params,
                    timeout=self.timeout,
                    headers={"Accept": "application/json"})
            except requests.RequestException:
                # A pooled keep-alive connection that has gone stale behind a
                # proxy shows up as a read timeout on the *first* request and
                # succeeds immediately afterwards. Drop the pool and retry.
                self._reset_session()
                time.sleep(self.backoff * (2 ** attempt))
                continue
            if resp.status_code == 200:
                try:
                    data = resp.json()
                except ValueError:
                    data = None
                self.stats["miss"] += 1
                self._write_cache(key, {"data": data})
                return data
            if resp.status_code == 404:
                # a definite "no such record" — cache it, it will not change
                self.stats["miss"] += 1
                self._write_cache(key, {"data": None})
                return None
            if resp.status_code == 429:
                wait = float(resp.headers.get("Retry-After",
                                              self.backoff * (2 ** attempt)))
                time.sleep(wait)
                continue
            if 500 <= resp.status_code < 600:
                time.sleep(self.backoff * (2 ** attempt))
                continue
            break
        self.stats["fail"] += 1
        return None

    def get(self, url, **kw):
        return self.request("GET", url, **kw)

    def post(self, url, **kw):
        return self.request("POST", url, **kw)


# ==========================================================================
# UniChem
# ==========================================================================

def unichem_compound(client: ApiClient, inchikey: str) -> dict | None:
    """Raw UniChem record for an InChIKey, or None when it resolves nowhere.

    UniChem answers 200 with an empty `compounds` list for a malformed or
    unknown key, so validity is decided by the payload, not the status code.
    """
    data = client.post(f"{UNICHEM}/compounds",
                       json_body={"type": "inchikey", "compound": inchikey})
    if not data or not data.get("compounds"):
        return None
    return data["compounds"][0]


def unichem_validate_inchikeys(client: ApiClient, inchikeys: Iterable[str],
                               concurrency: int = 32) -> dict[str, bool]:
    """{inchikey: resolves_to_a_real_compound}. Duplicates queried once.

    UniChem has no bulk endpoint for this, so validating a whole model is one
    POST per distinct key — around 800 for a mid-sized reconstruction. That
    endpoint answers in ~13 s through a proxy, so serially this single check
    is three hours and dominates everything else the review does.

    The requests are independent and are pure network wait, so they go out
    concurrently on an event loop (``references.unichem.bulk_fetch_json``).
    Measured on this model: 0.56 s/key at concurrency 20, against ~13 s/key
    serially. Results are written into the same on-disk cache the synchronous
    client uses, so a second run of the review costs nothing.
    """
    from unichem import bulk_fetch_json, run_async     # noqa: PLC0415

    keys = list(dict.fromkeys(k for k in inchikeys if k))
    if not keys:
        return {}

    # Serve whatever the shared cache already holds; only miss goes to the wire.
    out: dict[str, bool] = {}
    todo: list[dict] = []
    url = f"{UNICHEM}/compounds"
    for key in keys:
        body = {"type": "inchikey", "compound": key}
        ck = client.cache_key("POST", url, json_body=body, params=None)
        cached = client._read_cache(ck)
        if cached is not None:
            client.stats["hit"] += 1
            out[key] = bool((cached.get("data") or {}).get("compounds"))
        elif client.offline:
            out[key] = False
        else:
            todo.append({"key": key, "url": url, "method": "POST", "json": body})

    if todo:
        fetched = run_async(bulk_fetch_json(
            todo, concurrency=concurrency,
            rate_per_sec=RATE_LIMITS.get("www.ebi.ac.uk", 8.0),
            timeout=client.timeout, max_retries=client.retries,
            retry_backoff=client.backoff,
            progress_every=max(50, len(todo) // 10), label="unichem"))
        for spec in todo:
            key = spec["key"]
            data = fetched.get(key)
            out[key] = bool((data or {}).get("compounds"))
            if data is not None:
                client.stats["miss"] += 1
                client.cache_put("POST", url, spec["json"], None, data)
            else:
                client.stats["fail"] += 1   # not cached: may be transient
    return {k: out.get(k, False) for k in keys}


def unichem_xrefs(client: ApiClient, inchikey: str) -> dict[str, list[str]]:
    """{source_short_name: [ids]} for every database UniChem links this key to.

    This is the cheapest way to turn one structural key into ChEBI, KEGG,
    HMDB, LIPID MAPS, MetaCyc and PubChem ids in a single request.
    """
    comp = unichem_compound(client, inchikey)
    if not comp:
        return {}
    out: dict[str, list[str]] = {}
    for src in comp.get("sources", []):
        name = (src.get("shortName") or src.get("name") or "").lower()
        cid = src.get("compoundId")
        if name and cid:
            out.setdefault(name, []).append(str(cid))
    return out


def unichem_connectivity(client: ApiClient, inchikey: str) -> dict[str, list[str]]:
    """Xrefs for every compound sharing this key's *connectivity layer*.

    Use when the exact key misses: it finds the same skeleton in other
    protonation / stereo states, which is precisely the situation a charge
    review is in. Falls back to the exact-match xrefs if the endpoint fails.
    """
    data = client.post(f"{UNICHEM}/connectivity",
                       json_body={"type": "inchikey", "compound": inchikey,
                                  "searchComponents": False})
    if not data:
        return unichem_xrefs(client, inchikey)
    out: dict[str, list[str]] = {}
    for group in data.get("searchedCompound", []) or data.get("compounds", []):
        for src in (group.get("sources") or []):
            name = (src.get("shortName") or "").lower()
            cid = src.get("compoundId")
            if name and cid:
                out.setdefault(name, []).append(str(cid))
    return out


# ==========================================================================
# ChEBI
# ==========================================================================

# --------------------------------------------------------------------------
# Request specs
# --------------------------------------------------------------------------
# Each `spec_*` returns the exact request its sibling endpoint function will
# make, so `warm_cache` produces cache keys that the synchronous path then
# hits. Defining them beside the endpoint is deliberate: a URL that appears
# twice in this module will drift, and a drifted prefetch fails silently —
# every request simply goes out again serially, which is the bug the prefetch
# exists to prevent.

def spec_unichem_compound(inchikey: str) -> dict:
    return {"url": f"{UNICHEM}/compounds", "method": "POST",
            "json": {"type": "inchikey", "compound": inchikey}}


def spec_chebi_compound(chebi_id) -> dict | None:
    num = _chebi_num(chebi_id)
    return {"url": f"{CHEBI}/compound/{num}/"} if num else None


def spec_chebi_search(term: str, limit: int = 5) -> dict:
    return {"url": f"{CHEBI}/es_search/", "params": {"term": str(term),
                                                     "size": limit}}


def spec_pubchem_by_name(name: str) -> dict:
    from urllib.parse import quote
    return {"url": f"{PUBCHEM}/compound/name/{quote(str(name))}"
                   f"/property/{PUBCHEM_PROPS}/JSON"}


def spec_pubchem_by_inchikey(inchikey: str) -> dict:
    return {"url": f"{PUBCHEM}/compound/inchikey/{inchikey}"
                   f"/property/{PUBCHEM_PROPS}/JSON"}


def _chebi_num(chebi_id) -> str | None:
    m = re.search(r"(\d+)", str(chebi_id or ""))
    return m.group(1) if m else None


def chebi_compound(client: ApiClient, chebi_id) -> dict | None:
    """Normalised ChEBI record.

    Returns name / formula / charge / mass / inchikey / smiles / xrefs plus the
    conjugate acid-base neighbours, which is what makes charge review possible.
    """
    num = _chebi_num(chebi_id)
    if not num:
        return None
    raw = client.get(spec_chebi_compound(num)["url"])
    if not raw:
        return None
    chem = raw.get("chemical_data") or {}
    struct = raw.get("default_structure") or {}
    xrefs: dict[str, list[str]] = {}
    for group in (raw.get("database_accessions") or {}).values():
        for x in group:
            prefix = (x.get("prefix") or x.get("source_name") or "").lower()
            acc = x.get("accession_number")
            if prefix and acc:
                xrefs.setdefault(prefix, []).append(str(acc))
    rel = raw.get("ontology_relations") or {}
    conjugates = []
    for direction in ("incoming_relations", "outgoing_relations"):
        for x in rel.get(direction, []):
            rt = (x.get("relation_type") or "").lower()
            if "conjugate" not in rt:
                continue
            other = (x["init_id"] if direction == "incoming_relations"
                     else x["final_id"])
            other_name = (x["init_name"] if direction == "incoming_relations"
                          else x["final_name"])
            conjugates.append({"chebi": f"CHEBI:{other}", "name": other_name,
                               "relation": rt})
    return {
        "chebi": f"CHEBI:{num}",
        "name": raw.get("name"),
        "formula": chem.get("formula"),
        "charge": chem.get("charge"),
        "mass": chem.get("mass"),
        "inchikey": struct.get("standard_inchi_key") or struct.get("inchikey"),
        "smiles": struct.get("smiles"),
        "stars": raw.get("stars"),
        "xrefs": xrefs,
        "conjugates": conjugates,
    }


def chebi_charge_states(client: ApiClient, chebi_id,
                        max_depth: int = 2) -> list[dict]:
    """Every charge state of a species, by walking conjugate acid/base links.

    A metabolite is not restricted to one charge: malonyl-CoA exists in ChEBI
    as the neutral acid and as the 5- anion, and *only* those. A proposed
    charge should be checked for membership in this set. If the charge implied
    by a mass/charge imbalance is absent from it, the imbalance belongs to the
    reaction or to another metabolite — do not "fix" it here.

    Returns [{chebi, name, formula, charge}, ...] including the seed entry.
    """
    seed = chebi_compound(client, chebi_id)
    if not seed:
        return []
    seen = {seed["chebi"]}
    out = [{k: seed[k] for k in ("chebi", "name", "formula", "charge")}]
    frontier = [(seed, 0)]
    while frontier:
        node, depth = frontier.pop()
        if depth >= max_depth:
            continue
        for conj in node.get("conjugates", []):
            if conj["chebi"] in seen:
                continue
            seen.add(conj["chebi"])
            rec = chebi_compound(client, conj["chebi"])
            if not rec:
                continue
            out.append({k: rec[k] for k in ("chebi", "name", "formula", "charge")})
            frontier.append((rec, depth + 1))
    return out


def chebi_search(client: ApiClient, name: str, limit: int = 5) -> list[dict]:
    """Name search. A hint, not evidence — always confirm with a structure."""
    spec = spec_chebi_search(name, limit)
    data = client.get(spec["url"], params=spec["params"])
    if not data:
        return []
    hits = data.get("results") or data.get("hits") or []
    out = []
    for h in hits[:limit]:
        src = h.get("_source", h)
        out.append({"chebi": src.get("chebi_accession") or src.get("id"),
                    "name": src.get("name"), "score": h.get("_score")})
    return out


# ==========================================================================
# PubChem
# ==========================================================================

PUBCHEM_PROPS = ("MolecularFormula,Charge,IUPACName,InChIKey,"
                 "CanonicalSMILES,MonoisotopicMass")


def _pubchem_unwrap(data) -> dict | None:
    if not data:
        return None
    props = (data.get("PropertyTable") or {}).get("Properties") or []
    return props[0] if props else None


def _pubchem_data(client: ApiClient, spec: dict) -> dict | None:
    return _pubchem_unwrap(client.get(spec["url"]))


def _pubchem_props(client: ApiClient, route: str) -> dict | None:
    return _pubchem_unwrap(
        client.get(f"{PUBCHEM}/compound/{route}/property/{PUBCHEM_PROPS}/JSON"))


def pubchem_by_inchikey(client: ApiClient, inchikey: str) -> dict | None:
    return _pubchem_data(client, spec_pubchem_by_inchikey(inchikey))


def pubchem_by_cid(client: ApiClient, cid) -> dict | None:
    return _pubchem_props(client, f"cid/{cid}")


def pubchem_by_name(client: ApiClient, name: str) -> dict | None:
    """Name lookup. PubChem resolves synonyms generously — treat a hit as a
    hint and confirm the InChIKey before propagating anything from it."""
    return _pubchem_data(client, spec_pubchem_by_name(name))


def pubchem_xrefs(client: ApiClient, cid) -> dict[str, list[str]]:
    data = client.get(f"{PUBCHEM}/compound/cid/{cid}/xrefs/RegistryID,RN/JSON")
    if not data:
        return {}
    info = (data.get("InformationList") or {}).get("Information") or []
    return {k.lower(): [str(x) for x in v]
            for k, v in (info[0] if info else {}).items() if k != "CID"}


# ==========================================================================
# the thing a reviewer actually calls
# ==========================================================================

# Priority when sources disagree. ChEBI is first because it is manually
# curated and models charge states explicitly; PubChem is broader but
# auto-generated; UniChem contributes identity links rather than properties.
SOURCE_PRIORITY = ("chebi", "pubchem", "unichem")


def _same_name(a, b) -> bool:
    """Do two compound names denote the same thing?

    Strips the HTML markup ChEBI returns in search hits ("<em>a</em>"), then
    all punctuation and case. Deliberately strict — this decides whether a
    name-search hit is allowed to supply formula and charge, so a near miss
    must fail. ``plastoquinol`` matches ``Plastoquinol``; it does not match
    ``plastoquinol-9``.
    """
    def norm(s):
        s = re.sub(r"<[^>]+>", "", str(s or ""))
        return re.sub(r"[^a-z0-9]", "", s.lower())
    na, nb = norm(a), norm(b)
    return bool(na) and na == nb


def annotation_evidence(client: ApiClient, *, inchikey: str | None = None,
                        chebi_id: str | None = None,
                        name: str | None = None,
                        pubchem_cid: str | None = None) -> dict:
    """Gather independent formula/charge evidence for one metabolite.

    Propagation order follows the curation rule "identifiers before names":
    a supplied ChEBI id is used directly; otherwise the InChIKey is resolved
    through UniChem to a ChEBI id and a PubChem CID; a name search is the last
    resort and is reported as such.

    Returns
    -------
    dict with
      formulas      {source: formula}
      charges       {source: charge}
      charge_states [{chebi, name, formula, charge}, ...] — the legitimate
                    protonation family, empty when ChEBI could not be reached
      xrefs         {namespace: [ids]} harvested along the way
      sources       which services actually answered
      resolved_via  'chebi_id' | 'inchikey' | 'name' | 'offline'
    """
    ev = {"formulas": {}, "charges": {}, "charge_states": [], "xrefs": {},
          "sources": [], "resolved_via": None}

    if inchikey and not chebi_id:
        xr = unichem_xrefs(client, inchikey)
        if xr:
            ev["xrefs"].update(xr)
            ev["sources"].append("unichem")
            ev["resolved_via"] = "inchikey"
            for key in ("chebi", "chebi_id"):
                if xr.get(key):
                    chebi_id = xr[key][0]
                    break
            if not pubchem_cid and xr.get("pubchem"):
                pubchem_cid = xr["pubchem"][0]

    if chebi_id:
        rec = chebi_compound(client, chebi_id)
        if rec:
            ev["sources"].append("chebi")
            ev["resolved_via"] = ev["resolved_via"] or "chebi_id"
            if rec["formula"]:
                ev["formulas"]["chebi"] = rec["formula"]
            if rec["charge"] is not None:
                ev["charges"]["chebi"] = int(rec["charge"])
            for k, v in rec["xrefs"].items():
                ev["xrefs"].setdefault(k, []).extend(v)
            ev["charge_states"] = chebi_charge_states(client, chebi_id)

    if not chebi_id and name:
        # ChEBI's own name index resolves the biological vocabulary that
        # PubChem does not carry ("plastoquinol", "divinylprotochlorophyllide
        # a"). Accept only an exact match on the normalised name: a loose hit
        # would silently attach the wrong structure, which is worse than no
        # evidence at all. Markup and stereo-descriptor differences are
        # normalised away; anything else is rejected.
        for cand in chebi_search(client, name, limit=5):
            if _same_name(cand.get("name"), name):
                chebi_id = cand.get("chebi") or cand.get("id")
                ev["resolved_via"] = ev["resolved_via"] or "name"
                break
        if chebi_id:
            rec = chebi_compound(client, chebi_id)
            if rec:
                ev["sources"].append("chebi")
                if rec["formula"]:
                    ev["formulas"]["chebi"] = rec["formula"]
                if rec["charge"] is not None:
                    ev["charges"]["chebi"] = int(rec["charge"])
                for k, v in rec["xrefs"].items():
                    ev["xrefs"].setdefault(k, []).extend(v)
                ev["charge_states"] = chebi_charge_states(client, chebi_id)

    if pubchem_cid or inchikey or name:
        rec = None
        if pubchem_cid:
            rec = pubchem_by_cid(client, pubchem_cid)
        if rec is None and inchikey:
            rec = pubchem_by_inchikey(client, inchikey)
        if rec is None and name:
            rec = pubchem_by_name(client, name)
            if rec:
                ev["resolved_via"] = ev["resolved_via"] or "name"
        if rec:
            ev["sources"].append("pubchem")
            if rec.get("MolecularFormula"):
                ev["formulas"]["pubchem"] = rec["MolecularFormula"]
            if rec.get("Charge") is not None:
                ev["charges"]["pubchem"] = int(rec["Charge"])
            if rec.get("CID"):
                ev["xrefs"].setdefault("pubchem", []).append(str(rec["CID"]))

    if not ev["sources"]:
        ev["resolved_via"] = "offline" if client.offline else "unresolved"
    ev["xrefs"] = {k: sorted(set(v)) for k, v in ev["xrefs"].items()}
    return ev


def warm_cache(client: ApiClient, specs: list[dict], *,
               concurrency: int = 24, label: str = "warm") -> int:
    """Fetch a batch of independent requests concurrently into the cache.

    Returns the number of responses newly cached. Requests already cached are
    skipped, so this is cheap to call repeatedly and free on a re-run.

    This is the primitive the prefetch waves are built from: issue everything
    that can go out at once, then let the ordinary synchronous code paths run
    against a warm cache. Keeping it cache-shaped means there is exactly one
    implementation of the parsing, retry and priority logic — the async layer
    only moves bytes.
    """
    from unichem import bulk_fetch_json, run_async     # noqa: PLC0415

    if client.offline:
        return 0
    todo, seen = [], set()
    for spec in specs:
        method, url = spec.get("method", "GET"), spec["url"]
        body, params = spec.get("json"), spec.get("params")
        ck = client.cache_key(method, url, json_body=body, params=params)
        if ck in seen or client._read_cache(ck) is not None:
            continue
        seen.add(ck)
        todo.append({"key": ck, "url": url, "method": method,
                     "json": body, "params": params,
                     "_body": body, "_params": params})
    if not todo:
        return 0

    host = todo[0]["url"].split("/")[2]
    fetched = run_async(bulk_fetch_json(
        todo, concurrency=concurrency,
        rate_per_sec=RATE_LIMITS.get(host, 5.0), timeout=client.timeout,
        max_retries=client.retries, retry_backoff=client.backoff,
        progress_every=max(50, len(todo) // 10), label=label))

    n = 0
    for spec in todo:
        data = fetched.get(spec["key"])
        if data is not None:
            client.cache_put(spec.get("method", "GET"), spec["url"],
                             spec["_body"], spec["_params"], data)
            client.stats["miss"] += 1
            n += 1
    return n


def prefetch_annotation_evidence(client: ApiClient, rows: list[dict], *,
                                 id_key: str = "id",
                                 concurrency: int = 24) -> dict[str, dict]:
    """Gather evidence for many metabolites at once. {id: evidence}.

    ``annotation_evidence`` for one metabolite is a chain of dependent calls:
    resolve a name to a ChEBI id, read that compound, then walk its conjugate
    family. The chain cannot be collapsed — but each *stage* is independent
    across metabolites, so the whole batch advances one stage at a time with
    every request in a stage issued concurrently:

        wave 1  name search + PubChem name lookup + UniChem key lookup
        wave 2  ChEBI compound records for the ids wave 1 resolved
        wave 3  conjugate acid/base relations for those same ids

    Three concurrent round-trips replace 3N serial ones. Each wave only warms
    the shared cache; ``annotation_evidence`` is then run normally per
    metabolite and finds every response already local, so the parsing and
    source-priority logic exists in exactly one place.

    Async rather than threads throughout — see the note in ``unichem.py``.
    """
    todo = [r for r in rows if r.get(id_key)]
    if not todo or client.offline:
        return {r[id_key]: annotation_evidence(
                    client, name=r.get("name"), inchikey=r.get("inchikey"),
                    chebi=r.get("chebi"), formula=r.get("formula"))
                for r in todo}

    # ---- wave 1: everything reachable from the row as given ---------------
    ebi1, pc1 = [], []
    for r in todo:
        name, ik, ch = r.get("name"), r.get("inchikey"), r.get("chebi")
        if ik:
            ebi1.append(spec_unichem_compound(ik))
            pc1.append(spec_pubchem_by_inchikey(ik))
        if ch:
            s = spec_chebi_compound(ch)
            if s:
                ebi1.append(s)
        elif name:
            ebi1.append(spec_chebi_search(name))
        if name and not ik:
            pc1.append(spec_pubchem_by_name(name))
    warm_cache(client, ebi1, concurrency=concurrency, label="evidence/ebi")
    warm_cache(client, pc1, concurrency=8, label="evidence/pubchem")

    # ---- wave 2: ChEBI ids are known only once wave 1 resolved the names --
    # (chebi_search now reads from cache, so this loop makes no requests)
    chebi_ids = set()
    for r in todo:
        cid = r.get("chebi")
        if not cid and r.get("name"):
            for cand in chebi_search(client, r["name"], limit=5):
                if _same_name(cand.get("name"), r["name"]):
                    cid = cand.get("chebi") or cand.get("id")
                    break
        if cid:
            chebi_ids.add(cid)
    warm_cache(client, [s for s in (spec_chebi_compound(c) for c in chebi_ids)
                        if s], concurrency=concurrency, label="evidence/chebi")

    # ---- wave 3: conjugate acid/base neighbours of those compounds --------
    # chebi_charge_states walks one hop out; those records are independent of
    # each other, so fetch them all before the walk runs serially.
    neighbours = set()
    for cid in chebi_ids:
        rec = chebi_compound(client, cid)               # cached by wave 2
        for c in (rec or {}).get("conjugates", []):
            neighbours.add(c["chebi"])
    warm_cache(client, [s for s in (spec_chebi_compound(c) for c in neighbours)
                        if s], concurrency=concurrency, label="evidence/conj")

    # ---- assemble from the warm cache (no network left to do) -------------
    return {r[id_key]: annotation_evidence(
                client, name=r.get("name"), inchikey=r.get("inchikey"),
                chebi=r.get("chebi"), formula=r.get("formula"))
            for r in todo}


def consensus_charge(evidence: dict) -> tuple[int | None, str]:
    """Single best charge from gathered evidence, plus the source that gave it.

    Follows SOURCE_PRIORITY rather than voting: a manually curated value is
    worth more than agreement between two auto-generated ones.
    """
    for src in SOURCE_PRIORITY:
        if src in evidence.get("charges", {}):
            return evidence["charges"][src], src
    return None, "none"


def charge_is_supported(evidence: dict, charge: int | None) -> bool:
    """Is `charge` a state any database records for this species?

    True when it matches a reported charge or any member of the ChEBI
    conjugate family. A False here on a balance-derived correction is the
    signal to defer the fix to a reaction review instead of applying it.
    """
    if charge is None:
        return False
    if charge in evidence.get("charges", {}).values():
        return True
    return any(cs.get("charge") == charge
               for cs in evidence.get("charge_states", []))


def is_probable_fantasy(evidence: dict, *, model_xrefs: dict | None = None) -> bool:
    """No structural database knows this compound.

    True when nothing resolved through UniChem, ChEBI or PubChem and the model
    carries no structural xref of its own — only namespace ids (BiGG, SEED,
    MetaNetX) that a reconstruction pipeline can mint without evidence.
    """
    if evidence.get("sources"):
        return False
    structural = {"chebi", "inchikey", "smiles", "hmdb", "lipidmaps",
                  "kegg.compound", "metacyc.compound", "pubchem"}
    have = {k for k, v in (model_xrefs or {}).items() if v}
    return not (have & structural)


if __name__ == "__main__":
    c = ApiClient(cache_dir=None)
    print(unichem_validate_inchikeys(c, ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                                         "AAAAAAAAAAAAAA-AAAAAAAAAA-A"]))
    print(chebi_charge_states(c, "CHEBI:15531"))      # malonyl-CoA: 0 and -5
    print(pubchem_by_name(c, "acetyl-CoA"))
