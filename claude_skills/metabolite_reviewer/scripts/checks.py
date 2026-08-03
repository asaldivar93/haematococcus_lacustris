"""Deterministic metabolite checks for GSMM curation (BiGG namespace).

Design principle
----------------
Model annotation blocks (name / formula / charge / inchikey / all xrefs) are
almost always *derived* from the assigned MetaNetX id. Comparing them against
that same MNXM therefore cannot detect a wrong MNXM assignment — the block is
internally consistent by construction. Every check below is built to draw on a
source that is NOT downstream of the assigned MNXM:

  * reaction stoichiometry  (mass / charge balance, delta localisation)
  * the BiGG identifier and BiGG's own InChIKey
  * cross-compartment agreement of the same base id
  * duplicate structure within a compartment
  * the live chemistry APIs (UniChem / ChEBI / PubChem) in references/apis.py

Curation policy (set by the curator, round 2)
---------------------------------------------
1.  Database annotation outranks arithmetic. A formula or charge implied by
    balance is a *hypothesis*, never a fix. Only a value carried by a
    structural database may be written into the model.
2.  A metabolite may legitimately hold several charge states. Enumerate them
    (ChEBI conjugate acid/base relations) before calling a charge wrong.
3.  If no database value resolves an imbalance, the error is not in this
    metabolite: raise a reaction-reviewer flag instead of proposing a fix.
4.  Balance findings are held out of round 1 entirely — they go to a separate
    file and are only actioned in round 2, after annotations are settled.
5.  MetaNetX identifier problems are not errors; they are annotation work.
    They go to the annotation-update file for human review, not the worklist.
6.  A metabolite carrying any xref beyond the three namespaces a
    reconstruction pipeline can mint (bigg / metanetx / seed) is annotated
    enough — do not report missing annotations for it.
7.  Two species that are chemically identical but hold *distinct* ids in an
    independent namespace are distinct entities, not duplicates.
8.  Where two BiGG ids are both valid, prefer the one used by more models.
9.  An InChIKey must be validated against UniChem before being flagged.
10. Every finding that a metabolite-level fix cannot resolve carries a note
    for the downstream reaction_reviewer agent.

Usage
-----
    from checks import ReviewContext, run_all
    ctx = ReviewContext(model=..., mets=..., mnx=..., depr=..., bigg=...)
    findings, extras = run_all(ctx)
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from fractions import Fraction

import polars as pl

# --------------------------------------------------------------------------
# conventions
# --------------------------------------------------------------------------

#: Elements that are placeholders, not real atoms. Reactions containing a
#: metabolite with any of these cannot be mass-balanced and must be excluded
#: from balance-based checks rather than reported as errors.
PLACEHOLDER_ELEMENTS = {"R", "X", "Z", "*"}

#: Proton species are conventionally written with charge +1 and formula H in
#: BiGG/MetaNetX. Models carrying charge 0 produce a flood of spurious charge
#: imbalance; normalise before balancing and report the convention once.
PROTON_NAMES = {"h+", "h", "proton", "hydrogen ion"}

#: MetaNetX xref values for BiGG are lowercased and may carry an `m_`
#: prefix; BiGG ids use `__l`-style stereo suffixes where older releases
#: used `-l`.
_STEREO = [("-l", "__l"), ("-d", "__d"), ("-r", "__r"), ("-s", "__s")]

#: Names generated from a template rather than issued by a nomenclature body.
#: No database indexes every acyl permutation, every photon wavelength bin or
#: every carrier-protein conjugate, so a failed name search on one of these
#: says nothing about whether the compound is real. Without this guard the
#: fantasy check flags 63 legitimate species (plastocyanin, photons, ACP
#: conjugates, every phosphatidate) against 2 real hits.
SYSTEMATIC_NAME = re.compile(
    r"\d{1,2}:\d|\(\d+[ZE]|sn-glycerol|\[acyl-carrier|\bACP\b|^Photon|"
    r"^Protein |repeat units|Dummy|GlcNAc|n-C\d|Plastocyanin|flavoprotein|"
    r"octadec|hexadec|tetradec|dodec|geranylgeranyl", re.I)

ERROR_TAGS = (
    "deprecated_metanetx_id",
    "unknown_metanetx_id",
    "missing_metanetx_id",
    "multiple_metanetx_id",
    "formula_mismatch_reference",
    "charge_mismatch_reference",
    "inchikey_mismatch_reference",
    "invalid_inchikey",
    "inchikey_skeleton_conflict_bigg",
    "xref_disagreement",
    "bigg_id_not_in_namespace",
    "bigg_id_deprecated",
    "duplicate_metabolite",
    "chemically_identical_but_distinct_in_databases",
    "name_collision",
    "shadow_duplicate",
    "uncurated_name_stub",
    "uncurated_metabolite",
    "probable_fantasy_metabolite",
    "possible_fantasy_metabolite",
    "charge_implausible_for_formula",
    "reference_lacks_structure",
    "cross_compartment_inconsistency",
    "mass_imbalance_localised",
    "charge_imbalance_localised",
    "proton_charge_convention",
    "placeholder_formula",
    "orphan_metabolite",
    "dead_end_metabolite",
    "missing_annotations",
    "schema_anomaly",
)

CONFIDENCE = ("high", "medium", "low")


# --------------------------------------------------------------------------
# normalisation helpers
# --------------------------------------------------------------------------

def norm_formula(formula: str | None) -> str | None:
    """Canonical element-multiset string, so C6H12O6 == H12C6O6."""
    if not formula:
        return None
    counts: Counter[str] = Counter()
    for element, number in re.findall(r"([A-Z][a-z]?|\*)(\d*)", str(formula)):
        if not element:
            continue
        counts[element] += int(number) if number else 1
    if not counts:
        return None
    return "".join(f"{e}{counts[e]}" for e in sorted(counts))


def parse_formula(formula: str | None) -> dict[str, int]:
    if not formula:
        return {}
    counts: Counter[str] = Counter()
    for element, number in re.findall(r"([A-Z][a-z]?|\*)(\d*)", str(formula)):
        if not element:
            continue
        counts[element] += int(number) if number else 1
    return dict(counts)


def norm_name(name: str | None) -> str | None:
    if name is None:
        return None
    return re.sub(r"[^a-z0-9]", "", str(name).lower()) or None


def norm_bigg(value: str) -> set[str]:
    """All plausible spellings of one BiGG id (case, m_ prefix, stereo suffix)."""
    v = str(value).strip().lower()
    if v.startswith(("m_", "r_")):
        v = v[2:]
    out = {v}
    for old, new in _STEREO:
        if v.endswith(old):
            out.add(v[: -len(old)] + new)
        if v.endswith(new):
            out.add(v[: -len(new)] + old)
    return out


def split_multi(value) -> set[str]:
    """MetaNetX/BiGG multi-valued cells use `|` or `;`."""
    if value is None:
        return set()
    return {x.strip() for x in re.split(r"[|;]", str(value)) if x.strip()}


def inchikey_skeleton(key: str | None) -> str | None:
    """First InChIKey block — connectivity only, ignores stereo/protonation."""
    if not key:
        return None
    key = str(key).strip()
    return key.split("-")[0] if "-" in key else key[:14]


def is_proton(metabolite) -> bool:
    return (norm_formula(metabolite.formula) == "H1"
            and str(metabolite.name or "").strip().lower() in PROTON_NAMES)


def has_placeholder(metabolite) -> bool:
    return bool(set(parse_formula(metabolite.formula)) & PLACEHOLDER_ELEMENTS)


# --------------------------------------------------------------------------
# context
# --------------------------------------------------------------------------

@dataclass
class ReviewContext:
    """Everything the checks need: the model, its table, and the references."""

    model: object                      # cobra.Model
    mets: pl.DataFrame                 # metabolite table, one row per model met
    # Every reference is optional. A check whose reference is absent skips
    # rather than fails, so a review with no tables and no network still
    # produces the offline findings — `review` reports what it skipped.
    mnx: pl.DataFrame | None = None    # MetaNetX chem_prop + pivoted xrefs
    depr: pl.DataFrame | None = None   # MetaNetX chem_depr (old,new,version)
    bigg: pl.DataFrame | None = None   # BiGG universal metabolites
    bigg_rxns: pl.DataFrame | None = None   # BiGG universal reactions
    apis: object | None = None         # references.apis.ApiClient, or None
    id_col: str = "id"
    compartment_pattern: str = r"_([a-z][a-z0-9]?)$"
    xref_columns: tuple[str, ...] = (
        "kegg.compound", "chebi", "metacyc.compound", "seed.compound",
        "hmdb", "lipidmaps", "reactome", "rhea", "kegg.drug", "kegg.glycan",
    )
    #: Namespaces a reconstruction pipeline can mint by itself. A metabolite
    #: annotated *only* in these has not actually been curated against any
    #: structural database (policy 6).
    pipeline_namespaces: tuple[str, ...] = (
        "bigg.metabolite", "metanetx.chemical", "seed.compound",
    )
    #: Curated namespaces whose identifiers are assigned by humans. Two
    #: species holding *different* ids here are distinct entities (policy 7).
    independent_namespaces: tuple[str, ...] = (
        "kegg.compound", "chebi", "metacyc.compound", "seed.compound",
        "hmdb", "lipidmaps", "reactome", "kegg.glycan", "sabiork.compound",
    )
    #: base BiGG id -> number of published BiGG models using it (policy 8).
    #: Built from the *raw* BiGG distribution: the per-compartment rows must
    #: be unioned, since collapsing them loses most of the model breadth.
    bigg_model_count: dict = field(default_factory=dict, repr=False)
    _evidence: dict = field(default_factory=dict, repr=False)
    _joined: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def joined(self) -> pl.DataFrame:
        if self._joined is None:
            self._joined = self._build_joined()
        return self._joined

    def _build_joined(self) -> pl.DataFrame:
        m = self.mets.with_columns(
            pl.col(self.id_col).str.extract(self.compartment_pattern, 1).alias("compartment"),
            pl.col(self.id_col).str.replace(self.compartment_pattern, "").alias("base_id"),
        )
        rename = {"id": "metanetx.chemical", "name": "ref_name",
                  "formula": "ref_formula", "charge": "ref_charge",
                  "inchikey": "ref_inchikey", "bigg.metabolite": "ref_bigg",
                  "smiles": "ref_smiles"}
        # Each reference joins only if present. A missing table leaves its
        # ref_* columns absent, and every check that consumes them tests for
        # the column first, so an offline review degrades instead of failing.
        ref = None
        if self.mnx is not None:
            keep = [c for c in rename if c in self.mnx.columns]
            ref = self.mnx.select(keep).rename(
                {k: v for k, v in rename.items() if k in keep})
            for c in self.xref_columns:
                if c in self.mnx.columns:
                    ref = ref.join(
                        self.mnx.select(pl.col("id").alias("metanetx.chemical"),
                                        pl.col(c).alias("ref_" + c)),
                        on="metanetx.chemical", how="left")
        brename = {"id": "base_id", "name": "bigg_name",
                   "inchikey": "bigg_inchikey", "old_bigg_ids": "bigg_old_ids"}
        bg = None
        if self.bigg is not None:
            bkeep = [c for c in brename if c in self.bigg.columns]
            bg = (self.bigg.select(bkeep)
                  .rename({k: v for k, v in brename.items() if k in bkeep})
                  .unique("base_id"))

        j = m
        if ref is not None:
            j = j.join(ref, on="metanetx.chemical", how="left")
        if bg is not None:
            j = j.join(bg, on="base_id", how="left")
        # Keep the joined schema constant whether or not a reference was
        # supplied: absent tables contribute all-null columns. Checks then
        # find nothing to compare against and stay silent, instead of every
        # consumer needing its own hasattr-style guard.
        placeholders = ["ref_name", "ref_formula", "ref_charge", "ref_inchikey",
                        "ref_bigg", "ref_smiles", "bigg_name", "bigg_inchikey",
                        "bigg_old_ids"]
        placeholders += ["ref_" + c for c in self.xref_columns]
        missing = [p for p in placeholders if p not in j.columns]
        if missing:
            j = j.with_columns([pl.lit(None, dtype=pl.String).alias(p)
                                for p in missing])
        return j.with_columns(
            pl.col("formula").map_elements(norm_formula, return_dtype=pl.String).alias("_nf"),
            pl.col("ref_formula").map_elements(norm_formula, return_dtype=pl.String).alias("_nf_ref"),
            pl.col("name").map_elements(norm_name, return_dtype=pl.String).alias("_nn"),
            pl.col("ref_name").map_elements(norm_name, return_dtype=pl.String).alias("_nn_ref"),
            pl.col("inchikey").map_elements(inchikey_skeleton, return_dtype=pl.String).alias("_skel"),
            pl.col("bigg_inchikey").map_elements(inchikey_skeleton, return_dtype=pl.String).alias("_skel_bigg"),
        )

    def api_evidence(self, row: dict) -> dict:
        """Independent formula/charge evidence for one metabolite row.

        Thin adapter over references.apis.annotation_evidence: pulls the
        identifiers out of the model row in the propagation order the curator
        specified (identifiers before names) and hands them to the APIs.
        Results are memoised per metabolite id — a review run touches the same
        species from several checks.
        """
        from apis import annotation_evidence          # noqa: PLC0415
        mid = row[self.id_col]
        if mid in self._evidence:
            return self._evidence[mid]
        chebi = next(iter(sorted(split_multi(row.get("chebi")))), None)
        ev = annotation_evidence(
            self.apis, inchikey=row.get("inchikey"), chebi_id=chebi,
            name=row.get("name"))
        self._evidence[mid] = ev
        return ev

    def prefetch_evidence(self, rows) -> None:
        """Warm the evidence cache for many metabolites at once.

        Each ``api_evidence`` call is several dependent HTTP round-trips, so
        doing them one metabolite at a time inside a check loop makes API
        latency the entire runtime of the review. ``prefetch_annotation_
        evidence`` issues them as concurrent waves on an event loop — see
        ``references/unichem.py`` for why async and not threads.
        """
        if self.apis is None:
            return
        todo = [r for r in rows if r[self.id_col] not in self._evidence]
        if not todo:
            return
        from apis import prefetch_annotation_evidence   # noqa: PLC0415
        try:
            self._evidence.update(prefetch_annotation_evidence(
                self.apis, todo, id_key=self.id_col))
        except Exception:
            pass    # a check that needs the evidence will fetch it itself

    def apis_is_fantasy(self, row: dict) -> bool:
        """True when no structural database knows this compound (policy 1/9)."""
        from apis import is_probable_fantasy          # noqa: PLC0415
        model_xrefs = {c: row.get(c) for c in self.independent_namespaces
                       if c not in self.pipeline_namespaces}
        model_xrefs.update({k: row.get(k) for k in ("inchikey", "smiles")})
        return is_probable_fantasy(self.api_evidence(row),
                                   model_xrefs=model_xrefs)

    @property
    def deprecated_map(self) -> dict[str, list[str]]:
        if self.depr is None:
            return {}
        g = self.depr.group_by("old").agg(pl.col("new").unique())
        return dict(zip(g["old"].to_list(), g["new"].to_list()))

    @property
    def valid_mnx(self) -> set[str]:
        if self.mnx is None:
            return set()
        return set(self.mnx["id"].to_list())

    @property
    def bigg_ids(self) -> set[str]:
        if self.bigg is None:
            return set()
        return set(self.bigg["id"].to_list())

    @property
    def bigg_old_to_new(self) -> dict[str, set[str]]:
        if self.bigg is None or "old_bigg_ids" not in self.bigg.columns:
            return {}
        e = (self.bigg.select("id", "old_bigg_ids").drop_nulls()
             .with_columns(pl.col("old_bigg_ids").str.split(";")).explode("old_bigg_ids")
             .with_columns(pl.col("old_bigg_ids").str.strip_chars()))
        out: dict[str, set[str]] = defaultdict(set)
        for old, new in zip(e["old_bigg_ids"].to_list(), e["id"].to_list()):
            if old:
                out[old].add(new)
        return dict(out)


def bigg_model_counts(path) -> dict[str, int]:
    """base BiGG id -> number of distinct models using it.

    Read this from the raw ``bigg_models_metabolites.txt``: the file carries
    one row per (metabolite, compartment) and each row lists only the models
    using *that* compartment form. Deduplicating on ``universal_bigg_id``
    before counting — the obvious thing to do — silently keeps one arbitrary
    compartment and understates breadth by an order of magnitude (stcoa
    reads as 1 model instead of 87), which inverts the policy-8 preference.
    """
    raw = pl.read_csv(path, separator="\t", infer_schema_length=0,
                      quote_char=None)
    agg = (raw.with_columns(pl.col("model_list").fill_null("").str.split("; "))
           .group_by("universal_bigg_id")
           .agg(pl.col("model_list").flatten().unique().alias("models"))
           .with_columns(pl.col("models").list.len().alias("n")))
    return dict(zip(agg["universal_bigg_id"].to_list(), agg["n"].to_list()))


# --------------------------------------------------------------------------
# finding accumulator
# --------------------------------------------------------------------------

#: Where a finding goes. The curator's worklist must contain only findings a
#: human can act on *now*; everything else is routed to its own file.
#:   review       -> curator worklist (round 1)
#:   annotation   -> annotation-update file, for human review before round 2
#:   balance      -> held-out balance file, actioned in round 2 only
#:   observation  -> model-property report; true statements, not errors
ROUTES = ("review", "annotation", "balance", "observation")

#: Tags that describe a *property* of the model rather than a mistake in it.
#: A dead-end metabolite is usually a real gap in the reconstruction, not a bad
#: metabolite; a placeholder formula on a macromolecule is a modelling
#: convention; a reference entry with no structure is a gap in the reference,
#: not in the model. All four are true and none is a curator action, and
#: together they were 808 of the 1104 worklist rows on the test model — a
#: worklist that is 73% non-actionable is one a curator stops reading.
#:
#: Checked against the scored round-1 cases before demoting: no confirmed
#: spiked error is reachable *only* through one of these tags.
OBSERVATION_TAGS = frozenset({
    "dead_end_metabolite", "placeholder_formula", "reference_lacks_structure",
    "orphan_metabolite",
})

_FINDING_SCHEMA = {
    "id": pl.String, "error_tag": pl.String, "confidence": pl.String,
    "check": pl.String, "rationale": pl.String, "proposed_fix": pl.String,
    "fix_verb": pl.String, "fix_arg": pl.String,
    "evidence": pl.String, "route": pl.String, "reaction_flag": pl.String,
    "reaction_note": pl.String,
}


class Findings:
    """Accumulates metabolite findings plus the reaction-reviewer handoff.

    ``reaction_flag`` / ``reaction_note`` implement policy 10: when a
    metabolite-level fix cannot resolve a finding, the finding still records
    *which reactions* the downstream reaction_reviewer agent must look at and
    *why*, rather than silently proposing a fix that would be wrong.
    """

    def __init__(self) -> None:
        self.rows: list[dict] = []

    def add(self, met_id, tag, confidence, rationale, proposed_fix="",
            evidence="", check="", route="review",
            reaction_flag="", reaction_note="",
            fix_verb="", fix_arg=""):
        """Record a finding.

        ``proposed_fix`` is prose for the curator. ``fix_verb`` / ``fix_arg``
        are the same fix in the decision vocabulary the updates logbook
        speaks (``merge_into`` + surviving id, ``set_charge`` + value, ...).
        Both are kept: the curator reads the prose and may overrule it, but
        when they simply write ``accept`` the logbook needs a fix it can apply
        without parsing English. A check that cannot express its fix
        structurally leaves these blank, and ``accept`` on that row is
        reported back rather than guessed at.
        """
        assert tag in ERROR_TAGS, f"unknown tag {tag}"
        assert confidence in CONFIDENCE, f"bad confidence {confidence}"
        assert route in ROUTES, f"bad route {route}"
        # Observation tags are routed by tag, not by the call site, so a check
        # cannot accidentally put a model property on the worklist.
        if tag in OBSERVATION_TAGS and route == "review":
            route = "observation"
        self.rows.append({
            "id": met_id, "error_tag": tag, "confidence": confidence,
            "check": check, "rationale": rationale,
            "proposed_fix": proposed_fix,
            "fix_verb": fix_verb, "fix_arg": fix_arg, "evidence": evidence,
            "route": route,
            "reaction_flag": reaction_flag, "reaction_note": reaction_note,
        })

    def frame(self) -> pl.DataFrame:
        if not self.rows:
            return pl.DataFrame(schema=_FINDING_SCHEMA)
        return pl.DataFrame(self.rows, schema_overrides=_FINDING_SCHEMA)

    def reaction_handoff(self, model) -> pl.DataFrame:
        """One row per (reaction, reason) for the reaction_reviewer agent."""
        # Guard the argument type explicitly. This used to accept anything and
        # swallow AttributeError per metabolite, so passing a ReviewContext
        # instead of a model produced an empty handoff file and no error at
        # all — 424 rows of downstream work silently became zero.
        if not hasattr(model, "metabolites"):
            raise TypeError(
                "reaction_handoff needs a cobra model, got "
                f"{type(model).__name__}; pass ctx.model, not ctx")
        by_met: dict[str, list[dict]] = defaultdict(list)
        for r in self.rows:
            if r["reaction_flag"]:
                by_met[r["id"]].append(r)
        rows = []
        for met_id, findings in by_met.items():
            try:
                met = model.metabolites.get_by_id(met_id)
            except KeyError:
                continue        # finding is not about a species (e.g. schema)
            for rxn in sorted(met.reactions, key=lambda r: r.id):
                for fnd in findings:
                    rows.append({
                        "reaction": rxn.id, "reaction_name": rxn.name,
                        "reaction_string": rxn.reaction,
                        "metabolite": met_id,
                        "flag": fnd["reaction_flag"],
                        "source_tag": fnd["error_tag"],
                        "note": fnd["reaction_note"] or fnd["rationale"],
                    })
        # Sort on every column, not just (reaction, metabolite): ties on those
        # two are common (one metabolite, several flags, same reaction) and a
        # partial sort leaves the tied block in set-iteration order, so two
        # runs of the same model emit byte-different files.
        if not rows:
            return pl.DataFrame()
        out = pl.DataFrame(rows).unique()
        return out.sort(*out.columns)


# --------------------------------------------------------------------------
# CHECK 1 — MetaNetX identifier hygiene
# --------------------------------------------------------------------------

def check_metanetx_ids(ctx: ReviewContext, f: Findings) -> None:
    """MetaNetX identifier hygiene.

    Policy 5: these are *not* curator errors — a missing or stale MNXM is
    annotation work. Every finding here is routed to the annotation-update
    file for human review rather than onto the worklist.
    """
    if ctx.mnx is None:
        return
    valid, depr = ctx.valid_mnx, ctx.deprecated_map
    for row in ctx.joined.iter_rows(named=True):
        mid, raw = row[ctx.id_col], row.get("metanetx.chemical")
        if raw is None:
            f.add(mid, "missing_metanetx_id", "high",
                  "no metanetx.chemical annotation; identifier cannot be validated",
                  "map from InChIKey, then name, then formula+charge",
                  check="metanetx_ids", route="annotation")
            continue
        ids = split_multi(raw)
        if len(ids) > 1:
            f.add(mid, "multiple_metanetx_id", "medium",
                  f"{len(ids)} MetaNetX ids assigned to one metabolite",
                  "keep the id whose formula and charge match the model",
                  evidence="|".join(sorted(ids)), check="metanetx_ids",
                  route="annotation")
        for i in sorted(ids):
            if i in valid:
                continue
            if i in depr:
                f.add(mid, "deprecated_metanetx_id", "high",
                      f"{i} is deprecated in the current MetaNetX release",
                      f"replace with {','.join(depr[i][:3])}"
                      + (" (pick by formula/charge match)" if len(depr[i]) > 1 else ""),
                      evidence=f"{i} -> {','.join(depr[i])}",
                      check="metanetx_ids", route="annotation")
            else:
                f.add(mid, "unknown_metanetx_id", "high",
                      f"{i} is neither current nor in the deprecation table",
                      "re-map from InChIKey or name", evidence=i,
                      check="metanetx_ids", route="annotation")


# --------------------------------------------------------------------------
# CHECK 2 — annotation vs reference (weak: block is usually MNXM-derived)
# --------------------------------------------------------------------------

def check_reference_agreement(ctx: ReviewContext, f: Findings) -> None:
    for row in ctx.joined.iter_rows(named=True):
        mid = row[ctx.id_col]
        if row["_nf_ref"] and row["_nf"] and row["_nf"] != row["_nf_ref"]:
            f.add(mid, "formula_mismatch_reference", "high",
                  f"formula {row['formula']} != MetaNetX {row['ref_formula']}",
                  f"set formula to {row['ref_formula']}",
                  fix_verb="set_formula", fix_arg=str(row["ref_formula"]),
                  evidence=row.get("metanetx.chemical", ""),
                  check="reference_agreement")
        rc = row.get("ref_charge")
        if (rc is not None and row["charge"] is not None
                and int(rc) != int(row["charge"]) and row["_nf"] != "H1"):
            f.add(mid, "charge_mismatch_reference", "high",
                  f"charge {row['charge']} != MetaNetX {int(rc)}",
                  f"set charge to {int(rc)}",
                  fix_verb="set_charge", fix_arg=str(int(rc)),
                  evidence=row.get("metanetx.chemical", ""),
                  check="reference_agreement")
        if (row.get("ref_inchikey") and row.get("inchikey")
                and row["inchikey"] != row["ref_inchikey"]):
            f.add(mid, "inchikey_mismatch_reference", "medium",
                  f"InChIKey {row['inchikey']} != MetaNetX {row['ref_inchikey']}",
                  "adopt the reference InChIKey once the MNXM is confirmed",
                  fix_verb="set_annotation",
                  fix_arg=f"inchikey={row['ref_inchikey']}",
                  check="reference_agreement")


def check_xref_agreement(ctx: ReviewContext, f: Findings) -> None:
    j = ctx.joined
    for col in ctx.xref_columns:
        rcol = "ref_" + col
        if col not in j.columns or rcol not in j.columns:
            continue
        for row in j.iter_rows(named=True):
            a, b = split_multi(row.get(col)), split_multi(row.get(rcol))
            if not a or not b:
                continue
            if not ({x.lower() for x in a} & {x.lower() for x in b}):
                f.add(row[ctx.id_col], "xref_disagreement", "medium",
                      f"{col} {sorted(a)} not among MetaNetX values {sorted(b)[:4]}",
                      f"reconcile {col} against the assigned MNXM",
                      evidence=f"{col}={sorted(a)}", check="xref_agreement")


# --------------------------------------------------------------------------
# CHECK 3 — BiGG namespace (independent of the MNXM assignment)
# --------------------------------------------------------------------------

def check_bigg_namespace(ctx: ReviewContext, f: Findings) -> None:
    # No BiGG table means "cannot judge", not "nothing is in BiGG". Without
    # this the empty reference makes every id in the model look absent from
    # the namespace and the check emits one finding per metabolite.
    if ctx.bigg is None:
        return
    known, old2new = ctx.bigg_ids, ctx.bigg_old_to_new
    for row in ctx.joined.iter_rows(named=True):
        mid, base = row[ctx.id_col], row["base_id"]
        if base in known:
            continue
        if base in old2new:
            f.add(mid, "bigg_id_deprecated", "high",
                  f"{base} is an old BiGG id",
                  f"rename to {','.join(sorted(old2new[base]))}",
                  check="bigg_namespace")
        else:
            f.add(mid, "bigg_id_not_in_namespace", "low",
                  f"{base} is not in the BiGG universal namespace "
                  "(legitimate for organism-specific species)",
                  "confirm the id is intentional and documented",
                  check="bigg_namespace")


def check_bigg_structure_conflict(ctx: ReviewContext, f: Findings) -> None:
    """BiGG's InChIKey is not derived from the model's MNXM, so a skeleton
    conflict is genuine evidence that the id and the structure disagree."""
    for row in ctx.joined.iter_rows(named=True):
        s, sb = row.get("_skel"), row.get("_skel_bigg")
        if not s or not sb or s == sb:
            continue
        same_name = bool(row.get("_nn")) and row["_nn"] == norm_name(row.get("bigg_name"))
        f.add(row[ctx.id_col], "inchikey_skeleton_conflict_bigg",
              "low" if same_name else "medium",
              f"model InChIKey skeleton {s} != BiGG {sb} for id {row['base_id']}"
              + (" (names agree — likely a tautomer/anomer choice)" if same_name else ""),
              "decide whether the BiGG id or the structure is authoritative",
              evidence=f"{row.get('inchikey')} vs {row.get('bigg_inchikey')}",
              check="bigg_structure_conflict")


def check_inchikey_validity(ctx: ReviewContext, f: Findings) -> None:
    """Validate InChIKeys against UniChem before believing any of them.

    Policy 9. A model can carry an InChIKey that no chemistry database has
    ever issued — a typo, a truncation, or a key computed from a wrong
    structure. Every downstream check that groups or compares on InChIKey
    inherits that error silently, so the key is validated first.

    When the model key and the reference (MetaNetX/BiGG) key differ and both
    validate, the reference key wins: it is the one a curated database
    actually issued for this identifier.
    """
    if ctx.apis is None:
        return
    from apis import unichem_validate_inchikeys       # noqa: PLC0415

    j = ctx.joined
    keys = set()
    for col in ("inchikey", "ref_inchikey", "bigg_inchikey"):
        if col in j.columns:
            keys |= {k for k in j[col].drop_nulls().to_list() if k}
    if not keys:
        return
    try:
        valid = unichem_validate_inchikeys(ctx.apis, sorted(keys))
    except Exception:
        return                                    # network failure: stay silent
    ctx._inchikey_valid = valid

    for row in j.iter_rows(named=True):
        own = row.get("inchikey")
        if not own:
            continue
        reference = row.get("ref_inchikey") or row.get("bigg_inchikey")
        if valid.get(own):
            if reference and reference != own and valid.get(reference):
                f.add(row[ctx.id_col], "inchikey_mismatch_reference", "medium",
                      f"model InChIKey {own} and reference key {reference} both "
                      "validate in UniChem but differ — they are different "
                      "structures, not a formatting difference",
                      f"prefer the reference key {reference}; it is the one a "
                      "curated database issued for this identifier",
                      evidence=f"model={own} reference={reference}",
                      check="inchikey_validity", route="annotation")
            continue
        replacement = (reference if reference and valid.get(reference) else None)
        f.add(row[ctx.id_col], "invalid_inchikey", "high",
              f"InChIKey {own} is not known to UniChem — no chemistry database "
              "has issued this key, so every comparison based on it is void",
              (f"replace with the validated reference key {replacement}"
               if replacement else
               "recompute the key from the structure, or remove it"),
              fix_verb="set_annotation" if replacement else "",
              fix_arg=f"inchikey={replacement}" if replacement else "",
              evidence=f"inchikey={own}", check="inchikey_validity",
              route="annotation")


# --------------------------------------------------------------------------
# CHECK 4 — duplicates within a compartment
# --------------------------------------------------------------------------

def distinct_namespace_ids(ctx: ReviewContext, ids: list[str]) -> list[str]:
    """Namespaces in which the candidate duplicates hold *different* ids.

    Policy 7 — the single highest-yield false-positive suppressor found in
    round 1. Species sharing a formula, a name and even a MetaNetX id are
    still distinct entities if an independent curated namespace (KEGG,
    ChEBI, SEED, MetaCyc, ...) assigns them different identifiers: those
    databases made a deliberate distinction that the model should keep.

    Scored on round 1: every one of the 27 duplicate false positives carries
    such a conflict; not one of the 16 true positives does.
    """
    rows = {r[ctx.id_col]: r
            for r in ctx.joined.filter(pl.col(ctx.id_col).is_in(ids)).iter_rows(named=True)}
    conflicts = []
    for ns in ctx.independent_namespaces:
        if ns not in ctx.joined.columns:
            continue
        sets = [split_multi(rows.get(i, {}).get(ns)) for i in ids]
        sets = [s for s in sets if s]
        if len(sets) < 2:
            continue
        if not set.intersection(*sets):
            conflicts.append(f"{ns}=" + "|".join("/".join(sorted(s)) for s in sets))
    return conflicts


def check_duplicates(ctx: ReviewContext, f: Findings) -> None:
    """Species that are the same compound in the same compartment.

    Three grouping keys of decreasing strength. A pair matching on more than
    one of them must be reported *once*, on the strongest — round 1 emitted
    each such pair three times and inflated the worklist threefold.
    """
    j = ctx.joined
    keys = [("metanetx.chemical", "high", "same MetaNetX id"),
            ("inchikey", "high", "identical InChIKey"),
            ("_skel", "low", "same InChIKey skeleton "
                             "(stereo/tautomer differences expected)")]
    seen: set[tuple] = set()
    for key, conf, label in keys:
        if key not in j.columns:
            continue
        groups = (j.filter(pl.col(key).is_not_null())
                  .group_by([key, "compartment"])
                  .agg(pl.col(ctx.id_col).alias("ids"), pl.len().alias("n"))
                  .filter(pl.col("n") > 1))
        for row in groups.iter_rows(named=True):
            ids = row["ids"]
            signature = tuple(sorted(ids))
            if signature in seen:
                continue
            seen.add(signature)
            conflicts = distinct_namespace_ids(ctx, ids)
            for mid in ids:
                others = [x for x in ids if x != mid]
                if conflicts:
                    # Policy 7: chemically identical, but curated databases
                    # treat them as separate entities. Not a duplicate.
                    f.add(mid, "chemically_identical_but_distinct_in_databases",
                          "medium",
                          f"shares {label.split(' (')[0]} with {others} in "
                          f"compartment {row['compartment']}, but independent "
                          f"databases assign distinct ids: {'; '.join(conflicts)}",
                          "keep both species; do not merge. Verify the model "
                          "distinguishes them for a documented reason",
                          evidence="; ".join(conflicts),
                          check=f"duplicates:{key}")
                    continue
                prose, verb, arg = duplicate_fix(ctx, mid, others)
                f.add(mid, "duplicate_metabolite", conf,
                      f"{label} as {others} in compartment {row['compartment']}",
                      prose, fix_verb=verb, fix_arg=arg,
                      evidence=f"{key}={row[key]}", check=f"duplicates:{key}",
                      reaction_flag="duplicate_metabolite_merge",
                      reaction_note=f"{mid} is a duplicate of {others}; after "
                                    "the merge these reactions may themselves "
                                    "become duplicates and need review")


def duplicate_fix(ctx: ReviewContext, mid: str,
                  others: list[str]) -> tuple[str, str, str]:
    """Which of the duplicate ids to keep (policy 8).

    Returns ``(prose, fix_verb, fix_arg)``. When no candidate is a BiGG
    identifier there is no defensible survivor, so the verb is empty and the
    curator must name one — the alternative is picking arbitrarily and
    recording it in the logbook as though it were reasoned.

    Preference order, in decreasing strength:
      1. an id that is a real BiGG identifier over one that is not;
      2. among real BiGG ids, the one used by more published models.
    Reproduces the curator's three round-1 calls exactly: hicit over hisocit
    (only hicit is BiGG), stcoa over ocdccoa (87 vs 2 models), ocdcea over
    ocdce9a (95 vs 1).
    """
    cands = [mid] + list(others)
    scored = []
    for c in cands:
        base = re.sub(ctx.compartment_pattern, "", c)
        counts = ctx.bigg_model_count or {}
        scored.append((base in ctx.bigg_ids, counts.get(base, 0), c, base))
    scored.sort(reverse=True)
    keep = scored[0]
    drop = [s[2] for s in scored[1:]]
    if not keep[0]:
        return ("merge the species and redirect the reactions; no candidate is "
                "a BiGG identifier, so the surviving id needs a curator "
                "decision", "", "")
    why = (f"{keep[3]} is a BiGG identifier"
           + (f" used by {keep[1]} models" if keep[1] else ""))
    rival = [s for s in scored[1:] if s[0]]
    if rival:
        why += f" vs {rival[0][3]} in {rival[0][1]}"
    prose = f"keep {keep[2]}, redirect the reactions of {drop} to it ({why})"
    # The finding is filed against `mid`; if `mid` is the survivor there is
    # nothing to do to it, and the merge is recorded on the row for the id
    # that disappears.
    verb, arg = (("merge_into", keep[2]) if keep[2] != mid else ("", ""))
    return prose, verb, arg


# --------------------------------------------------------------------------
# CHECK 4b — duplicates that carry NO shared identifier
#
# The MNXM/InChIKey duplicate check above is blind to the most common way a
# duplicate enters a model: a new species is added with an ad-hoc id and no
# annotation at all, so it shares no key with the species it duplicates. Two
# signals find these without any shared identifier.
# --------------------------------------------------------------------------

def _degree(ctx: ReviewContext) -> dict[str, int]:
    return {mt.id: len([r for r in mt.reactions if not r.boundary])
            for mt in ctx.model.metabolites}


def check_name_collisions(ctx: ReviewContext, f: Findings) -> None:
    """Two species in one compartment whose names normalise identically.

    Independent of every identifier, so it catches unannotated duplicates.
    Asymmetric reaction degree (one species wired into the network, the other
    into one or two reactions) is strong evidence of which is the intruder.
    """
    deg = _degree(ctx)
    j = ctx.joined.with_columns(
        pl.col("name").map_elements(norm_name, return_dtype=pl.String).alias("_nk"))
    groups = (j.filter(pl.col("_nk").is_not_null())
              .group_by(["_nk", "compartment"])
              .agg(pl.col(ctx.id_col).alias("ids"),
                   pl.col("formula").unique().alias("formulas"),
                   pl.len().alias("n"))
              .filter(pl.col("n") > 1))
    for row in groups.iter_rows(named=True):
        ids = row["ids"]
        degrees = {i: deg.get(i, 0) for i in ids}
        hub = max(degrees, key=degrees.get)
        same_formula = len(row["formulas"]) == 1
        for mid in ids:
            others = [x for x in ids if x != mid]
            asymmetric = degrees[mid] <= 2 and degrees[hub] >= 5
            conf = "high" if (same_formula and asymmetric) else (
                "medium" if same_formula else "low")
            note = ""
            if mid != hub and asymmetric:
                note = (f"; it carries {degrees[mid]} reaction(s) against "
                        f"{degrees[hub]} for {hub}, so {mid} is the likely intruder")
            f.add(mid, "name_collision", conf,
                  f"name normalises identically to {others} in compartment "
                  f"{row['compartment']}"
                  + ("; same formula" if same_formula else "; formulas differ")
                  + note,
                  (f"merge into {hub} and redirect its reactions"
                   if mid != hub and asymmetric
                   else "disambiguate the names or merge the species"),
                  evidence=f"name_key={row['_nk']}", check="name_collision")


def check_shadow_duplicates(ctx: ReviewContext, f: Findings,
                            max_degree: int = 2, hub_degree: int = 5) -> None:
    """Unannotated, barely-connected species that duplicate an annotated one.

    A species with no MetaNetX id, at most `max_degree` reactions, and the same
    molecular formula as a well-connected annotated species in the same
    compartment is almost always an accidental second copy.
    """
    deg = _degree(ctx)
    j = ctx.joined
    by_key: dict[tuple, list[dict]] = defaultdict(list)
    for row in j.iter_rows(named=True):
        by_key[(norm_formula(row["formula"]), row["compartment"])].append(row)
    for row in j.filter(pl.col("metanetx.chemical").is_null()).iter_rows(named=True):
        mid = row[ctx.id_col]
        if deg.get(mid, 0) > max_degree:
            continue
        peers = [p for p in by_key[(norm_formula(row["formula"]), row["compartment"])]
                 if p[ctx.id_col] != mid and p.get("metanetx.chemical")
                 and deg.get(p[ctx.id_col], 0) >= hub_degree]
        if not peers:
            continue
        hub = max(peers, key=lambda p: deg.get(p[ctx.id_col], 0))
        f.add(mid, "shadow_duplicate", "high",
              f"unannotated species with {deg.get(mid, 0)} reaction(s) has the "
              f"same formula ({row['formula']}) and compartment as "
              f"{hub[ctx.id_col]} ('{hub['name']}', "
              f"{deg.get(hub[ctx.id_col], 0)} reactions), which is annotated",
              f"merge into {hub[ctx.id_col]} and redirect its reactions",
              fix_verb="merge_into", fix_arg=hub[ctx.id_col],
              evidence=f"hub={hub[ctx.id_col]} mnxm={hub['metanetx.chemical']}",
              check="shadow_duplicate")


def parse_bigg_reaction(reaction_string: str) -> dict[str, float] | None:
    """BiGG reaction string -> {metabolite: signed coefficient}."""
    s = str(reaction_string)
    for arrow in ("<->", "-->", "<--", "->", "<-", "="):
        if arrow in s:
            left, right = s.split(arrow, 1)
            break
    else:
        return None
    out: dict[str, float] = defaultdict(float)
    for side, sign in ((left, -1.0), (right, 1.0)):
        for term in side.split(" + "):
            parts = term.split()
            if not parts:
                continue
            try:
                coefficient, met = float(parts[0]), parts[1]
            except (ValueError, IndexError):
                coefficient, met = 1.0, parts[0]
            out[met] += sign * coefficient
    return dict(out)


def align_stoichiometry(model_reaction, bigg_string):
    """Direction-agnostic alignment -> (shared, model_only, bigg_only)."""
    a = {m.id: c for m, c in model_reaction.metabolites.items()}
    b = parse_bigg_reaction(bigg_string)
    if not b:
        return None
    best = None
    for flip in (1.0, -1.0):
        bb = {k: flip * v for k, v in b.items()}
        shared = {k for k in a if k in bb and abs(a[k] - bb[k]) < 1e-6}
        cand = (len({k for k in a if k not in shared})
                + len({k for k in bb if k not in shared}),
                shared,
                {k: v for k, v in a.items() if k not in shared},
                {k: v for k, v in bb.items() if k not in shared})
        if best is None or cand[0] < best[0]:
            best = cand
    return best[1:]


def check_reaction_anchored_duplicates(ctx: ReviewContext, f: Findings,
                                       fuzzy_cutoff: float = 0.92) -> None:
    """Shadow duplicates found through the *reaction*, not the metabolite.

    The hardest duplicates carry a plausible chemical name, a plausible
    formula and a normal reaction degree, so every metabolite-level signal
    (name collision, formula grouping, degree asymmetry) misses them. What
    gives them away is the reaction: a reaction whose id is not a BiGG id but
    whose *name* matches a real BiGG reaction, whose stoichiometry then
    aligns one-for-one except for a single species — a non-BiGG metabolite
    standing where a real BiGG metabolite belongs.

    The discriminator matters as much as the match: the model-side species
    must be outside the BiGG namespace while its counterpart is inside it.
    Without it the check floods with cofactor variants (NAD/NADP, ATP/GTP)
    where both sides are legitimate BiGG ids.

    Recovers all three cases the curator flagged as undetectable by name or
    reaction context in round 1: hmcit_L_m -> hcit_m, gal1p_L_c -> gal1p_c,
    2oglutm_h -> 2ogm_h.
    """
    if ctx.bigg_rxns is None:
        return
    import difflib

    by_name: dict[str, list[str]] = defaultdict(list)
    strings: dict[str, str] = {}
    id_col = "id" if "id" in ctx.bigg_rxns.columns else "bigg_id"
    str_col = ("reaction" if "reaction" in ctx.bigg_rxns.columns
               else "reaction_string")
    for row in ctx.bigg_rxns.iter_rows(named=True):
        by_name[norm_name(row["name"]) or ""].append(row[id_col])
        strings[row[id_col]] = row[str_col]
    known_rxn = set(strings)
    name_keys = list(by_name)
    model_mets = {m.id for m in ctx.model.metabolites}

    for rxn in ctx.model.reactions:
        if rxn.boundary or rxn.id in known_rxn:
            continue
        key = norm_name(rxn.name) or ""
        candidates, fuzzy = by_name.get(key, []), False
        if not candidates:
            near = difflib.get_close_matches(key, name_keys, n=1,
                                             cutoff=fuzzy_cutoff)
            candidates = by_name.get(near[0], [])[:1] if near else []
            fuzzy = bool(candidates)
        for bigg_id in candidates[:2]:
            aligned = align_stoichiometry(rxn, strings[bigg_id])
            if not aligned:
                continue
            shared, model_only, bigg_only = aligned
            if not shared or not (0 < len(model_only) == len(bigg_only) <= 2):
                continue
            for mid, coefficient in model_only.items():
                partners = [k for k, v in bigg_only.items()
                            if abs(coefficient - v) < 1e-6]
                if len(partners) != 1:
                    continue
                target = partners[0]
                mid_base = re.sub(ctx.compartment_pattern, "", mid)
                tgt_base = re.sub(ctx.compartment_pattern, "", target)
                if mid_base in ctx.bigg_ids or tgt_base not in ctx.bigg_ids:
                    continue
                in_model = target in model_mets
                f.add(mid, "shadow_duplicate", "medium" if fuzzy else "high",
                      f"reaction {rxn.id} is not a BiGG reaction but its name "
                      f"matches BiGG {bigg_id}"
                      + (" (fuzzy name match)" if fuzzy else "")
                      + f"; the two stoichiometries agree on {len(shared)} "
                      f"species and differ only in {mid} standing where BiGG "
                      f"has {target}. {mid_base} is not a BiGG identifier, "
                      f"{tgt_base} is"
                      + ("; the target already exists in the model"
                         if in_model else ""),
                      f"substitute {target} for {mid}"
                      + (" and merge" if in_model else " (add the species)")
                      + f"; then reconcile {rxn.id} with {bigg_id}",
                      fix_verb="merge_into" if in_model else "replace_with",
                      fix_arg=target,
                      evidence=f"rxn={rxn.id} bigg_rxn={bigg_id} target={target}",
                      check="reaction_anchored_duplicate",
                      reaction_flag="non_bigg_reaction_duplicates_bigg",
                      reaction_note=f"{rxn.id} has no BiGG id and no xrefs but "
                                    f"its name matches BiGG {bigg_id} "
                                    f"({strings[bigg_id]}); confirm whether it "
                                    "should be replaced by the BiGG reaction")


def check_uncurated_name_stubs(ctx: ReviewContext, f: Findings) -> None:
    """A name that restates the id *with the compartment appended*
    (`hmcit-L[m]` for `hmcit_L_m`) means the species was never curated — no
    human ever gave it a chemical name.

    Note the deliberate narrowness: a name equal to the bare base id is NOT
    flagged, because legitimate abbreviations do that constantly (`NADPH` for
    `nadph_c`, `ATP` for `atp_h`). Only the id-plus-compartment form is
    diagnostic — no chemical name contains its own compartment tag.
    """
    for row in ctx.joined.iter_rows(named=True):
        base = re.sub(r"[^a-z0-9]", "", str(row["base_id"]).lower())
        name = re.sub(r"[^a-z0-9]", "", str(row["name"] or "").lower())
        comp = row["compartment"] or ""
        if name and comp and name == base + comp:
            # Scored round 1: all three stubs turned out to be shadow
            # duplicates whose partner was only findable through the
            # reaction. The stub itself is a weak signal; its value is that
            # it marks a species whose *reactions* need review.
            f.add(row[ctx.id_col], "uncurated_name_stub", "medium",
                  f"name '{row['name']}' only restates the identifier, so this "
                  "species was added without curation; in round 1 every such "
                  "stub proved to be a duplicate of a curated species reached "
                  "through the reaction rather than the metabolite",
                  "give it a chemical name and a MetaNetX id, or merge it into "
                  "the species it duplicates",
                  check="uncurated_name_stub",
                  reaction_flag="uncurated_species_in_reaction",
                  reaction_note=f"{row[ctx.id_col]} was never curated; check "
                                "whether these reactions duplicate a BiGG "
                                "reaction using the curated species instead")


def check_fantasy_metabolites(ctx: ReviewContext, f: Findings) -> None:
    """Species no structural database has ever heard of.

    A reconstruction pipeline can mint a BiGG id, a MetaNetX id and a SEED id
    for anything it likes — those three namespaces are *derived*, so their
    presence is not evidence that the compound exists. Evidence of existence
    means a curated structural record: ChEBI, KEGG, MetaCyc, HMDB, LipidMaps,
    an InChIKey or a SMILES.

    Two tiers, following the curator's two round-1 cases:

    probable_fantasy_metabolite
        annotated *only* in the pipeline namespaces, and the live APIs
        (UniChem / ChEBI / PubChem) return nothing for its name either.
        Round-1 case: cacoa_c, present in one model, no structural xref, no
        API hit for any variant of "coniferonyl-CoA".

    possible_fantasy_metabolite
        no structural xref and no API confirmation, *and* none of the
        reactions it participates in is a real BiGG reaction with a
        cross-reference of its own — so nothing anchors it to a real network.
        Round-1 cases: pchol1801835z9z12z_c, pchol1801845z9z12z15z_c.

    Both tiers also carry ``uncurated_metabolite`` and a reaction-reviewer
    flag: if the species is not real, neither are its reactions.
    """
    if ctx.apis is None or ctx.bigg is None:
        return
    structural = [c for c in ctx.independent_namespaces
                  if c not in ctx.pipeline_namespaces] + ["inchikey", "smiles",
                                                          "inchi"]
    present = [c for c in structural if c in ctx.joined.columns]
    bigg_rxn_ids = set()
    if ctx.bigg_rxns is not None:
        col = "id" if "id" in ctx.bigg_rxns.columns else "bigg_id"
        bigg_rxn_ids = set(ctx.bigg_rxns[col].to_list())

    # Warm the API cache in parallel before the serial loop: only the rows
    # that will actually be queried (no structural xref, non-systematic name).
    candidates = [r for r in ctx.joined.iter_rows(named=True)
                  if not any(r.get(c) for c in present)
                  and not SYSTEMATIC_NAME.search(r.get("name") or "")]
    ctx.prefetch_evidence(candidates)

    for row in ctx.joined.iter_rows(named=True):
        mid, name = row[ctx.id_col], row.get("name") or ""
        if any(row.get(c) for c in present):
            continue
        if set(parse_formula(row.get("formula"))) & PLACEHOLDER_ELEMENTS:
            continue  # generic/polymeric species are legitimately structureless

        base = re.sub(ctx.compartment_pattern, "", mid)
        in_namespace = base in ctx.bigg_ids
        n_models = ctx.bigg_model_count.get(base, 0)
        met = ctx.model.metabolites.get_by_id(mid)
        anchored = [r.id for r in met.reactions if r.id in bigg_rxn_ids]

        f.add(mid, "uncurated_metabolite", "high" if not in_namespace else "medium",
              "no cross-reference to any structural database (ChEBI, KEGG, "
              "MetaCyc, HMDB, LipidMaps), no InChIKey and no SMILES — it is "
              "annotated only in namespaces a reconstruction pipeline mints "
              "for itself",
              "assign a structural identifier, or delete the species if it "
              "cannot be substantiated",
              evidence=f"in_bigg={in_namespace} n_models={n_models}",
              check="fantasy", route="annotation" if in_namespace else "review")

        # Systematic / templated names are unsearchable by construction: no
        # database indexes every acyl permutation, every photon wavelength bin
        # or every ACP conjugate. Absence of a name hit for these is evidence
        # about the naming convention, not about the compound.
        if SYSTEMATIC_NAME.search(name):
            api_hit, api_note = None, "name is systematic; not name-searchable"
        elif ctx.apis is None:
            api_hit, api_note = None, "no API client configured"
        else:
            try:
                api_hit = not ctx.apis_is_fantasy(row)
                api_note = ("a structural database knows this compound by name"
                            if api_hit else
                            "no exact name match in ChEBI, PubChem or UniChem")
            except Exception as exc:                       # network/API failure
                api_hit, api_note = None, f"API lookup failed ({type(exc).__name__})"
        if api_hit:
            continue

        # Tier 1 — a trivial name that no database has ever issued, on a
        # species no more than one published model has ever used. Round-1
        # case: cacoa_c ("Coniferonyl-coa"), in BiGG but used by 1 model and
        # unknown to every structural database.
        if api_hit is False and n_models <= 1:
            f.add(mid, "probable_fantasy_metabolite", "high",
                  f"{api_note}, and only {n_models} published BiGG model uses "
                  "this identifier — the name is an ordinary chemical name, so "
                  "a real compound should have been found. It may not exist",
                  "confirm against the primary literature, or delete the "
                  "species and the reactions unique to it",
                  evidence=f"name={name!r} n_models={n_models} "
                           f"reactions={[r.id for r in met.reactions][:5]}",
                  check="fantasy",
                  reaction_flag="metabolite_may_not_exist",
                  reaction_note=f"{mid} has no structural database record; "
                                "these reactions are unverifiable until it is "
                                "confirmed or removed")
            continue

        # Tier 2 — outside the BiGG namespace with no BiGG reaction anchoring
        # it. Nothing ties the species to a curated network, but its name is
        # templated so absence of a database record proves less. Round-1
        # cases: pchol1801835z9z12z_c, pchol1801845z9z12z15z_c.
        if not in_namespace and not anchored:
            f.add(mid, "possible_fantasy_metabolite", "medium",
                  f"{base} is not a BiGG identifier, {api_note}, and none of "
                  f"its {len(met.reactions)} reaction(s) is a BiGG reaction — "
                  "nothing anchors this species to a curated network",
                  "confirm the species exists and is produced by this "
                  "organism, or delete it with its unique reactions",
                  evidence=f"reactions={[r.id for r in met.reactions][:5]}",
                  check="fantasy",
                  reaction_flag="metabolite_may_not_exist",
                  reaction_note=f"{mid} is outside BiGG and has no curated "
                                "reaction; confirm these reactions are real")


# --------------------------------------------------------------------------
# CHECK 4c — formula/charge plausibility, independent of the assignment
#
# Covers the blind spot where the assigned MNXM exists but has no formula or
# charge in chem_prop: the reference-agreement check silently skips those rows,
# and the modeller had to hand-fill the values (charge 0 is the usual default).
# --------------------------------------------------------------------------

def check_charge_plausibility(ctx: ReviewContext, f: Findings,
                              min_entries: int = 3) -> None:
    if ctx.mnx is None:
        return
    ref = (ctx.mnx.filter(pl.col("formula").is_not_null(),
                          pl.col("charge").is_not_null())
           .with_columns(pl.col("formula")
                         .map_elements(norm_formula, return_dtype=pl.String)
                         .alias("_k"))
           .group_by("_k").agg(pl.col("charge").unique().alias("charges"),
                               pl.len().alias("n")))
    charges = dict(zip(ref["_k"].to_list(), ref["charges"].to_list()))
    counts = dict(zip(ref["_k"].to_list(), ref["n"].to_list()))
    for row in ctx.joined.iter_rows(named=True):
        if set(parse_formula(row["formula"])) & PLACEHOLDER_ELEMENTS:
            continue
        key = norm_formula(row["formula"])
        allowed = charges.get(key)
        if allowed is None or row["charge"] is None:
            continue
        n = counts[key]
        if n < min_entries or float(row["charge"]) in {float(x) for x in allowed}:
            continue
        no_ref = row.get("ref_charge") is None
        allowed_sorted = sorted(float(x) for x in allowed)
        # Policy 1: a charge inferred from other entries sharing a formula is
        # still a calculated value. It may only be *applied* once the species
        # itself carries a database annotation; otherwise the fix is deferred
        # to the round that follows identifier assignment.
        has_identity = bool(row.get("metanetx.chemical")) and not no_ref
        api_charge, api_src = None, None
        if ctx.apis is not None and not has_identity:
            try:
                from apis import consensus_charge      # noqa: PLC0415
                api_charge, api_src = consensus_charge(ctx.api_evidence(row))
            except Exception:
                api_charge = None
        if api_charge is not None:
            fix = (f"set charge to {api_charge:g}, from {api_src}"
                   if api_charge != row["charge"] else "")
            if not fix:
                continue                      # the API vindicates the model
            deferral = ""
        elif has_identity:
            fix = (f"set charge to {allowed_sorted[0]:g} after confirming the "
                   "protonation state")
            deferral = ""
        else:
            fix = ""
            deferral = (" No database carries a charge for this species and it "
                        "has no usable identifier, so the correction is "
                        "deferred: a human must assign the identifier first, "
                        "then formula and charge are taken from it "
                        f"(arithmetic would suggest {allowed_sorted[0]:g})")
        f.add(row[ctx.id_col], "charge_implausible_for_formula",
              "high" if (no_ref and n >= 10) else "medium",
              f"charge {row['charge']} appears on no MetaNetX entry with formula "
              f"{row['formula']}; all {n} entries carry {allowed_sorted}"
              + (" — and the assigned MetaNetX id has no charge of its own, so "
                 "this value was hand-filled" if no_ref else "")
              + deferral,
              fix,
              evidence=f"formula={row['formula']} n_ref_entries={n}",
              check="charge_plausibility",
              route="review" if fix else "annotation",
              reaction_flag="" if fix else "charge_unresolved_no_reference",
              reaction_note="" if fix else
                            (f"{row[ctx.id_col]} carries an implausible charge "
                             "that no database can arbitrate; reactions using "
                             "it may appear charge-imbalanced for that reason"))


def check_reference_without_structure(ctx: ReviewContext, f: Findings) -> None:
    """MNXM entries lacking formula/charge make every reference comparison
    vacuous — flag them so the curator knows the row is unverified, not clean."""
    for row in ctx.joined.iter_rows(named=True):
        if not row.get("metanetx.chemical"):
            continue
        if row.get("ref_formula") is None and row.get("ref_charge") is None:
            f.add(row[ctx.id_col], "reference_lacks_structure", "low",
                  f"assigned {row['metanetx.chemical']} carries no formula or "
                  "charge in MetaNetX, so formula/charge here are unverified "
                  "rather than confirmed",
                  "verify formula and charge against MetaCyc/KEGG or a "
                  "structure source",
                  check="reference_without_structure")


# --------------------------------------------------------------------------
# CHECK 5 — cross-compartment consistency
# --------------------------------------------------------------------------

def check_cross_compartment(ctx: ReviewContext, f: Findings) -> None:
    grouped = (ctx.joined.group_by("base_id").agg(
        pl.col(ctx.id_col).alias("ids"),
        pl.col("formula").unique().alias("formulas"),
        pl.col("charge").unique().alias("charges"),
        pl.col("metanetx.chemical").unique().alias("mnxms"),
        pl.col("name").unique().alias("names"),
        pl.len().alias("n")).filter(pl.col("n") > 1))
    for row in grouped.iter_rows(named=True):
        problems = []
        if len(row["formulas"]) > 1:
            problems.append(f"formula {row['formulas']}")
        if len(row["charges"]) > 1:
            problems.append(f"charge {row['charges']}")
        if len(row["mnxms"]) > 1:
            problems.append(f"metanetx {row['mnxms']}")
        if len({norm_name(x) for x in row["names"]}) > 1:
            problems.append(f"name {row['names']}")
        if not problems:
            continue
        conf = "high" if any(p.startswith(("formula", "charge")) for p in problems) else "medium"
        for mid in row["ids"]:
            f.add(mid, "cross_compartment_inconsistency", conf,
                  f"the same base id {row['base_id']} differs across "
                  "compartments: " + "; ".join(problems),
                  "make all compartment instances of one species identical",
                  evidence=str(row["ids"]), check="cross_compartment")


# --------------------------------------------------------------------------
# CHECK 6 — balance-based localisation (the amplifier)
# --------------------------------------------------------------------------

def _imbalance(reaction) -> dict[str, float]:
    v: dict[str, float] = defaultdict(float)
    for met, coef in reaction.metabolites.items():
        for element, n in parse_formula(met.formula).items():
            v[element] += coef * n
        v["charge"] += coef * (met.charge or 0)
    return {k: round(x, 6) for k, x in v.items() if abs(x) > 1e-6}


def normalise_protons(model):
    """Return (model_copy, fixed_ids) with every proton charge set to +1."""
    fixed = []
    m = model.copy()
    for met in m.metabolites:
        if is_proton(met) and (met.charge or 0) != 1:
            met.charge = 1
            fixed.append(met.id)
    return m, fixed


def check_proton_convention(ctx: ReviewContext, f: Findings) -> list[str]:
    bad = [mt.id for mt in ctx.model.metabolites
           if is_proton(mt) and (mt.charge or 0) != 1]
    for mid in bad:
        f.add(mid, "proton_charge_convention", "high",
              "proton species carries charge != +1, which makes every reaction "
              "it appears in look charge-unbalanced",
              "set charge to +1 and formula to H",
              fix_verb="set_charge", fix_arg="1", check="proton_convention")
    return bad


def check_placeholder_formulas(ctx: ReviewContext, f: Findings) -> None:
    for mt in ctx.model.metabolites:
        used = set(parse_formula(mt.formula)) & PLACEHOLDER_ELEMENTS
        if used:
            f.add(mt.id, "placeholder_formula", "low",
                  f"formula {mt.formula} contains placeholder element(s) "
                  f"{sorted(used)}; reactions using it cannot be mass-balanced",
                  "acceptable for macromolecules/pseudo-metabolites — document it",
                  check="placeholder_formula")


def localise_imbalance(model, min_coefficient: float = 0.05,
                       min_reactions: int = 2) -> list[dict]:
    """Metabolites whose own formula/charge error explains the imbalance exactly.

    If a single metabolite carries error delta D, every reaction it appears in
    with coefficient c is imbalanced by exactly c*D. Requiring one consistent
    integer D = v_r / c_r across all of its reactions localises the error to
    that metabolite.

    Reactions containing placeholder-formula metabolites are excluded — their
    imbalance is a modelling convention, not an error. Boundary reactions are
    excluded for the same reason.
    """
    unbalanced = {}
    for r in model.reactions:
        if r.boundary:
            continue
        v = _imbalance(r)
        if v:
            unbalanced[r.id] = v

    out = []
    for met in model.metabolites:
        reactions = [r for r in met.reactions
                     if not r.boundary
                     and not any(has_placeholder(x) for x in r.metabolites)
                     and abs(r.metabolites[met]) >= min_coefficient]
        if len(reactions) < min_reactions:
            continue
        deltas = []
        for r in reactions:
            c = Fraction(r.metabolites[met]).limit_denominator(1000)
            if c == 0:
                deltas = None
                break
            v = unbalanced.get(r.id, {})
            deltas.append({k: Fraction(x).limit_denominator(1000) / c
                           for k, x in v.items()})
        if not deltas:
            continue
        base = deltas[0]
        if not base or any(d != base for d in deltas[1:]):
            continue
        if any(x.denominator != 1 for x in base.values()):
            continue
        # Sort the delta keys: this dict is rendered into human-readable
        # rationale strings, and unordered keys make two runs of the same
        # model produce byte-different files that the curator has to diff.
        out.append({"id": met.id, "name": met.name, "formula": met.formula,
                    "charge": met.charge, "n_reactions": len(reactions),
                    "delta": {k: int(base[k]) for k in sorted(base)}})
    return sorted(out, key=lambda h: h["id"])


def apply_delta(formula: str | None, charge, delta: dict[str, int]):
    """Formula/charge implied by subtracting the localised delta."""
    counts = parse_formula(formula)
    for element, d in delta.items():
        if element == "charge":
            continue
        counts[element] = counts.get(element, 0) - d
    if any(v < 0 for v in counts.values()):
        return None, None
    counts = {k: v for k, v in counts.items() if v}
    new_formula = "".join(
        f"{e}{counts[e] if counts[e] != 1 else ''}"
        for e in sorted(counts, key=lambda x: (x != "C", x != "H", x)))
    return new_formula, (charge or 0) - delta.get("charge", 0)


def check_balance_localisation(ctx: ReviewContext, f: Findings,
                               proton_fixed: list[str]) -> list[dict]:
    """Localise imbalance to single metabolites, then let the databases judge.

    Policy 1-4. The arithmetic tells you *which* metabolite is implicated and
    *what* formula/charge would close the balance. It does not tell you that
    the metabolite is wrong. Three things can be true:

      * the database confirms the implied value -> a real, applicable fix;
      * the database gives a different value, or the implied charge is not a
        state the compound can hold -> the metabolite is right and the error
        is in the *reaction*: flag it for the reaction_reviewer, propose
        nothing (the malonyl-CoA case: balance implies -4, ChEBI says the
        species is only ever 0 or -5);
      * no database value is available -> the fix is deferred until a human
        assigns a correct identifier (the gal1p_L_c case).

    Two further guards learned from round 1:

      * paired metabolites carrying equal and opposite imbalance (a redox
        couple such as ficytc_m / focytc_m) cancel: applying both fixes
        reproduces the original imbalance. These are reaction-level errors
        by construction and are never proposed as metabolite fixes.
      * every finding here is routed to ``balance``, held out of the round-1
        worklist entirely.
    """
    model = ctx.model
    if proton_fixed:
        model, _ = normalise_protons(model)
    hits = localise_imbalance(model)
    ref = {r[ctx.id_col]: r for r in ctx.joined.iter_rows(named=True)}
    cancelling = _cancelling_pairs(hits)
    ctx.prefetch_evidence([ref[h["id"]] for h in hits
                           if h["id"] in ref and h["id"] not in cancelling])

    for h in hits:
        mid = h["id"]
        delta = h["delta"]
        mass_part = {k: v for k, v in delta.items() if k != "charge"}
        new_formula, new_charge = apply_delta(h["formula"], h["charge"], delta)
        row = ref.get(mid, {})
        tag = "mass_imbalance_localised" if mass_part else "charge_imbalance_localised"
        base_rationale = (
            f"all {h['n_reactions']} balanceable reactions of this metabolite "
            f"are imbalanced by exactly its stoichiometric coefficient times "
            f"{delta}, which localises the arithmetic to this metabolite")
        evidence = f"delta={delta}; current={h['formula']}/{h['charge']}"

        if mid in cancelling:
            partner = cancelling[mid]
            f.add(mid, tag, "high",
                  base_rationale + f", but {partner} carries the equal and "
                  "opposite imbalance: correcting both would reproduce the "
                  "same imbalance, so the error is in the reaction, not in "
                  "either metabolite",
                  "", evidence=evidence + f"; cancels_with={partner}",
                  check="balance_localisation", route="balance",
                  reaction_flag="imbalance_cancels_between_metabolites",
                  reaction_note=f"{mid} and {partner} carry equal and opposite "
                                f"imbalance {delta}; the stoichiometry of these "
                                "reactions is the likely error")
            continue

        # --- what do the databases say? -----------------------------------
        verdict, source, detail = _database_verdict(
            ctx, row, new_formula, new_charge, bool(mass_part))

        if verdict == "confirms":
            f.add(mid, tag, "high",
                  base_rationale + f", and {source} independently carries the "
                  "same value",
                  f"formula -> {new_formula}, charge -> {new_charge} "
                  f"(confirmed by {source})",
                  evidence=evidence + f"; {detail}",
                  check="balance_localisation", route="balance")
        elif verdict == "contradicts":
            f.add(mid, tag, "high",
                  base_rationale + f", but {detail}. The database annotation "
                  "outranks the arithmetic, so this metabolite is not the "
                  "error — the reaction is",
                  "", evidence=evidence + f"; {detail}",
                  check="balance_localisation", route="balance",
                  reaction_flag="imbalance_not_explained_by_metabolite",
                  reaction_note=f"balance implies {new_formula}/{new_charge} for "
                                f"{mid}, which {detail}; check the "
                                "stoichiometry and the other participants")
        else:  # no database value to arbitrate with
            f.add(mid, tag, "medium",
                  base_rationale + ", but no structural database carries a "
                  "formula or charge for it, so the implied value cannot be "
                  "confirmed",
                  f"defer: assign a structural identifier first, then take "
                  f"formula/charge from it (arithmetic would give "
                  f"{new_formula}/{new_charge})",
                  evidence=evidence + "; no reference value",
                  check="balance_localisation", route="balance",
                  reaction_flag="imbalance_unresolved_no_reference",
                  reaction_note=f"{mid} has no reference formula/charge; "
                                "the imbalance in these reactions cannot be "
                                "attributed until it is annotated")
    return hits


def _cancelling_pairs(hits: list[dict]) -> dict[str, str]:
    """Metabolites whose localised imbalance is equal and opposite.

    Correcting both members reproduces the original imbalance, so neither is
    the error. Round-1 case: ficytc_m / focytc_m.
    """
    by_delta: dict[tuple, list[str]] = defaultdict(list)
    for h in hits:
        by_delta[tuple(sorted(h["delta"].items()))].append(h["id"])
    out: dict[str, str] = {}
    for key, ids in by_delta.items():
        opposite = tuple(sorted((k, -v) for k, v in key))
        for other in by_delta.get(opposite, []):
            for mid in ids:
                if mid != other:
                    out.setdefault(mid, other)
    return out


def _database_verdict(ctx: ReviewContext, row: dict, new_formula, new_charge,
                      is_mass: bool) -> tuple[str, str, str]:
    """Does any database confirm the balance-implied formula/charge?

    Returns (verdict, source, detail) where verdict is one of
    ``confirms`` / ``contradicts`` / ``silent``. Consults the cached MetaNetX
    reference first, then the live APIs when a client is configured —
    following the curator's rule that annotation propagates through the
    external identifiers and their APIs before anything is calculated.
    """
    if is_mass and row.get("ref_formula"):
        if norm_formula(new_formula) == norm_formula(row["ref_formula"]):
            return "confirms", "the MetaNetX reference", \
                   f"ref_formula={row['ref_formula']}"
        return "contradicts", "the MetaNetX reference", \
               (f"MetaNetX carries formula {row['ref_formula']}, not "
                f"{new_formula}")
    if not is_mass and row.get("ref_charge") is not None:
        if int(row["ref_charge"]) == new_charge:
            return "confirms", "the MetaNetX reference", \
                   f"ref_charge={row['ref_charge']}"

    if ctx.apis is None:
        return "silent", "", ""
    try:
        ev = ctx.api_evidence(row)
    except Exception as exc:
        return "silent", "", f"API lookup failed ({type(exc).__name__})"
    if not ev.get("sources"):
        return "silent", "", ""

    from apis import charge_is_supported              # noqa: PLC0415
    if is_mass:
        formulas = ev.get("formulas", {})
        for src, value in formulas.items():
            if norm_formula(value) == norm_formula(new_formula):
                return "confirms", src, f"{src} formula={value}"
        if formulas:
            src, value = next(iter(formulas.items()))
            return "contradicts", src, \
                   f"{src} carries formula {value}, not {new_formula}"
        return "silent", "", ""

    if charge_is_supported(ev, new_charge):
        return "confirms", "|".join(ev["sources"]), \
               f"charge {new_charge} is a recorded state"
    states = sorted({cs["charge"] for cs in ev.get("charge_states", [])
                     if cs.get("charge") is not None}
                    | set(ev.get("charges", {}).values()))
    if states:
        return "contradicts", "|".join(ev["sources"]), \
               (f"the databases record only charge state(s) {states} for this "
                f"species, and {new_charge} is not among them")
    return "silent", "", ""


def imbalance_residue(ctx: ReviewContext, localised: list[dict]) -> pl.DataFrame:
    """Reaction-level imbalance that no single metabolite explains."""
    model, _ = normalise_protons(ctx.model)
    explained = {h["id"] for h in localised}
    rows = []
    for r in model.reactions:
        if r.boundary:
            continue
        v = _imbalance(r)
        if not v:
            continue
        rows.append({
            "reaction": r.id, "name": r.name, "imbalance": str(v),
            "classes": ",".join(sorted(v)),
            "has_placeholder": any(has_placeholder(x) for x in r.metabolites),
            "localised_metabolites": ",".join(
                sorted(m.id for m in r.metabolites if m.id in explained)),
        })
    return pl.DataFrame(rows) if rows else pl.DataFrame()


# --------------------------------------------------------------------------
# CHECK 7 — network topology
# --------------------------------------------------------------------------

def check_topology(ctx: ReviewContext, f: Findings) -> None:
    for met in ctx.model.metabolites:
        reactions = list(met.reactions)
        if not reactions:
            f.add(met.id, "orphan_metabolite", "medium",
                  "metabolite participates in no reaction",
                  "remove it, or add the reaction that was intended",
                  check="topology")
            continue
        produced = consumed = False
        for r in reactions:
            c = r.metabolites[met]
            if (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0):
                produced = True
            if (c < 0 and r.upper_bound > 0) or (c > 0 and r.lower_bound < 0):
                consumed = True
        if not (produced and consumed):
            f.add(met.id, "dead_end_metabolite", "low",
                  "only " + ("produced" if produced else "consumed")
                  + " — no route in the other direction, so it cannot carry flux",
                  "add the missing reaction or a boundary exchange",
                  check="topology")


# --------------------------------------------------------------------------
# CHECK 8 — annotation completeness and schema
# --------------------------------------------------------------------------

def check_annotation_completeness(ctx: ReviewContext, f: Findings,
                                  core=("kegg.compound", "chebi",
                                        "metacyc.compound", "inchikey"),
                                  threshold: int = 3) -> None:
    """Metabolites too thinly annotated to be verifiable.

    Policy 6: a metabolite carrying *any* cross-reference outside the three
    namespaces a reconstruction pipeline can mint for itself (bigg,
    metanetx, seed) has been checked against something real and is not the
    curator's problem. Suppressing on that condition removes the bulk of
    round 1's low-value noise without losing a single scored true positive.

    Findings here are annotation work, not model errors, so they route to
    the annotation-update file (policy 5).
    """
    j = ctx.joined
    present = [c for c in core if c in j.columns]
    external = [c for c in ctx.independent_namespaces
                if c in j.columns and c not in ctx.pipeline_namespaces]
    external += [c for c in ("inchikey", "smiles", "inchi") if c in j.columns]
    for row in j.iter_rows(named=True):
        if any(row.get(c) for c in external):
            continue
        missing = [c for c in present if row.get(c) is None]
        if len(missing) < threshold:
            continue
        fillable = ([c for c in missing if row.get("ref_" + c) is not None]
                    if row.get("metanetx.chemical") else [])
        f.add(row[ctx.id_col], "missing_annotations", "low",
              f"missing {len(missing)}/{len(present)} core cross-references "
              f"({','.join(missing)}) and no cross-reference outside the "
              "namespaces a reconstruction pipeline mints for itself",
              ("fill " + ",".join(fillable) + " from the assigned MNXM"
               if fillable else "no reference values available"),
              check="annotation_completeness", route="annotation")


def check_schema(ctx: ReviewContext, f: Findings) -> None:
    """Misspelled or duplicated annotation columns silently lose data."""
    known = set(ctx.xref_columns) | {
        "id", "name", "formula", "charge", "inchikey", "inchi", "smiles",
        "bigg.metabolite", "metanetx.chemical", "sabiork.compound", "vmhm",
        "slm", "compartment", "base_id"}
    for c in ctx.mets.columns:
        if c in known or c.startswith(("valid_", "_")):
            continue
        near = [k for k in known
                if k != c and len(k) == len(c)
                and sum(a != b for a, b in zip(k, c)) <= 1]
        if near:
            n_filled = ctx.mets.height - ctx.mets[c].null_count()
            f.add("<schema>", "schema_anomaly", "high",
                  f"column '{c}' looks like a misspelling of '{near[0]}' and "
                  f"carries {n_filled} values that no check will read",
                  f"rename '{c}' to '{near[0]}' and merge the values",
                  evidence=c, check="schema")


# --------------------------------------------------------------------------
# runner
# --------------------------------------------------------------------------

def run_all(ctx: ReviewContext) -> tuple[pl.DataFrame, dict]:
    """Run every check and return (all findings, side outputs).

    Findings carry a ``route`` that decides which file they land in — the
    curator worklist holds only what a human can act on now; MetaNetX and
    annotation-completeness findings go to the annotation-update file, and
    balance findings are held out of round 1 entirely (policies 4 and 5).
    Split them with ``Findings.by_route`` / ``routed_frames``.
    """
    f = Findings()
    check_schema(ctx, f)
    check_metanetx_ids(ctx, f)
    check_reference_agreement(ctx, f)
    check_xref_agreement(ctx, f)
    check_bigg_namespace(ctx, f)
    check_inchikey_validity(ctx, f)          # policy 9 — before any key is used
    check_bigg_structure_conflict(ctx, f)
    check_duplicates(ctx, f)
    check_name_collisions(ctx, f)
    check_shadow_duplicates(ctx, f)
    check_reaction_anchored_duplicates(ctx, f)
    check_uncurated_name_stubs(ctx, f)
    check_fantasy_metabolites(ctx, f)
    check_charge_plausibility(ctx, f)
    check_reference_without_structure(ctx, f)
    check_cross_compartment(ctx, f)
    proton_fixed = check_proton_convention(ctx, f)
    check_placeholder_formulas(ctx, f)
    localised = check_balance_localisation(ctx, f, proton_fixed)
    residue = imbalance_residue(ctx, localised)
    check_topology(ctx, f)
    check_annotation_completeness(ctx, f)
    frame = f.frame()
    return frame, {"localised": localised, "imbalance_residue": residue,
                   "proton_fixed": proton_fixed,
                   "routes": routed_frames(frame),
                   "reaction_handoff": f.reaction_handoff(ctx.model)}


def routed_frames(frame: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """Split findings by route into the three deliverable files."""
    return {r: frame.filter(pl.col("route") == r) for r in ROUTES}
