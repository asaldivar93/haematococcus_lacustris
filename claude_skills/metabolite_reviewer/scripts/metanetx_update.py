"""Build the MetaNetX annotation file the curator reviews between rounds.

Most of what a metabolite review can detect depends on the metabolite being
correctly mapped to a reference entry: a wrong MNXM makes every formula,
charge and structure comparison downstream compare against the wrong
compound, and a missing one makes them impossible. So identifier work is not
reported as curator errors at all — it is collected here, the curator applies
it to the model, and the *second* review round runs against corrected
annotations.

The file has two sheets, because the two cases need different judgements:

``to_add``
    The metabolite carries no ``metanetx.chemical``. Candidates are proposed
    from whatever the species already offers, strongest evidence first:
    InChIKey (exact structure), then cross-references the curator already
    trusts (ChEBI, KEGG, MetaCyc, SEED), then name, and finally
    formula+charge. The evidence tier is written out beside the candidate so
    the curator can see *why* it was proposed rather than being asked to
    trust a bare identifier.

``to_review``
    The metabolite has an identifier that is deprecated, unrecognised, or
    disagrees with the model's own chemistry. The current value and the
    proposed replacement sit side by side with the disagreement spelled out.

Both sheets carry an ``accept`` column the curator fills, in the same
TRUE/blank convention as the review workbook.
"""
from __future__ import annotations

import polars as pl

#: Evidence tiers, strongest first. The tier is recorded on every proposal:
#: a structure match and a name match are not the same kind of claim, and the
#: curator should not have to infer which one they are looking at.
TIERS = ("inchikey", "xref", "name", "formula_charge")

XREF_PRIORITY = ("chebi", "kegg.compound", "metacyc.compound", "seed.compound")

#: Above this many formula+charge matches, the tier proposes nothing. Set from
#: the observed distribution: genuine one-or-two-way ambiguity is worth showing
#: a curator, a 90-way tie between lipid isomers is not.
FORMULA_CHARGE_CAP = 3


def _mnx_by_inchikey(mnx: pl.DataFrame) -> dict[str, list[str]]:
    if "inchikey" not in mnx.columns:
        return {}
    g = (mnx.filter(pl.col("inchikey").is_not_null())
            .group_by("inchikey").agg(pl.col("id").unique().alias("ids")))
    return dict(zip(g["inchikey"].to_list(), g["ids"].to_list()))


def _mnx_by_name(mnx: pl.DataFrame) -> dict[str, list[str]]:
    if "name" not in mnx.columns:
        return {}
    g = (mnx.filter(pl.col("name").is_not_null())
            .with_columns(pl.col("name").str.to_lowercase().str.strip_chars()
                          .alias("_n"))
            .group_by("_n").agg(pl.col("id").unique().alias("ids")))
    return dict(zip(g["_n"].to_list(), g["ids"].to_list()))


def _mnx_by_formula_charge(mnx: pl.DataFrame) -> dict[tuple, list[str]]:
    cols = {"formula", "charge"}
    if not cols <= set(mnx.columns):
        return {}
    g = (mnx.filter(pl.col("formula").is_not_null()
                    & pl.col("charge").is_not_null())
            .group_by(["formula", "charge"]).agg(pl.col("id").unique()
                                                 .alias("ids")))
    return {(f, c): ids for f, c, ids in
            zip(g["formula"].to_list(), g["charge"].to_list(),
                g["ids"].to_list())}


def propose_for_row(row: dict, idx: dict) -> tuple[list[str], str, str]:
    """(candidate MNXM ids, tier, evidence). Empty list when nothing matches."""
    ik = row.get("inchikey")
    if ik and ik in idx["inchikey"]:
        return idx["inchikey"][ik], "inchikey", f"inchikey={ik}"

    for col in XREF_PRIORITY:
        val = row.get(col)
        if val and idx["xref"].get(col, {}).get(str(val)):
            return idx["xref"][col][str(val)], "xref", f"{col}={val}"

    name = (row.get("name") or "").strip().lower()
    if name and name in idx["name"]:
        return idx["name"][name], "name", f"name={row.get('name')}"

    # Formula+charge is the weakest tier and is only a proposal when it is
    # nearly unique. On this model it matched a median of 93 MetaNetX entries
    # per species — every positional isomer of a lipid shares one formula, so
    # a list of 93 ids is not a candidate set, it is the curator doing the
    # search by hand. Above the cap the metabolite is reported as unresolved
    # with the count, which is the honest statement.
    f, c = row.get("formula"), row.get("charge")
    if f and c is not None:
        hit = idx["formula_charge"].get((f, int(c))) or []
        if 0 < len(hit) <= FORMULA_CHARGE_CAP:
            return hit, "formula_charge", f"formula={f} charge={int(c)}"
        if hit:
            return [], "formula_charge_ambiguous", (
                f"formula={f} charge={int(c)} matches {len(hit)} MetaNetX "
                "entries; needs a structure or a name to resolve")
    return [], "", ""


def build_metanetx_update(path, mets: pl.DataFrame, findings: pl.DataFrame,
                          mnx: pl.DataFrame, id_col: str = "id") -> dict:
    """Write the two-sheet MetaNetX update file. Returns a row-count report."""
    idx = {
        "inchikey": _mnx_by_inchikey(mnx),
        "name": _mnx_by_name(mnx),
        "formula_charge": _mnx_by_formula_charge(mnx),
        "xref": {},
    }
    for col in XREF_PRIORITY:
        if col in mnx.columns:
            g = (mnx.filter(pl.col(col).is_not_null())
                    .group_by(col).agg(pl.col("id").unique().alias("ids")))
            idx["xref"][col] = dict(zip(g[col].cast(pl.String).to_list(),
                                        g["ids"].to_list()))

    # ---- to_add: no identifier at all ------------------------------------
    add_rows = []
    for row in mets.filter(pl.col("metanetx.chemical").is_null()) \
                   .iter_rows(named=True):
        cands, tier, evidence = propose_for_row(row, idx)
        add_rows.append({
            "id": row[id_col], "name": row.get("name"),
            "formula": row.get("formula"), "charge": row.get("charge"),
            "inchikey": row.get("inchikey"),
            "proposed_metanetx": "|".join(cands[:5]),
            "n_candidates": len(cands),
            "evidence_tier": tier or "none",
            "evidence": evidence,
            "note": ("unambiguous" if len(cands) == 1 else
                     "ambiguous — pick one" if cands else
                     "formula+charge alone cannot resolve this; supply a "
                     "structure or a curated name"
                     if tier == "formula_charge_ambiguous" else
                     "no candidate; the species may not exist in MetaNetX"),
            "accept": "",
        })

    # ---- to_review: an identifier exists but is suspect -------------------
    suspect_tags = ("deprecated_metanetx_id", "unknown_metanetx_id",
                    "multiple_metanetx_id", "formula_mismatch_reference",
                    "charge_mismatch_reference", "inchikey_mismatch_reference",
                    "invalid_inchikey")
    sus = findings.filter(pl.col("error_tag").is_in(suspect_tags))
    met_lookup = {r[id_col]: r for r in mets.iter_rows(named=True)}
    review_rows = []
    for row in sus.iter_rows(named=True):
        met = met_lookup.get(row[id_col], {})
        review_rows.append({
            "id": row[id_col], "name": met.get("name"),
            "current_metanetx": met.get("metanetx.chemical"),
            "issue": row["error_tag"], "detail": row["rationale"],
            "proposed_fix": row.get("proposed_fix"),
            "evidence": row.get("evidence"),
            "confidence": row.get("confidence"),
            "accept": "",
        })

    # Sort by id: the rows are accumulated in check order, which is not
    # stable across runs, and a curator diffing round N against round N-1
    # should see only real changes.
    add = (pl.DataFrame(add_rows).sort("id") if add_rows else pl.DataFrame(
        {"id": pl.Series([], dtype=pl.String)}))
    rev = (pl.DataFrame(review_rows).sort("id") if review_rows else pl.DataFrame(
        {"id": pl.Series([], dtype=pl.String)}))

    import xlsxwriter as xw
    wb = xw.Workbook(str(path))
    for name, frame in (("to_add", add), ("to_review", rev)):
        frame.write_excel(workbook=wb, worksheet=name, autofit=True,
                          header_format={"bold": True, "bg_color": "#DDDDDD"},
                          freeze_panes=(1, 0))
        ws = wb.get_worksheet_by_name(name)
        i = frame.columns.index("accept") if "accept" in frame.columns else None
        if i is not None:
            ws.data_validation(1, i, max(frame.height + 1, 2), i,
                               {"validate": "list",
                                "source": ["TRUE", "FALSE"]})
    wb.close()

    by_tier = (add.group_by("evidence_tier").len()
               .sort(["len", "evidence_tier"], descending=[True, False])
               if add_rows else pl.DataFrame())
    return {"to_add": add.height, "to_review": rev.height,
            "by_tier": by_tier.to_dicts() if add_rows else []}
