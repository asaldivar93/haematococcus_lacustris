"""Build the curator review workbook from a findings frame.

Sheet layout
------------
review          one row per (metabolite, error_tag) — the worklist, sorted so
                high-confidence rows come first; `curator_decision` and
                `curator_note` are left blank for the human to fill
by_metabolite   one row per flagged metabolite with every tag it carries and
                the current annotation, so a curator can judge a species once
                instead of once per finding
summary         counts per tag x confidence
imbalance       reaction-level imbalance that no single metabolite explains
context         the reactions each high-confidence metabolite participates in
legend          what every tag means, what evidence it rests on, and how the
                confidence tier was assigned
"""

from __future__ import annotations

import polars as pl

CONF_ORDER = {"high": 0, "medium": 1, "low": 2}

TAG_LEGEND = [
    ("deprecated_metanetx_id", "identifier", "MetaNetX chem_depr",
     "high — the deprecation table is authoritative"),
    ("unknown_metanetx_id", "identifier", "absent from chem_prop and chem_depr",
     "high — the id resolves nowhere"),
    ("missing_metanetx_id", "identifier", "no annotation present",
     "high — nothing can be validated for this species"),
    ("multiple_metanetx_id", "identifier", "more than one id in one cell",
     "medium — often a legitimate ambiguity"),
    ("formula_mismatch_reference", "annotation", "MetaNetX chem_prop formula",
     "high — but note this check is blind when the annotation was derived "
     "from the same MNXM"),
    ("charge_mismatch_reference", "annotation", "MetaNetX chem_prop charge",
     "high — same caveat"),
    ("inchikey_mismatch_reference", "annotation", "MetaNetX chem_prop inchikey",
     "medium"),
    ("inchikey_skeleton_conflict_bigg", "structure",
     "BiGG's own InChIKey, which is independent of the MNXM assignment",
     "medium; low when the names agree, since that indicates a tautomer or "
     "anomer choice rather than a wrong species"),
    ("xref_disagreement", "annotation", "MetaNetX chem_xref",
     "medium — release skew produces false positives"),
    ("bigg_id_not_in_namespace", "identifier", "BiGG universal namespace",
     "low — organism-specific ids are legitimate in a reconstruction"),
    ("bigg_id_deprecated", "identifier", "BiGG old_bigg_ids",
     "high — the rename target is explicit"),
    ("duplicate_metabolite", "duplicate", "shared MNXM / InChIKey / skeleton",
     "high for a shared MNXM or full InChIKey; low for a shared skeleton, "
     "where isomer families are expected"),
    ("name_collision", "duplicate",
     "normalised name, independent of every identifier",
     "high when the formulas also match and reaction degree is asymmetric; "
     "medium when only the name matches"),
    ("shadow_duplicate", "duplicate",
     "unannotated + low reaction degree + formula twin of an annotated hub",
     "high — this combination has no benign explanation"),
    ("uncurated_name_stub", "annotation",
     "name equals the identifier with the compartment appended",
     "high — no chemical name contains its own compartment tag"),
    ("charge_implausible_for_formula", "annotation",
     "the charge distribution of every MetaNetX entry sharing this formula",
     "high when the assigned MNXM itself has no charge (so the value was "
     "hand-filled) and at least 10 reference entries agree; else medium"),
    ("reference_lacks_structure", "coverage",
     "assigned MNXM has no formula or charge in chem_prop",
     "low — not an error, but marks the row as unverified rather than clean"),
    ("cross_compartment_inconsistency", "consistency",
     "the same base id in other compartments",
     "high when formula or charge differ; medium for name or MNXM only"),
    ("mass_imbalance_localised", "stoichiometry",
     "one consistent integer delta across every reaction of the metabolite",
     "high when the implied correction matches the MetaNetX reference, or the "
     "metabolite appears in at least 4 reactions; else medium"),
    ("charge_imbalance_localised", "stoichiometry", "same, charge only",
     "same rubric"),
    ("proton_charge_convention", "convention", "BiGG/MetaNetX proton convention",
     "high — mechanical, and it must be fixed before balance means anything"),
    ("placeholder_formula", "convention", "R/X/Z/* in the formula",
     "low — expected for macromolecules and pseudo-metabolites"),
    ("orphan_metabolite", "topology", "no reaction references the species",
     "medium"),
    ("dead_end_metabolite", "topology", "produced or consumed but not both",
     "low — common and often intentional in a draft reconstruction"),
    ("missing_annotations", "coverage", "core xref columns are empty",
     "low — completeness, not correctness"),
    ("schema_anomaly", "schema", "column name one edit away from a known one",
     "high — a misspelled column silently hides its values from every check"),
]


# --------------------------------------------------------------------------
# Curator decision vocabulary
# --------------------------------------------------------------------------
# The curator writes one of these verbs in `curator_decision`; the row is
# actioned only when `fix_accepted` is TRUE. Each verb names the sheet of
# updates_to_model.xlsx it becomes, so a filled review sheet converts into a
# model update with no further interpretation. Anything outside this
# vocabulary is left for a human to read and is not written to the logbook.
#
# The argument goes in `curator_note` where the verb needs one — the target
# id for a merge or replace, the value for a set_*.
DECISION_VOCABULARY = (
    ("accept", "apply the proposed_fix as written",
     "annotation (or changes, if the fix is structural)"),
    ("set_formula", "curator_note holds the correct formula", "annotation"),
    ("set_charge", "curator_note holds the correct charge", "annotation"),
    ("set_name", "curator_note holds the correct name", "annotation"),
    ("set_annotation", "curator_note holds field=value", "annotation"),
    ("merge_into", "curator_note holds the surviving id; this species is "
                   "deleted and its reactions redirected", "changes"),
    ("replace_with", "curator_note holds the id that supersedes this one; "
                     "the new id is created if absent",
     "changes (+ new_metabolites)"),
    ("delete", "remove the species and every reaction that needs it",
     "changes"),
    ("keep_both", "the candidates are genuinely distinct entities; record "
                  "the decision so the pair is not re-flagged", "annotation"),
    ("defer_to_reaction_review", "the metabolite is correct and the "
                                 "surrounding reaction is not", "none"),
    ("reject", "not an error", "none"),
)

DECISION_VERBS = tuple(d for d, _, _ in DECISION_VOCABULARY)


def build_workbook(path, findings: pl.DataFrame, mets: pl.DataFrame,
                   model, imbalance_residue: pl.DataFrame | None = None,
                   id_col: str = "id") -> None:
    import xlsxwriter  # noqa: F401  (polars uses it as the writer engine)

    review = (findings
              .with_columns(pl.col("confidence").replace_strict(CONF_ORDER, default=9)
                            .alias("_o"),
                            pl.lit("").alias("fix_accepted"),
                            pl.lit("").alias("curator_decision"),
                            pl.lit("").alias("curator_note"))
              # Sort on the full row after the priority keys: ties on
              # (confidence, tag, id) are routine — one metabolite carries
              # several duplicate findings — and an unbroken tie leaves those
              # rows in check-execution order, so two runs of the same model
              # produce workbooks the curator cannot diff.
              .sort(["_o", "error_tag", id_col, "check", "rationale"])
              .drop("_o"))

    annotation_cols = [c for c in
                       ["name", "formula", "charge", "metanetx.chemical",
                        "bigg.metabolite", "kegg.compound", "chebi",
                        "metacyc.compound", "inchikey"] if c in mets.columns]
    per_met = (findings.group_by(id_col).agg(
        pl.col("error_tag").unique().sort().str.join("; ").alias("tags"),
        pl.col("confidence").replace_strict(CONF_ORDER, default=9).min()
          .alias("_best"),
        pl.len().alias("n_findings"))
        .with_columns(pl.col("_best")
                      .replace_strict({v: k for k, v in CONF_ORDER.items()},
                                      default="low")
                      .alias("max_confidence")).drop("_best")
        .join(mets.select([id_col] + annotation_cols), on=id_col, how="left")
        .sort(["max_confidence", "n_findings", id_col],
              descending=[False, True, False]))

    summary = (findings.group_by(["error_tag", "confidence"]).agg(pl.len().alias("n"))
               .with_columns(pl.col("confidence").replace_strict(CONF_ORDER, default=9)
                             .alias("_o"))
               .sort(["_o", "n", "error_tag"],
                     descending=[False, True, False]).drop("_o"))

    # The `context` and `imbalance` sheets were dropped at the curator's
    # request: reaction context is faster to read in the model itself, and
    # balance findings do not belong in round 1 at all (they are written to
    # their own file and actioned in round 2).

    legend = pl.DataFrame(
        [{"error_tag": t, "category": c, "evidence_source": e,
          "confidence_rubric": r} for t, c, e, r in TAG_LEGEND])

    decisions = pl.DataFrame(
        [{"curator_decision": d, "means": m, "writes_to": w}
         for d, m, w in DECISION_VOCABULARY])

    sheets = {"review": review, "by_metabolite": per_met, "summary": summary,
              "legend": legend, "decisions": decisions}

    import xlsxwriter as xw
    wb = xw.Workbook(str(path))
    for name, frame in sheets.items():
        frame.write_excel(workbook=wb, worksheet=name, autofit=True,
                          header_format={"bold": True, "bg_color": "#DDDDDD"},
                          freeze_panes=(1, 0))

    # Dropdowns on the two columns the curator fills. Typed free text in
    # `curator_decision` cannot be actioned automatically, and a typo would be
    # silently skipped when the logbook is written — so constrain it here
    # rather than validating after the fact.
    ws = wb.get_worksheet_by_name("review")
    cols = review.columns
    last = review.height + 1
    for col_name, options in (("fix_accepted", ["TRUE", "FALSE"]),
                              ("curator_decision", list(DECISION_VERBS))):
        i = cols.index(col_name)
        ws.data_validation(1, i, max(last, 2), i,
                           {"validate": "list", "source": options,
                            "input_title": col_name,
                            "input_message": "pick one"})
    wb.close()
