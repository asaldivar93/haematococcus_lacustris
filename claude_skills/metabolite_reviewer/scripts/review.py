"""One call that runs a review round and writes every deliverable.

    from review import review
    report = review(model, mets, out_dir="review_output", **refs)

Having a single entry point matters more than it looks. The round produces
five files that must agree with each other — the curator worklist, the
MetaNetX update, the held-out balance findings, the annotation updates, and
the reaction-reviewer handoff — and they are all views of one findings frame
under different routes. Assembling them by hand in a notebook is how a round
ends up with a worklist from one run and a handoff from another.

What comes out:

  <stem>_metabolite_review.xlsx   curator worklist (round 1 findings only)
  <stem>_metanetx_update.xlsx     identifiers to add / to review
  <stem>_balance_round2.csv       balance findings, held for round 2
  <stem>_annotation_updates.csv   annotation-completeness findings
  <stem>_reaction_handoff.csv     what the reaction_reviewer must look at

Then the curator fills `fix_accepted` and `curator_decision` in the workbook,
and ``updates_logbook.update_logbook`` turns those into model updates.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from checks import ReviewContext, run_all
from metanetx_update import build_metanetx_update
from workbook import build_workbook


def review(model, mets: pl.DataFrame, out_dir="review_output",
           stem: str | None = None, id_col: str = "id",
           **context_kwargs) -> dict:
    """Run one review round and write its deliverables. Returns a report.

    ``context_kwargs`` are passed to :class:`ReviewContext` — the reference
    tables (``mnx``, ``depr``, ``bigg``, ``bigg_rxns``, ``bigg_model_count``)
    and the API client (``apis``). Everything the checks need is optional;
    a check whose reference is missing skips rather than failing, so a review
    with no network still produces the offline findings.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = stem or getattr(model, "id", None) or "model"

    ctx = ReviewContext(model=model, mets=mets, id_col=id_col,
                        **context_kwargs)

    # Report which references are absent. Checks skip rather than fail when a
    # reference is missing, which keeps an offline review working — but a
    # silent skip is indistinguishable from a clean result. Omitting
    # `bigg_rxns` alone drops the shadow-duplicate check and, with it, half
    # the machine-applicable fixes on a model we had already scored.
    degraded = {name: cost for name, cost in (
        ("mnx", "reference agreement, identifier hygiene, charge plausibility"),
        ("depr", "deprecated-identifier detection"),
        ("bigg", "BiGG namespace and structure-conflict checks"),
        ("bigg_rxns", "shadow-duplicate detection"),
        ("bigg_model_count", "fantasy-metabolite detection"),
        ("apis", "all database-verdict checks (offline mode)"),
    ) if context_kwargs.get(name) is None}
    for name, cost in degraded.items():
        print(f"[review] no {name!r} reference: {cost} will not run")

    findings, side = run_all(ctx)
    routes = side["routes"]

    paths = {}
    wb_path = out / f"{stem}_metabolite_review.xlsx"
    build_workbook(wb_path, routes["review"], mets, model, id_col=id_col)
    paths["workbook"] = wb_path

    mnx = context_kwargs.get("mnx")
    if mnx is not None:
        p = out / f"{stem}_metanetx_update.xlsx"
        mnx_report = build_metanetx_update(p, mets, findings, mnx,
                                           id_col=id_col)
        paths["metanetx_update"] = p
    else:
        mnx_report = None

    for key, frame, suffix in (
            ("balance", routes["balance"], "balance_round2"),
            ("observations", routes["observation"], "observations"),
            ("annotation", routes["annotation"], "annotation_updates"),
            ("reaction_handoff", side["reaction_handoff"], "reaction_handoff"),
            ("imbalance_residue", side["imbalance_residue"],
             "imbalance_residue")):
        if frame is not None and frame.height:
            p = out / f"{stem}_{suffix}.csv"
            frame.write_csv(p)
            paths[key] = p

    return {
        "paths": {k: str(v) for k, v in paths.items()},
        "degraded": degraded,
        "n_findings": findings.height,
        "n_review": routes["review"].height,
        "n_review_metabolites": routes["review"][id_col].n_unique(),
        "n_high_confidence": routes["review"].filter(
            pl.col("confidence") == "high").height,
        "n_balance_held": routes["balance"].height,
        "n_observations": routes["observation"].height,
        "n_annotation": routes["annotation"].height,
        "n_reaction_handoff": side["reaction_handoff"].height,
        "metanetx": mnx_report,
        "by_tag": routes["review"].group_by("error_tag").len()
                  .sort("len", descending=True).to_dicts(),
    }
