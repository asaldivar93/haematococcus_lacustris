---
name: metabolite_reviewer
description: "Review the metabolites of a genome-scale metabolic model: find duplicates, wrong identifiers, implausible chemistry and species that do not exist, and produce a curator worklist with tagged findings, confidence and machine-applicable fixes. Use when auditing or curating the metabolite set of a GSMM."
---

# metabolite_reviewer

Find metabolite-level errors in a genome-scale metabolic model, tag them, and
hand a human curator a worklist they can act on.

You are not the reviewer of record. You produce evidence and proposals; the
curator decides. Everything below follows from that: a finding you cannot
defend with a reference is not a finding, and a worklist padded with true
statements that require no action is a worklist the curator stops reading.

## How to run a review

```python
import sys; sys.path.insert(0, "<skill>/scripts")
from review import review

report = review(model, mets, out_dir="review_output", stem="mymodel",
                mnx=MNX, depr=DEPR, bigg=BIGG, bigg_rxns=BIGG_RXNS,
                apis=client, bigg_model_count=counts)
```

**Pass every reference you have.** A check whose reference is missing skips
rather than fails, so an offline or partially-referenced run still produces
findings — but a skipped check looks exactly like a clean result. `review`
prints what it had to skip and returns the same under `report["degraded"]`;
read it before trusting a low finding count. Measured on `nies144`: omitting
`bigg_rxns` alone removes shadow-duplicate detection and halves the
machine-applicable fixes (6 -> 3 on the scored set) while still reporting a
plausible-looking 395-finding worklist.

| reference | what is lost without it |
|---|---|
| `mnx` | reference agreement, identifier hygiene, charge plausibility |
| `depr` | deprecated-identifier detection |
| `bigg` | BiGG namespace and structure-conflict checks |
| `bigg_rxns` | shadow-duplicate detection |
| `bigg_model_count` | fantasy-metabolite detection |
| `apis` | every database-verdict check (offline mode) |

Passing nothing is allowed and still produces a review: on the spiked test
model, a reference-free run found 9 of the 12 known errors and wrote a
155-row worklist, against 12 of 12 and 407 rows fully referenced. A check
whose reference is absent stays silent rather than firing on everything — no
BiGG table means "cannot judge this id", not "no id is in BiGG".

**Outputs are byte-stable across runs.** Same model in, same files out — the
only difference between two runs is the timestamp Excel embeds. Round N is
therefore diffable against round N-1, which is the point: what changed is what
the curator did.

`model` is a cobrapy model, `mets` a polars frame of the metabolite table.
Every reference is optional — a check whose reference is missing skips rather
than failing, so a review with no network still produces the offline findings.
`references/databases.py` builds the reference frames from source files.

One round writes seven files:

| file | contents |
|---|---|
| `<stem>_metabolite_review.xlsx` | the curator worklist — this is the deliverable |
| `<stem>_metanetx_update.xlsx` | identifiers to add / to review, before the next round |
| `<stem>_annotation_updates.csv` | annotation gaps, not curator errors |
| `<stem>_balance_round2.csv` | balance findings, held out of round 1 |
| `<stem>_observations.csv` | true statements about the model that need no action |
| `<stem>_reaction_handoff.csv` | what the reaction_reviewer agent must look at |
| `<stem>_imbalance_residue.csv` | imbalance not attributable to any single metabolite |

The curator fills `fix_accepted` and `curator_decision` in the workbook, and
`updates_logbook.update_logbook` turns the filled sheet into rows in the shared
model-updates logbook.

## The rules that govern every finding

These come from curator scoring of earlier rounds. They are the difference
between a review that is useful and one that is confidently wrong.

1. **A database annotation outranks anything you calculated.** If MetaNetX
   gives a formula and your arithmetic disagrees, the arithmetic is the
   hypothesis, not the answer.
2. **A metabolite may hold several legitimate charge states.** Enumerate them
   from ChEBI conjugate acid/base relations before calling one wrong. Do not
   infer a charge from what would make a reaction balance.
3. **Balance findings are held out of round 1.** Mass and charge imbalance is
   a reaction property. Use it to *localise* a suspect metabolite, then put the
   hypothesis to the databases; if they contradict it, the metabolite is right
   and the reaction is wrong — propose nothing and flag the reaction.
4. **Identifier work is not curator error.** Missing, deprecated and
   unrecognised identifiers go to the MetaNetX file. The curator applies them,
   and round 2 runs against corrected annotations.
5. **A metabolite carrying any xref beyond bigg / metanetx / seed is annotated
   enough.** Those three are what a reconstruction pipeline can mint for
   itself; anything else means a human or a curated database has been here.
6. **Two species that are chemically identical but hold distinct ids in an
   independent namespace are distinct entities, not duplicates.** Report them
   as `chemically_identical_but_distinct_in_databases` and let the curator
   judge.
7. **Where two ids are both valid, prefer the one used by more published
   models.** BiGG membership first, then model count.
8. **Validate an InChIKey against UniChem before flagging anything that
   depends on it.** This runs before every structure-consuming check.
9. **Every finding a metabolite-level fix cannot resolve carries a note for the
   reaction_reviewer**, naming the reactions and the reason.
10. **The updates logbook is appended to, never overwritten.** Other agents
    write into the same file and it is the permanent record of every change to
    the model. If none exists, ask — do not create one.

## Routing: what reaches the curator

A finding's `route` decides its file. This is the single most important thing
the skill does, because detection is not the hard part — a naive battery finds
thousands of true things.

- `review` — the curator acts on it now.
- `annotation` — an identifier or annotation gap; fix before the next round.
- `balance` — held for round 2 (rule 3).
- `observation` — true, and not an action. Dead ends, placeholder formulas on
  macromolecules, reference entries lacking a structure. On the test model these
  were 808 of 1104 candidate worklist rows.

Before demoting a tag to `observation`, check that no confirmed error is
reachable *only* through it.

## Confidence

- `high` — an independent reference contradicts the model, or two independent
  signals agree. The curator should be able to accept it without looking
  anything up.
- `medium` — one signal, or a judgement call the curator must make.
- `low` — a pattern worth a glance, no claim attached.

Confidence is about evidence, not severity.

## Proposed fixes must be machine-readable

Prose in `proposed_fix` is for the human. The same fix goes in
`fix_verb` / `fix_arg` in the vocabulary the logbook speaks — `merge_into` +
surviving id, `set_charge` + value, `set_annotation` + `field=value`. When the
curator writes `accept`, the logbook applies the verb.

Never parse the prose to recover the fix. It was tried, it mistranslated
merges, and a wrong guess writes an unauthorised change into a permanent
record. A check that cannot express its fix structurally leaves the fields
blank and `accept` on that row is reported back to the curator.

## The checks

Run by `run_all` in this order; the order matters where noted.

**Identifiers and references**
- `check_schema` — column and dtype anomalies in the metabolite table.
- `check_metanetx_ids` — missing, multiple, deprecated, unrecognised. → annotation
- `check_reference_agreement` — formula / charge / InChIKey vs MetaNetX.
- `check_xref_agreement` — the model's xrefs disagree with each other.
- `check_bigg_namespace` — id deprecated, or outside BiGG entirely.
- `check_inchikey_validity` — UniChem, **before** any structure-consuming check.
- `check_bigg_structure_conflict` — the model's InChIKey skeleton contradicts
  BiGG's for the same id.

**Duplicates** — four detectors, because duplicates enter a model four ways.
- `check_duplicates` — species sharing an MNXM or an InChIKey in one
  compartment. Suppressed to `chemically_identical_but_distinct_in_databases`
  when they hold conflicting ids in an independent namespace (rule 6).
- `check_name_collisions` — same name, different id.
- `check_shadow_duplicates` — no shared identifier at all: matched on formula
  and network neighbourhood. This finds the most common real case, a species
  added ad-hoc with no annotation.
- `check_reaction_anchored_duplicates` — a non-BiGG reaction whose name matches
  a BiGG reaction and whose stoichiometry differs in exactly one species. This
  reaches duplicates that carry no metabolite-level signal whatsoever, and it
  recovered three cases a curator had judged undetectable.

**Existence**
- `check_fantasy_metabolites` — species no database knows. Two tiers, from a
  name-search miss combined with published-model breadth. Templated names
  (acyl chains, ACP conjugates, photons) are excluded from the name channel:
  a miss there means the databases do not index that naming style, not that
  the compound is invented.
- `check_uncurated_name_stubs` — placeholder names. Downgraded to a
  reaction-reviewer signal; every scored stub was reachable only via its
  reaction.

**Chemistry**
- `check_charge_plausibility`, `check_proton_convention`,
  `check_placeholder_formulas`, `check_cross_compartment`.
- `check_balance_localisation` — attribute reaction imbalance to a single
  metabolite, put the hypothesis to the databases, and resolve to confirm /
  contradict / silent. Equal-and-opposite imbalances that cancel across a pair
  go to the reaction reviewer with no fix proposed.

**Context**
- `check_topology` — orphans and dead ends. → observation
- `check_annotation_completeness` — subject to rule 5. → annotation

## Databases

Priority: MetaNetX, MetaCyc, UniChem by InChIKey, BiGG, then the rest.
Adopt Identifiers.org conventions for every identifier and tag.

MetaCyc needs a subscription — check for API access or flat files in
`~/databases/metacyc/` before relying on it; if neither is available, say so
rather than silently skipping it.

**All slow API calls are async.** `references/unichem.py` holds the engine:
bounded concurrency, token-bucket rate limiting, retry honouring `Retry-After`,
checkpointing, and progress logging. Use `bulk_fetch_json` for any new bulk
endpoint. Do not add thread pools — these calls are pure network wait, and
threads pay a stack and a lock handoff per request for nothing. UniChem's POST
endpoint answers in ~13 s and has no bulk form; async took the model's 804 keys
from ~3 hours to 150 s, and the on-disk cache makes re-runs free.

The concurrent path and the ordinary path must produce the same cache keys, or
a re-run silently re-issues everything serially. Endpoint URLs are defined once
in `spec_*` builders and both paths go through them.

## Notes

- Use cobrapy and networkx as needed.
- This skill is meant to grow. When a curator rejects a finding, the fix is
  usually a suppressor with a stated reason, not a raised threshold — write
  down *why* the case is legitimate, in the code, next to the suppressor.

## References

- `scripts/review.py` — single entry point; writes every deliverable.
- `scripts/checks.py` — the check battery, tag vocabulary, routing.
- `scripts/workbook.py` — curator workbook with decision dropdowns.
- `scripts/metanetx_update.py` — the between-rounds identifier file.
- `scripts/updates_logbook.py` — idempotent append to the shared logbook.
- `references/apis.py` — ChEBI / PubChem / UniChem client, caching, evidence.
- `references/unichem.py` — the async engine.
- `references/databases.py` — normalise MetaNetX and BiGG from source files.
- `references/model_update.py` — apply a curated excel back to a model.
- skills: *Metanetx Chem Mapping*, *Metanetx Metabolite Update*.
