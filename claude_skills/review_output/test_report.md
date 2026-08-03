# metabolite_reviewer — test report

Detection of the spiked errors in `nies144`, pass 1 (draft skill, run verbatim)
versus pass N (polished skill). The 12 scored cases are the spiked errors the
curator confirmed after scoring round 1; the agent did not know them in advance.

![detection progress]({{artifact:art_651028fe-2c96-4e0c-b7e6-fec7d69bfb7c}})

## Headline

| measure | pass 1 (draft) | pass N (polished) |
|---|---|---|
| spiked errors **flagged** | 10 / 12 | **12 / 12** |
| flagged at **high** confidence | — (no confidence model) | **12 / 12** |
| **diagnosed** (specific defect named) | 0 / 12 | **12 / 12** |
| carries a **machine-applicable fix** | 0 / 12 | 6 / 12 (+6 survivors, nothing to apply) |
| worklist findings the curator must read | 128 | 407 |
| metabolites on the worklist | 138 | 275 |
| wall clock, warm cache | n/a (no API layer) | 17 s |
| outputs byte-stable across runs | not applicable | yes |

Pass 1 is generous to the draft: it counts every metabolite the draft's four
executable steps touch (deprecated/missing MetaNetX id, and duplicate groups by
MetaNetX id, full InChIKey, or InChIKey skeleton), whether or not the draft says
what to do about it.

## The distinction that matters: flag vs. diagnose

Pass 1's 10 hits are not 10 answers. Six of them are three duplicate *pairs*, and
the draft reports the pair — it never says which member is the error. Four more
are "no MetaNetX id", which is true of 38 species in the model and says nothing
about which of those are wrong. A curator handed the pass-1 output still has to
do the whole diagnosis by hand.

Pass N names the defect and, where the fix is unambiguous, writes it in the
decision vocabulary the updates logbook already speaks:

| metabolite | pass 1 signal | pass N high-confidence tag | pass N fix |
|---|---|---|---|
| `hicit_m` | inchikey + metanetx.chemical + skeleton | duplicate_metabolite | `— (merge survivor)` |
| `hisocit_m` | inchikey + metanetx.chemical + skeleton | duplicate_metabolite | `merge_into hicit_m` |
| `stcoa_c` | inchikey + metanetx.chemical + skeleton | duplicate_metabolite | `— (merge survivor)` |
| `ocdccoa_c` | inchikey + metanetx.chemical + skeleton | duplicate_metabolite | `merge_into stcoa_c` |
| `ocdcea_c` | inchikey + metanetx.chemical + skeleton | duplicate_metabolite | `— (merge survivor)` |
| `ocdce9a_c` | inchikey + metanetx.chemical + skeleton | duplicate_metabolite | `merge_into ocdcea_c` |
| `cacoa_c` | — | probable_fantasy_metabolite | curator decides |
| `hmcit_L_m` | missing_mnx | probable_fantasy_metabolite|shadow_duplicate|uncurated_metabolite | `replace_with hcit_m` |
| `gal1p_L_c` | — | probable_fantasy_metabolite|uncurated_metabolite | `merge_into gal1p_c` |
| `2oglutm_h` | missing_mnx | probable_fantasy_metabolite|shadow_duplicate|uncurated_metabolite | `merge_into 2ogm_h` |
| `pchol1801835z9z12z_c` | missing_mnx | uncurated_metabolite | curator decides |
| `pchol1801845z9z12z15z_c` | missing_mnx | uncurated_metabolite | curator decides |

The three merge survivors (`hicit_m`, `stcoa_c`, `ocdcea_c`) carry no verb by
design — the fix is applied to the species that disappears, and the survivor is
the argument of that fix. The two phosphatidylcholines are flagged as uncurated
with no automatic proposal: no database indexes that acyl-chain naming style, so
proposing a mapping would be a guess.

## The two cases pass 1 could not reach

`cacoa_c` and `gal1p_L_c` carry a complete, self-consistent annotation block —
formula, charge, InChIKey, cross-references — all of it derived from the assigned
MetaNetX id. Every pass-1 check compares the model against that same assignment,
so a wrong assignment is invisible: it agrees with itself. This was the pivotal
negative result of the baseline run — **0 formula mismatches out of 1464
comparable rows**. Detection had to come from evidence not downstream of the
assignment:

- `cacoa_c` — flagged by BiGG-namespace breadth (present in BiGG but used by a
  single published model) combined with a failed name lookup in both ChEBI and
  PubChem. Neither signal alone is sufficient; the conjunction is.
- `gal1p_L_c` — flagged by name/id form (`gal1p-L[c]` style stub, meaning nobody
  ever named the species) plus a well-connected `gal1p_c` hub carrying the same
  chemistry.

## Why the worklist grew and then shrank

Pass 2 (the check battery) produced 1491 findings — real detection, unusable
volume. Two subsequent passes cut it without losing a single scored case:

1. **Curator rules** removed whole finding classes from the worklist rather than
   downweighting them: identifier hygiene became annotation work, balance
   findings were held for round 2, and chemically-identical species holding
   distinct ids in an independent namespace stopped counting as duplicates.
2. **Routing** gave four tags that describe properties of the model rather than
   mistakes in it — dead ends, placeholder formulas on macromolecules, reference
   entries lacking structure, orphans — their own observations file.

Before demoting anything, every scored case was checked for reachability at high
confidence only. None was reachable only at low confidence, so the demotion is
safe against the known truth set.

## What the report cannot claim

- **Precision is unmeasured.** 407 findings were scored against 12 known spikes;
  the rest are unlabelled. A finding not in the truth set may be a real defect in
  the underlying reconstruction or a false positive, and this test cannot tell
  them apart. Only a curator pass over the full worklist can.
- **One model, one spiking pass.** The suppressors were tuned against these 12
  cases. Two of them — the BiGG-breadth threshold for probable-fantasy species
  and the formula+charge candidate cap in the MetaNetX file — are thresholds
  fitted on a single example each, and should be re-checked on the next spiked
  model before being trusted.
- **Detection plateaued at pass 2 on the flag measure** (12/12) and every later
  pass moved the *quality* measures — diagnosis, fix applicability, worklist
  size. A further pass would need a new truth set to move anything.


## Two defects this test itself surfaced

Neither was in the spiked set; both were found by running the finished skill
twice and comparing.

1. **Non-deterministic output.** Five of the six output files differed between
   two runs on identical input — same rows, different order — because three
   code paths iterated cobra sets and a dict of element deltas. The curator
   diffs round N against round N-1, so ordering noise is indistinguishable from
   real change. Fixed by sorting at every write point; two runs now differ only
   in the timestamp Excel embeds inside the file.
2. **Silent degradation on a missing reference.** Omitting `bigg_rxns` from the
   call disabled shadow-duplicate detection and cut machine-applicable fixes on
   the scored set from 6 to 3 — while still reporting a plausible 395-finding
   worklist with all 12 spikes flagged at high confidence. Nothing warned. The
   skip-not-fail policy is right for offline use, but it needs to be visible:
   `review` now prints and returns `degraded`, naming each absent reference and
   the detection it costs.

The second is the more serious of the two. A run that looks clean because a
reference was missing is worse than a run that fails.

## Reproducing

```python
from review import review
report = review(model, mets,                 # cobra model + polars metabolite table
                out_dir="review_output",
                mnx=mnx, depr=depr, bigg=bigg,
                bigg_model_count=bigg_model_count, apis=apis)
```

Every reference argument is optional: a check whose reference is missing skips
rather than failing, so a review with no network still produces the offline
findings. Measured warm re-run on `nies144`: **17 s** for 1751 findings, of which
407 reach the curator worklist across 275 metabolites. The warm path re-reads the
on-disk API cache and issues no network calls, so re-runs after a code change are
effectively free; the first run on a new model pays the API cost once.

## Post-hoc: the offline path

Packaging the review as a loadable skill forced a code path the scored runs
never touched — a review with no reference tables at all. Both the skill
document and the code comments said a check with a missing reference skips
rather than fails; the entry point in fact raised `TypeError` before any check
ran, because three reference frames were declared required. The documented
offline mode had never been executed.

Making them optional surfaced the more serious half. The first offline run
produced **more** findings than the fully-referenced one — 2480 worklist rows
at 292 high confidence, against 407 and 82 — because an empty BiGG table made
every identifier in the model look absent from the namespace. A missing
reference was being read as a negative verdict rather than as no verdict, which
is the precise failure the skill's first rule exists to prevent.

| run | findings | worklist | high confidence | spiked errors found |
|---|---|---|---|---|
| fully referenced | 1751 | 407 | 82 | 12 / 12 |
| offline, before fix | 6773 | 2480 | 292 | n/a (unusable) |
| offline, after fix | 2663 | 155 | 22 | 9 / 12 |

Every check that consumes a reference now returns early when it is absent, and
the reference joins contribute all-null placeholder columns so the joined
schema does not depend on which tables were supplied. The three spiked errors
lost offline are the ones that genuinely need a reference to see. The
referenced run is byte-identical to before the change — graceful degradation
cost nothing when the references are present.
