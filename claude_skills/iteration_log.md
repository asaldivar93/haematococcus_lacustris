# metabolite_reviewer — iteration log

Test model: `test_model/nies144.xlsx` + `nies144.xml` (1785 metabolites, 2222
non-boundary reactions, 6 compartments c/h/m/u/x/e). XML and XLSX verified
identical (same metabolite set, no formula/charge disagreement) — either is a
valid single source of truth.

References: MetaNetX 4.5 (`chem_prop` 1,495,668 rows, `chem_xref`, `chem_depr`
92,378 rows), BiGG universal namespace (9,090 metabolites). MetaCyc flat files
are NOT present on this machine — the draft skill's MetaCyc step is
unexecutable here.

## Round 0 — baseline (draft skill, verbatim, no added checks)

What the draft skill instructs, and what it actually yields:

| draft step | executable? | result |
|---|---|---|
| update deprecated MNXM ids | yes | 9 deprecated, 38 missing, 0 unknown of 1785 |
| review/fill external annotations | yes | 0 disagreements across 10 xref namespaces |
| check name/formula/charge | yes | **0 formula mismatches**, 6 charge (all H+), 6 name (all H+) |
| duplicates | partially — no method given | 11 groups by MNXM, 10 by full InChIKey, 43 by skeleton |
| MetaCyc cross-check | **no** | database absent |
| trace downstream effects | no method given | not done |
| emit curator workbook | no spec given | not done |

### The pivotal negative result

`formula` vs MetaNetX reference formula: **0 mismatches out of 1464 comparable
rows.** Charge: 6 mismatches, all of them the H+ species. InChIKey: 0
mismatches. Every xref namespace: 0 disagreements.

The model's `name`, `formula`, `charge`, `inchikey` and all xref columns were
*derived from* the assigned MNXM. An incorrectly assigned MNXM therefore drags
its whole annotation block with it and stays perfectly self-consistent. **Any
check that compares model annotation against the annotation of its own assigned
MNXM is structurally blind to identifier-assignment errors.**

This kills most of the draft skill's checking strategy and reframes the problem:
discriminating signal must come from sources *not* downstream of the assigned
MNXM —
1. the **BiGG identifier** (independent of the MNXM assignment),
2. **reaction stoichiometry** (mass/charge balance),
3. **cross-compartment** agreement of the same base id,
4. **duplicate structure** within a compartment.

### Ambiguities / gaps found in the draft skill

- No mass/charge balance step at all — the single most informative source.
- No duplicate-detection method (which key? within or across compartments?).
- No handling of the H+ charge convention; naive charge checking flags 6 H+
  species and naive balance checking then reports 1286/2222 unbalanced.
- No handling of placeholder elements (R, X, Z, `*`) which make balance
  meaningless — 9 photon pseudo-metabolites carry formula `Z`.
- MetaCyc step has no fallback when the database is absent.
- "confidence" tiers named but no rubric for assigning them.
- Output workbook has no column spec.
- `bigg.metabolite` xrefs in MetaNetX are lowercased and `m_`-prefixed;
  comparing raw gives 433 false disagreements vs 36 after normalisation.
- BiGG's own `metanetx.chemical` column is MetaNetX-3.x era: 875/1655 disagree
  with the model's 4.5-era ids purely from release skew. Unusable as a direct
  check; the **InChIKey skeleton** from BiGG is the usable independent signal
  (41 mismatches).

## Round 1 — deterministic check battery (`scripts/checks.py`)

20 checks, 26 error tags, three-tier confidence. **1491 findings, 130
high-confidence across 100 metabolites** (baseline: 421 findings, 90
actionable).

New checks that the draft skill had no equivalent of, and what each one is
for:

| check | independent evidence it uses | why it was needed |
|---|---|---|
| `check_balance_localisation` | reaction stoichiometry | a metabolite wrong by Δ makes *every* reaction it appears in imbalanced by `coefficient × Δ`; requiring one consistent integer Δ localises the error to that metabolite. 56 localised. |
| `check_proton_convention` | BiGG/MetaNetX H+ convention | must be normalised *before* balancing: naive balance reports 1286/2222 unbalanced, 6 of which is the real story |
| `check_name_collisions` | normalised name only | finds duplicates that share **no identifier** — invisible to MNXM/InChIKey grouping |
| `check_shadow_duplicates` | degree asymmetry + formula | unannotated 1-reaction species duplicating a well-connected hub |
| `check_uncurated_name_stubs` | name vs id form | `hmcit-L[m]` for `hmcit_L_m` means nobody ever named the species |
| `check_charge_plausibility` | charge distribution of *all* MetaNetX entries sharing the formula | covers the blind spot below |
| `check_reference_without_structure` | — | 283 assigned MNXMs carry **no formula/charge in chem_prop**, so the reference-agreement check silently skipped them: those rows were *unverified*, not clean |
| `check_bigg_structure_conflict` | BiGG's own InChIKey | not downstream of the MNXM assignment |
| `check_schema` | column-name edit distance | `keeg.glycan` typo hid 31 values from every check |

### False positives caught and suppressed during round 1

- `uncurated_name_stub` first matched `name == base_id` and flagged **165**
  metabolites (`NADPH` for `nadph_c` is legitimate). Narrowed to the
  `base_id[comp]` bracket form only → **3**, all genuine. No chemical name
  contains its own compartment tag.
- Formula-only duplicate grouping within a compartment: 40+ groups, mostly
  genuine isomers/anomers. Unusable without a name or degree tie-breaker —
  which is why `shadow_duplicate` requires *unannotated + low-degree + hub
  peer* jointly.
- BiGG's `metanetx.chemical` column is MetaNetX-3.x era → dropped as a
  detector, its InChIKey kept.

### Round 1 high-confidence breakdown (submitted for scoring)

| tag | n |
|---|---|
| duplicate_metabolite | 43 |
| missing_metanetx_id | 38 |
| deprecated_metanetx_id | 9 |
| charge_imbalance_localised | 7 |
| proton_charge_convention | 6 |
| charge_implausible_for_formula | 6 |
| bigg_id_deprecated | 5 |
| mass_imbalance_localised | 5 |
| uncurated_name_stub | 3 |
| cross_compartment_inconsistency | 3 |
| shadow_duplicate | 2 |
| name_collision | 2 |
| schema_anomaly | 1 |

**Status: paused awaiting curator scoring.** Files staged for review in
`claude_skills/review_output/`:
`high_confidence_for_scoring.csv` (130 rows, blank `verdict_TP_FP` column),
`nies144_metabolite_review.xlsx`, `iteration1_findings.csv`,
`baseline_findings.csv`, `imbalance_residue.csv`.

### Baseline output

`baseline_findings.csv` — 421 rows / 353 metabolites, but 331 of those are
low-value `missing_annotations`. Actionable core: 9 deprecated + 38 missing
MNXM + 43 duplicate-group members.

## Round 2 — async API layer

The check battery could not finish: the first full run was killed at 1526 s,
all of it inside `unichem_validate_inchikeys`. UniChem's compound-search POST
answers in ~13 s through the proxy and has no bulk endpoint, so validating a
model's 804 distinct InChIKeys serially is about three hours.

First attempt used a thread pool. That was the wrong tool and the curator
said so: `references/unichem.py` already had a working async bulk pattern.
Rewritten on that basis.

* `unichem.py` — the async engine is now generic. `bulk_fetch_json(specs)`
  runs any list of independent HTTP calls with bounded concurrency, a shared
  token bucket, retry/backoff honouring `Retry-After`, and progress logging.
  `run_async()` runs a coroutine from synchronous code and raises rather than
  deadlocking if a loop is already running. Two defects fixed in the existing
  bulk function: it returned a list of records while its annotation and
  docstring both promised `{inchikey: bool}`, and being `async def` it could
  not reach its own `asyncio.run` fallback — the blocking entry point is now
  `check_inchikeys_unichem_bulk_sync`.
* `apis.py` — `unichem_validate_inchikeys` consults the on-disk cache first
  and sends only the misses through the async engine, writing results back
  under the same cache keys the synchronous client uses.
  `prefetch_annotation_evidence` advances the whole batch one dependent stage
  at a time (name resolution, then ChEBI records, then conjugate
  neighbours), each stage issued concurrently; the per-metabolite assembly
  then runs against a warm cache, so parsing and source-priority logic stay
  in one place. `spec_*` builders sit beside each endpoint function and both
  the endpoint and the prefetch use them, because a duplicated URL drifts and
  a drifted prefetch fails silently — every request just goes out serially
  again.
* All thread-pool code removed from `checks.py` and `apis.py`.

Measured on nies144:

| stage | before | after |
|---|---|---|
| 804 InChIKeys, cold | ~3 h (est., killed at 25 min) | 150 s |
| `run_all`, cold | did not finish | 181 s |
| `run_all`, warm cache | — | 15 s |

Detection unchanged: all 12 scored cases still flagged. The worklist is
275 high/medium findings over 192 metabolites (round 1: 1491 findings, 130
high-confidence), the reduction coming from the routing and suppressor work.

## Round 4 — determinism and reference-degradation reporting

Found by running the finished skill twice on identical input rather than by
scoring against the spiked set.

- **Output was not reproducible.** `reaction_handoff`, `balance_round2` and
  both workbooks differed between runs: same rows, different order. Three
  causes — cobra's `met.reactions` is a set, the element-delta dict was
  rendered into prose with unordered keys, and four `.sort()` calls had ties
  they did not break. Every write point now sorts on enough columns to be
  total. Two runs now differ only in the timestamp inside the xlsx container.
- **A missing reference degraded detection silently.** Omitting `bigg_rxns`
  dropped shadow-duplicate detection: 407 -> 395 findings, machine-applicable
  fixes on the scored set 6 -> 3, and all 12 spikes still flagged at high
  confidence, so nothing about the output looked wrong. `review` now reports
  `degraded` — which references were absent and what detection each one cost.

Scored set after round 4: 12/12 flagged, 12/12 high confidence, 12/12
diagnosed, 6 carrying a machine-applicable fix (the other 6 are three merge
survivors, which take no action, and the two phosphatidylcholines and
`cacoa_c`, where no database supports an automatic proposal).

## Round 5 — packaging as a reusable skill, and the offline path

Published the review as `metabolite-review`. Loading it into a clean kernel and
re-running the referenced review reproduced the folder run exactly: 1751
findings, 407 worklist rows, 82 high confidence, all five CSVs identical, the
review sheet equal frame-for-frame, 12/12 scored cases high with 6 fix verbs.

Packaging exposed a defect the referenced runs could never show. Both skill.md
and the code claimed a missing reference "skips rather than fails", but
`ReviewContext` declared mnx/depr/bigg as required fields, so a review with no
references raised TypeError before any check ran. The offline mode was
documented but had never been executed.

Making the references optional turned up a second, worse problem: the first
working offline run emitted *more* findings than the fully-referenced one —
6773 findings, 2480 worklist, 292 high confidence, against 1751/407/82. An
empty BiGG set made every id look absent from the namespace, so
`bigg_id_not_in_namespace` fired 1785 times. The absent reference was being
read as a negative verdict instead of as no verdict, which is exactly the
failure rule 1 exists to prevent: a finding no reference supports.

Fixes: each reference join is now conditional and contributes all-null
placeholder columns so the joined schema is constant either way; every check
that consumes a reference returns early when it is absent; and the fallible
derived properties yield empty containers. Offline now gives 155 worklist rows
and finds 9 of the 12 spiked errors — the three it misses are the ones that
need a reference to see. The referenced run is unchanged and byte-identical to
before the patch, which is the check that matters: degrading gracefully must
not cost anything when the references are there.
