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
