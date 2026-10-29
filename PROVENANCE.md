# Provenance of the frozen protocol

This repository is a curated copy of the project's development repository, prepared for publication.
Its history starts at the copy; the development repository is not public. This repository was created on GitHub at
2026-10-26T02:28Z, about four hours after the confirmatory analysis was aggregated (2026-10-25T22:40Z), so its own
history cannot show when the protocol was frozen. The freeze is documented by the commit identifiers below, by the
content checksums that let anyone verify that the files here are the frozen files, and by the job records of the
computing cluster described in the paper; an independently time-stamped archive of the development history is a
decision for the authors (see the paper, Appendix A).

## Freeze

- Protocol `docs/protocol/preregistration_v2.md` frozen at **2026-10-25T20:58:37Z**: freeze-target commit
  `e62592444170907104d8375473d40b41cc053097`, freeze commit `1c1faf6841eef2b732df552ee31622a81f422fbc`,
  tag `prereg-v2-frozen` (development repository).
- Before the freeze, no model had been fitted to the observed flood labels. The labels were used only as a
  label-only incidence table and as AOI-level flooded fractions to which the synthetic labels of the pilots
  and the power analysis were calibrated. Through the latter, the pilots showed that the regional
  composition of names (AOI means of the name features) tracks AOI-level flood rates; this association
  between real flood rates and real names, seen before the freeze, is why the primary estimand holds the
  composition fixed (see the paper's Sect. 2 and Appendix A).
- The freeze was recorded in version control and in the job records of the computing cluster, which ran
  the confirmatory analysis from a container image built after the freeze. **The protocol was not deposited
  with a third-party registry before the analysis**; the paper therefore describes the test as protocol-based.

## Verifying the frozen content

The files in this repository are byte-identical to their frozen versions:

| File | sha256 |
|---|---|
| `docs/protocol/preregistration_v2.md` (as of tag `prereg-v2-frozen`, unchanged since) | `4d9870ad888ebf88457cb0c8c1f2283264168fbab59ebad3db9599008bfd6dde` |
| `docs/protocol/v2/lexicon_hydro_v1.yaml` | `6c9999921956ba13f929b18ef921afcd0e05334293ef0989d9357c6a8d59c58d` |
| `docs/data/ems_inventory_20261025.json` | `f5c55eff3e470bce1fee7ae228641d1d3b21361f30b6d5ac9b87eef438028699` |

```bash
shasum -a 256 docs/protocol/preregistration_v2.md docs/protocol/v2/lexicon_hydro_v1.yaml \
  docs/data/ems_inventory_20261025.json
```

The digests of the other frozen inputs (GeoNames dump, coastline, incidence table, language-model
outputs, power-analysis runs) are listed in §10 of the protocol and in
`docs/verification/v33_freeze_record.json`.

## Changes after the freeze

- `docs/protocol/amendments/B-1_s6_group_buffer_fix.md`: an implementation defect in one sensitivity analysis
  (the location-group buffer did not reach the fold assignment) was found after its first result, corrected
  and re-run; the verdict did not change. Both runs are reported in the paper.
- `docs/protocol/amendments/B-2_posthoc_review_analyses.md` and its four addenda: post hoc analyses added after
  the confirmatory results in response to five rounds of critiques by language-model agents with access to the
  project records, which are not peer review (alternative yardsticks and a two-way
  bootstrap, recalibration-only reference models, grouping by flood event, reassigned folds, a label-free positive
  control, distance bands, results by class of name, within-area odds ratios, the mapping sensor, effect-size
  translation, power for the secondary splits). Most were planned, with their computations and reporting rules,
  before they were run; the paper marks those that were not, and those that the agents had computed first. All are labelled post hoc and change no pre-specified
  verdict. Scripts: `backend/verification/r_okauchi/v36`–`v40`; outputs: `docs/verification/`.
- `docs/protocol/amendments/B-3_gfm_independent_labels.md` and its two addenda, written in response to the fifth of
  those critiques and the authors' decision: independent flood labels from Copernicus Global Flood Monitoring
  (Sentinel-1), with the rules for making them, the comparisons (agreement of the
  label sets near settlements, flooding near water-side names, rescoring and refitting the primary comparison) and
  their interpretation fixed before the labels were made. Addendum 1 records that images with an unreadable layer
  were left out; addendum 2 is a check added after one result was seen and is reported as such. Scripts: `v41`,
  `v45`; label builder: `backend/packages/topsem/scripts/world_gfm.py`.
- `docs/protocol/amendments/B-4_lexicon_norway_yaml_fix.md`: an implementation defect found by a language-model
  agent that critiqued the manuscript and read this repository. The YAML parser read the country code of Norway, `NO`, as the Boolean false, so the lexicon was not
  applied to Norwegian names in the confirmatory analysis. The frozen lexicon was left unchanged, the parser was
  corrected and the confirmatory analysis was run again; the amendment fixed in advance that both runs are reported
  and that the frozen-code run stays primary if both reach the same verdict. Size of the defect: `v46`.
- `docs/protocol/amendments/B-5_blind_review_analyses.md` and its four addenda: post hoc analyses in response to
  critiques by language-model agents given only the manuscript, which could also read this repository; these
  critiques are not peer review (the distribution of the increment over 200 reassignments of names,
  second-stage coefficients, thirteen assignments of places to folds, bounds on missed effects by test inversion, an
  oracle control, within-area metrics, two variants of the labels, element groups of the lexicon, exclusion of the
  Americas and Oceania, and checks of the name readings against etymologies recorded on Wikidata and Wikipedia).
  Each was fixed before it was computed; the Wikipedia sample was enlarged before any excerpt was coded (addendum 4).
  Scripts: `v42`–`v44`, `v47`; the coded excerpts with sources and quotations: `docs/verification/v44_etym_codes.csv`.
- `docs/protocol/amendments/B-6_blind_review2_analyses.md` and its three addenda: post hoc analyses in response to a
  second round of such critiques: the largest increase in the pooled measure that perfect predictions of the ground
  around settlements could produce, bounds on missed effects recomputed by injecting effects into the observed labels
  (plasmode) with the odds ratio detectable with 80 % power, controls that carry flood information (the flooded share
  around each settlement in the same event, and surface water around each settlement from JRC Global Surface Water),
  diagnostics of the GFM labels with a variant that makes seasonal water assessable, documentation of the lexicon and
  a strict reading of it, further sensitivity runs (a cleaned gazetteer, a terrain ladder without river distances and
  water fractions, AOIs without a hydrography layer or a flood type left out, the language-model reading on the
  primary split, reassignments of names on the split grouped by activation), and the descriptions and p-values
  requested in addendum 1. Each was fixed before it was computed. Addenda 1 to 3 record the implementation events,
  none of which changed a computation: cluster jobs stopped before computing (an input file missing from the
  container image, because the Dockerfile did not copy it, which took three attempts to diagnose and fix, and a
  repository-root lookup), and the terrain-ladder arms without river
  distances were first skipped by a check that also required a feature built inside each fold, and were run again
  after the check was corrected. Scripts: `v48`–`v54`; the lexicon tables and the sample and guide for an independent
  coding of etymologies: `docs/verification/v51_*`.

- `docs/protocol/amendments/B-7_blind_review3_analyses.md` and its three addenda: post hoc analyses in response to a
  third round of such critiques: the permutation distribution rebuilt with the regional composition recomputed from
  the reassigned names, the plasmode with every model refitted on the injected labels, the activation-grouped split
  with the corrected code, a control for recurrent flooding from other EMS events (assigned to cells directly and
  through the nearest settlement), a floodplain reading of the lexicon (riverine and wetland classes), partly flooded
  cells, odds ratios by event magnitude and adjusted for built-up and tree cover, within-area discrimination on the
  GFM labels, a pairs bootstrap and joint tests of the near-settlement odds ratios, and descriptions (excluded AOIs,
  digraphs in ASCII names, recall by language, residual correlation between AOIs of one activation, Europe alone, the
  slope of the terrain ladder, fixed-threshold CSI). Each was fixed before it was computed; addendum 1 records the
  implementation choices of the analyses computed from existing artefacts after they were run, and addendum 2 the
  splitting of the permutations into seed ranges and the resubmission of jobs that had not started with lower resource
  requests; addendum 3 fixed, before it was computed, a joint test of the coefficients of the lexicon features in a
  model of the held-out flood labels given the held-out terrain prediction. By the rule fixed in B-7, the reading that
  real names cost less than reassigned ones (B-5) was withdrawn once the permutations were rebuilt. Scripts: `v55`,
  `v56`.
- Amendment B-8 followed the fourth round of critiques, made on the PDFs alone. It corrected the reporting (the abstract
  and the conclusions were brought into line with the reporting rules of B-2, B-3 and B-6) and fixed, before they were
  computed: the no-effect reference of the plasmode bounds (the rebuilt permutations of B-7), with Monte Carlo
  intervals and 40 more replicates (C1); the achieved equivalence bound; the location groups needed to settle the
  inconclusive comparisons; a permutation test; the second stage refitted inside the bootstrap (C2); the operating
  characteristics of the decision rule at the margin (C3); recurrent surface water of each cell (C4); a cleaned lexicon
  (C5); estimates near settlements by class and variant (the uncapped sample, a misclassification correction, river
  size, the two label sets on the same cells); heights of settlements above the main rivers; calibration and skill; an
  audit of matched forms; a reading by nearby river and lake names; and a second coder of the etymology excerpts (the
  default model of OpenAI's Codex). Addenda 1 and 2 record implementation choices before the results; addendum 3 the
  failure of the cluster's storage and the move of C1 to C5 to a local machine (runs marked `_local`); addendum 4 that
  item B10 was not counted, the tolerance of the reproduction check of C1 and details of B1 and C2 to C4; addendum 5
  the return of C1 and C3 to the cluster when the storage was back; addenda 6 and 7 two requests missing from the plan
  (estimates near settlements without Pakistan and Bangladesh, and the precision of the strict reading by language);
  addendum 8 the rerun on the cluster of one C1 job whose output failed to save; and addendum 9, written after the
  values had been seen, the withdrawal of the odds ratios read from the B-7 simulation with every model refitted,
  whose grid has no effect small enough to read them without bias (a straight line from zero across a convex gain
  curve understated them). Scripts: `v57`, `v57b`, `v58`, `v59`, `v60`, `v61`, `v62`.
- Amendment B-9 followed the fifth round of critiques. It corrected the reporting and fixed, before they were computed,
  bounds read with refitted intervals and over fold assignments, a finer grid of small injected effects, effects that
  remove flooding, effects that differ between settlements, simulated power with every model refitted, the Murphy
  decomposition, an evaluation within AOIs, the skill of flood maps, one definition of permanent water for all AOIs,
  other subsampling seeds and a bare-earth terrain model (FABDEM). Addendum 1 records the implementation before any
  cluster result. Scripts: `v63`, `v64`, `v65`; FABDEM tiles: `backend/packages/topsem/scripts/world_fabdem.py`.
- Amendment B-10 followed the sixth round of critiques. It reframed the paper (what the pooled test bounds and cannot
  detect), corrected the reported numbers, reported planned quantities that had been left out, and fixed, before they
  were computed: the islands-and-meadows group at settlement sites, the floodplain reading by country and form,
  detectable effects and jackknife and wild-bootstrap intervals near settlements, misreading over a grid of rates,
  validation by country, sensors and product generations, the activation-merged split in full, the within-AOI
  evaluation on GFM labels, the near-settlement estimator on injected labels, one definition of the HAND-only map,
  the between-area relation within language areas, cells read in another language, the class distribution of the
  language-model reading, distance and continent interactions, a penalty tuned by grouped cross-validation, country
  indicators, spatially clustered injection, and JRC flood hazard maps in the baseline. Addendum 1 records the
  implementation before any cluster result. Addendum 2, written after the within-AOI evaluation on GFM labels had
  returned an interval above zero, added the placebo comparison that the decision table requires for such an interval
  (and the same comparison on the EMS labels); the paper reports it as added after the result. Scripts: `v66`, `v67`;
  hazard tiles: `backend/packages/topsem/scripts/world_jrc.py`.
- Amendment B-11 followed the seventh round of critiques. It refocused the paper on the scheme that was tested (name
  readings assigned to the nearest 250 m cells) and the within-area estimand, moved the ceiling of local effects into
  the design, recorded that the interpretation plan was abandoned after the result, corrected the reporting, and fixed,
  before they were computed: the change of the activation-merged split between the frozen and the corrected code, by
  location group; odds ratios near settlements on every evaluable cell (no cap) on the EMS and GFM labels and on the
  AOIs whose EMS maps were not drawn from Sentinel-1; the bounds on the scale of the near-settlement estimator; a test
  of the GFM domain with the EMS-trained predictions; the denominator of the flooded fraction; a corrected
  repeated-cross-validation interval and an interval that combines refitting, two-way resampling and the denominator; a
  census of the exploratory interval estimates (`docs/verification/v68_estimates.csv`); the gazetteer's feature codes
  and rounded coordinates; the pilots; a new sample and guide for a toponymic check
  (`docs/verification/v68_audit_sample.csv`, `v68_audit_guide.md`; the study's readings of that sample are withheld
  until the check is done); and, on the cluster, reassignments with rebuilt composition on the country split. The
  uncapped near-settlement estimates were made the main ones after the capped and an uncapped estimate had both been
  seen. Script: `v68`; cluster run `world_v2_b11ctryperm`.
- Amendment B-12 followed the eighth round of critiques. It restructured the paper (the pre-specified test alone in
  Sect. 3; the uncertainty that the protocol's interval leaves out; what assignment through settlements can carry; the
  limits of the readings and of the labels near settlements), removed the odds ratios near settlements from the abstract,
  and fixed, before they were computed: the operating characteristics of the protocol's decision rule and of a rule on
  the interval of the ratio; the composition of the uncapped sample near settlements and an estimate without Pakistan and
  Bangladesh; a recount of the cells whose country's first language no lexicon block covers; the inclusion criteria in
  order; the outlying null simulations; highland names and GFM's domain; the recurrence control on the pooled scale;
  flexible single-reading comparisons; a second coder's intervals; the leak of the primary split and the reference model;
  a penalty on the names alone; the GFM time window; flood typing; one HAND-only map; constant forecasts; built-up land
  near settlements; a specification curve of every version of the odds ratio near settlements; and a revised audit
  sample of 300 names that mixes names read and not read as water-side (`docs/verification/v69_audit_sample.csv`,
  `v69_audit_guide.md`; whether the study read each name as water-side is withheld until the check is done), and, on the
  cluster, an injection within 1 km of settlements. Addendum 1 corrected the near-settlement flags of the census of B-11
  (an implementation error; `docs/verification/v69_estimates.csv`) and added Benjamini–Yekutieli values; addenda 2 and 3
  corrected the label flags and the versions taken by the specification curve. Tables that neither the paper nor the
  Supplement cites are in `docs/paper/paper_2_nhess/supplement_archive.pdf`. Script: `v69`; cluster run
  `world_v2_b12plas1k`.
- Amendment B-13 followed the ninth round of critiques. It restructured the paper into four parts (the protocol-based
  test with all of its uncertainty; designated checks; what assignment through settlements can carry; the validity of
  the readings and of the labels near settlements), kept the protocol's verdict as the record, and fixed, before they
  were computed: an equivalence test on the ratio that carries the uncertainty of the baseline gain (the Fieller-type
  limits of the change in Brier score minus the margin times the baseline gain, from the same resamples), with the
  models fixed, with the second stage refitted, with activations resampled within countries and with the spread over
  fold assignments, and its operating characteristics; the odds ratios that each interval's upper limit corresponds to;
  the second coder's rates in the dilution examples; size strata near settlements; the reference model without the
  character model's area mean; influence; other denominators; a joint test of the name coefficients; further checks of
  the odds ratios near settlements; an audit of homographs, language blocks and Unicode names; the language-model
  reading by block and country; both coders' codes (`docs/verification/v70_two_coders.csv`, leaving out the names of
  every audit sample); AOI-balanced metrics; a diagnosis of the GFM labels; pairs of settlements that differ by a lower
  or upper qualifier; and, on the cluster, records known before each event in place of the names (surface water on the
  AOIs of events from 2022, flooding in earlier events only, the JRC hazard maps through settlements), a lexicon with
  the errors named in the critiques corrected (`docs/verification/v70_lexicon_corrected.yaml`, derived mechanically from
  the frozen lexicon, which is unchanged), and early stopping on held-out location groups. Addenda 1 to 5 record the
  implementation before the results, including one correction (the run on events from 2022 keeps the cells of the main
  table). Scripts: `v70`, `v70b`, `v70_lexicon_corrected`, `v71`, `v72`; cluster runs `world_v2_b13gsw22`,
  `world_v2_b13recprior`, `world_v2_b13jrcset`, `world_v2_b13lexfix`, `world_v2_b13gres`.
- Amendment B-14 followed the tenth round of critiques. It kept the protocol's verdict as the record, split Sect. 3 into
  the protocol-based test and a post hoc reanalysis whose result is the one-sided upper limit (non-superiority), a choice
  made after every interval had been seen, and fixed, before they were computed: the equivalence test on the ratio in the
  order its intervals were defined, with percentile and normal limits, TOST p values and Monte Carlo errors; jackknife,
  studentized, weighted and pigeonhole versions and the effective number of clusters; operating characteristics with
  activations resampled within countries; the decomposition of the reference model and a comparison without the regional
  composition; the share of each record's value kept through settlements; the direction of the name coefficients on
  both label sets; tables of the lexicon's edits (`docs/verification/v73_lexicon_edits.csv`), of the elements it leaves
  out and of the gazetteer; elevation and river size of qualifier pairs; further checks of the odds ratios near
  settlements; a check of the positional join of the cell tables; and, on the cluster, labels on observable ground, the
  relief of settlement sites assigned through settlements as a positive control, and country indicators for the
  language-model reading. Addenda 1 to 3 record the implementation and the aggregation before the results. Addendum 4,
  written after the pigeonhole result had been seen, adds descriptions of the overlap of location groups and activations
  and reports both upper limits. Addendum 5 fixed versions of the odds ratio near settlements without EMSR052 and EMSR130
  before they were computed; addendum 6 records that its premise was wrong (amendment B-13 had computed these versions,
  which v75 reproduced) and that the only new fact is that neither GFM table holds an AOI of either activation. Scripts:
  `v73`, `v74`, `v75`; cluster runs `world_v2_b14obs`, `world_v2_b14oracle`, `world_v2_b14h4k`.

English translations of the protocol and of amendments B-1 to B-14 are in `docs/protocol/en/`; they were drafted with a
language model after the analysis and have not yet been checked by the authors, and the Japanese files are the record
and govern where the two differ. Amendments B-2, B-3 and B-5 to B-14 followed critiques by language-model agents; the
paper lists the rounds and their times. Dates written
inside the repository files are local (Japan Standard Time, UTC+9) unless marked; the paper gives UTC. An interpretation
plan, written after the freeze and before any confirmatory result, is in
`docs/research/2026-10-26_v2_interpretation_plan.md`.

The GFM labels are built from the images listed in the public STAC catalogue of EODC at the time of retrieval.
Rebuilding them on 2026-10-27 for the variant of amendment B-6 reproduced the labels of amendment B-3 exactly in 189 of
192 AOIs; in the other three, 23 images could not be read on the rebuild, which changed 16 of 6.76 million cells.
An earlier version of the example figure (Fig. 3) marked GeoNames features of every class within the AOI, not only the
populated places that the analysis uses; the figure was corrected, and no analysis was affected.

The confirmatory run records label the protocol "DRAFT until frozen" because the code wrote a fixed string
there; the status string is now read from the protocol's status line. The power-analysis runs record the
lexicon under its draft version label (sha256 prefix `4f9203535d64`); its entries are identical to the frozen
lexicon, which differs only in a status comment and the version label.
