# B-7 Post hoc analyses in response to the third round of LLM-agent critiques given only the manuscript (NHESS version) (2026-10-27)

English translation of `docs/protocol/amendments/B-7_blind_review3_analyses.md`; the Japanese file is the record and
governs where the two differ.

Status: **planned (2026-10-27, before any of the analyses below was computed)**. The confirmatory verdicts are not
changed. Everything below is **post hoc**. It responds to an editor and three referees — language-model agents that
read only the PDFs of the manuscript and its Supplement for the first time (reports in
`docs/research/reviews/2026-10-27_paper2_blind3_*.md`, commit c8db9c6). The paper describes these as post hoc analyses
fixed in this document after the confirmatory results were seen and before they were computed. The word "registered"
is used only for registration with a third party.
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924). Code: with the correction of B-4
(the corrected code). Comparisons are set beside the corrected run (world_v2_lexfix).

Values seen when this document was written: all confirmatory and B-2 to B-6 results (those in the paper version of
commit b31dc7b) and the values the reports pointed out or estimated (for example the synthetic-label translation of
B-2, OR 4.12; the robust standard deviation of the null replicates, 0.00210; R = 0.004 on the activation-grouped split
with the corrected code). None of the quantities below has been computed.

## Analyses (the numbers are used in reporting)

### 1. Permutation distribution rebuilt (composition recomputed from the reassigned names)

The permutations of B-5 item 1 (world_v2_b5perm) combined the reassigned name features with the regional composition
(c_*) of the **real** names. Because names move between settlements nearest to different numbers of cells, the AOI
means of the reassigned features differ from c_*, so every reassignment carries an AOI-level component that can put it
at a disadvantage against the real names (referee 2).

- **Permutations with recomputed composition**: in each fold of the primary split, for each of the reassignments of
  B-5 item 1 (seeds 1–200, names reassigned among the settlements of the same AOI), form c(q), the AOI means of the
  reassigned lexicon features and of the N-learn score refitted on the reassigned names, and fit the reference
  M1 + c(q) and the names model M1 + c(q) + q in the second stage (base models fitted once, as in the confirmatory run).
  Increment per reassignment ΔBS_j = BS(M1 + c(q)) − BS(M1 + c(q) + q). The real increment is M1c → M3n (for real names
  c(real) = c_*, so it has the same form).
- **Permutations of AOI-centred features**: for each reassignment, centre the reassigned features within each AOI
  (q − AOI mean) and add them to M1s (recalibration of M1 only). The real counterpart is M1s → M3d (M3d of B-2).
- **Reporting**: the mean and 95 % range (on the R scale) of both distributions, and the percentile of the real
  increment. Also the percentile of the confirmatory placebo (M3np) in the B-5 distribution (referee 2, minor point 10).
- **Rule**: unless the real increment lies at or above the 95th percentile of the distribution with recomputed
  composition, the paper withdraws its reading that real names cost less than reassigned ones (the reading of B-5) and
  states that the B-5 percentile was an artefact of how the reassignments were built. If it lies at or above the 95th
  percentile, the reading stays, with the values of the rebuilt distribution.

### 2. Plasmode with refitted base models

Item 2 of B-6 fitted the base models (M0, M1) once, before injection. The injected effect is then orthogonal to the
terrain prediction and cannot be absorbed by the terrain, position or settlement-layout features, so the gain may be
inflated and the odds-ratio translations too low (referee 2).

- **Computation**: form A (all cells nearest a water-side name), δ ∈ {0.5, 1.0, 1.5} (OR 1.65, 2.72, 4.48),
  repetitions 0–2. The injected labels are those of B-6 item 2 (same seeds and the same injection probabilities, i.e.
  the M1 prediction of the uninjected run). Refit everything on the injected labels — N-learn, M0, M1 (cross-fitted)
  and the second stage (M1c, M3n) — and record Brier scores on the injected held-out labels.
- **Reporting**: for the same δ and repetitions, the gain with fixed base models and with refitted base models, and
  their ratio. From the refitted gain curve (δ = 0 and the three points, linear interpolation in δ), the odds ratio that
  the margin corresponds to and the upper bound (mean increment of reassigned names as reference, U90); beyond the grid
  "> 4.48". Also the yardstick (M0 → M1) in the injected world.
- **Rule**: the abstract and text give the translation of the margin as the **range** over three methods — the
  synthetic-label model of B-2 (OR 4.12), the plasmode with fixed base models of B-6 (2.34) and this plasmode with
  refitted base models — and do not pick one value.

### 3. Activation-grouped split with the corrected code

The activation-grouped split of B-2 (world_v2_actgrp) was run with the frozen code. Run the same split with the
corrected code (the confirmatory arms and second stage) and set its verdict, R and intervals beside those of the
frozen-code run. Use the change in the Brier scores of M1c and M3n to describe why the correction of the Norwegian
reading moved R so much on this split.

### 4. Control for recurrent flooding (not using the labels of the event itself)

- **Quantity**: for each analysed cell, among the other EMS activations that evaluated it (events more than 30 days
  from the date of the target event; all labelled AOIs; the confirmatory label and evaluability rules), the share in
  which the cell was flooded (recurrence), and the number of such events. Cell locations are matched to each AOI's grid
  by longitude and latitude.
- **(a) Assigned to each cell**: features = [recurrence (0 where there is no other event), whether there is another
  event]. Second stage M1cEc → M3oEc (regional composition = AOI means in both).
- **(b) Assigned through the nearest settlement**: for each settlement, the mean recurrence of the cells within 500 m of
  its point that have other events, assigned to cells through the nearest settlement as names are (0 and an indicator
  where there is none). Second stage M1cE → M3oE.
- **Reporting**: R, 90 % interval and verdict of (a) and (b) on the primary split. R restricted to AOIs with other
  events (yardstick of that subset, cluster bootstrap over places), descriptive. Coverage (share of cells with other
  events, number of AOIs).
- **Reading** (as proposed by referee 1, fixed before the results): if (a) adds much and (b) little, the assignment is
  the bottleneck; if neither adds anything, a design that averages over whole areas cannot see recurrent flooding; if (b)
  exceeds the margin, recurrence information is registered even when assigned through settlements.

### 5. A reading of names that denote floodplain ground

A reading that uses only two classes of the frozen lexicon, riverine (river islands, river meadows, floodplains, banks,
dikes, polders) and wetland (the lexicon is unchanged; this is a subset of its classes). Toponymic rationale: as the
Introduction states, names of water bodies, crossings and lowlands can denote dry sites or relative position, whereas
riverine and wetland elements denote floodplain ground itself.

- Features: indicators of a riverine or wetland element in the nearest name, distance-weighted counts over the three
  nearest names, and an indicator of either (built as for the strict reading); placebo versions too.
- M1cF → M3nF on the primary split (placebo M3nFp), verdict.
- Odds ratios near settlements (within 500 m, given terrain, EMS and GFM) and the difference in HAND at the
  settlement (as in the positive control of B-2).

### 6. Partially flooded cells

- **(a) Distribution**: the flood cover (f_flood_any) of cells with y = 0, overall and within 500 m of settlements.
- **(b) Re-run without them**: drop cells with y = 0 and 0.05 < f_flood_any (partially flooded) and run the confirmatory
  comparison on the primary split (M1c → M3n, placebo), verdict.
- **(c) Odds ratios**: the EMS odds ratios near settlements (i) without the partially flooded cells and (ii) with the
  riverine flood cover f_riverine as a fractional outcome (quasi-binomial logistic regression with AOI intercepts and
  the logit of M1).

### 7. Odds ratios near settlements by event magnitude

Magnitude per AOI = share of flooded cells / share of cells with HAND ≤ 2 m (hand_min). Split AOIs into terciles and
estimate the EMS odds ratio within 500 m (given terrain) in each.

### 8. Adjustment for built-up and tree cover

Recompute the difference in HAND at settlements (positive control) and the odds ratios near settlements (EMS, GFM)
with the CCI built-up and tree fractions added as covariates (to address the upward bias of HAND from a surface model
under buildings and canopy; a rebuild on a bare-earth elevation model is not done in this revision because of the data
volume).

### 9. Within-AOI discrimination on the GFM labels

On the GFM-assessable cells and GFM labels, the within-AOI AUC of the held-out predictions of M1, M1c and M3n, and the
CSI when, in each AOI, as many cells as were observed flooded are marked flooded in order of prediction. Median and
interquartile range of M3n − M1c.

### 10. Intervals of the odds ratios near settlements

Report a pairs cluster bootstrap over places (2000 replicates) as the main interval for the EMS odds ratios (water-side,
by class, by distance band), because there are few places; the existing cluster-robust values are kept beside them.

### 11. Descriptions from existing artefacts

- (a) The effect-size translations of B-2 (`docs/verification/v37_positive_control.json`, effect_size_translation:
  margin, upper limit, 500 m, range over fold assignments, activation split) as a Supplement table.
- (b) The robust diagnostics of the interval coverage of the null replicates (v48: robust standard deviation, standard
  deviation of standardised deviations) in the text, with the character of the outlying replicates (their yardsticks).
- (c) A table of the Wikidata P138 check (v43, with counts).
- (d) The distribution of the 25 m flooded area of the excluded AOIs (by criterion).
- (e) The share of ASCII names of the analysed settlements with digraphs (ae, oe, ue).
- (f) Precision of the etymology check by language block (codes by the language model, with counts).
- (g) The correlation of held-out M1 residuals between AOIs of the same activation that lie in different places (by
  distance band).
- (h) R for the European AOIs alone (European yardstick, cluster bootstrap over European places) and the flooded share
  by region.
- (i) A bootstrap test of the H3 slope (linear slope of R of the ladder T0–T4 on the rung index, cluster bootstrap over
  places, 4000 replicates).
- (j) The CSI of M1 at fixed probability thresholds (0.5 and the training-side flooded share), for comparison with the
  literature.
- (k) A list of estimate categories (Supplement): every reported estimate with its category (pre-specified; post hoc
  planned before it was computed; computed before it was planned; not planned), amendment and table.

## Reporting rules

- All quantities above are reported whatever their results (in the text or the Supplement). The confirmatory verdicts
  are not changed.
- All post hoc intervals are descriptive, and none is interpreted on its own without corroboration (stated in the text
  as the policy on multiplicity).
- Follow the rule of item 1 and the rule of item 2 (the abstract gives the translation as a range).
- The order of computation, thresholds and group definitions are not changed after the results are seen. Any change is
  added to this document, and the results before the change are reported too.

## Addendum 1 (2026-10-27T05:20Z, before items 12–14 below were computed): implementation record and descriptions added in response to the editor's decision

Values seen when this addendum was written: the local computations of items 5–9 and 11 above (v55; the pairs bootstrap of
item 10 was still running and had not been seen) and the fold-0-only smoke-test outputs of the cluster run modes (items
1, 2, 4, 5 and 6(b)), run to check the code (the Brier scores of M1c and M3n matched the corrected run and every arm was
produced). The results of items 1–4, the verdicts of items 5 and 6(b), and items 12–14 had not been computed.

### Implementation record (items 5–11; written into the v55 code and run; this addendum records them afterwards)

- 5: The floodplain features are built from the riverine and wetland classes of the frozen-lexicon match columns (lex_*).
  The near-settlement odds ratios use the same table as the odds ratios of B-2 to B-6 (frozen-lexicon matching); the
  difference in HAND at the settlement uses the estimator of the positive control of v37.
- 6(c)(ii): The fractional outcome is the riverine flood cover (plus the untyped cover in AOIs where untyped flooding is
  counted as riverine, as in the label rule). The point estimate solves the logistic score equations (as quasi-binomial);
  the interval is cluster-robust over places.
- 7: Magnitude is computed on the analysed cells. An AOI without cells with HAND ≤ 2 m would go to the top tercile (none
  occurred).
- 8: The difference in HAND at settlements is estimated by AOI fixed-effects regressions with and without the covariates
  on the same settlements (cluster bootstrap over places, 2,000 replicates), beside the v37 estimator.
- 9: Predictions of the corrected run are the main ones; the frozen run is set beside them.
- 10: The pairs bootstrap of v50 (places drawn more than once are refitted as separate places and AOIs), 2,000
  replicates, seed 20261028.
- 11(b): An outlying replicate lies more than three robust standard deviations from the median.
- 11(d): Flooded area = sum over assessable cells of flood cover (riverine by the label rule) × 0.0625 km².
- 11(e): Shares by settlement (and weighted by cells); countries with at least 50 settlements.
- 11(f): The per-language precision of v44 is re-reported (the four languages whose coded counts reached the minimum).
- 11(g): AOI means of the held-out M1 residuals of the frozen run; correlations over pairs of AOIs of the same activation
  in different location groups, by distance between centres (0–50, 50–100, 100–250, 250–500, > 500 km); pairs of
  different activations as reference.
- 11(h): Europe = cells whose nearest settlement lies in a country of GeoNames continent EU (as in the paper's
  by-continent table). The corrected run is the main one; the frozen run and the non-European cells are set beside it.
- 11(i): Least-squares slope of R on the (centred) rung index, frozen run (as the paper's ladder), with the 4,000
  bootstrap replicates of the confirmatory aggregation.
- 11(j): M1 of the frozen run (as the CSI in the paper), thresholds 0.5 and the fold's training-side flooded share.

### Added descriptions (items 12–14; the confirmatory verdicts are not changed)

- **12. Joint tests of the coefficients of the near-settlement models** (editor's minor point 12, as recommended by Pepe
  et al. 2013): in the EMS near-settlement logistic regressions (the form of v39), Wald tests with the cluster-robust
  covariance over places of the hypotheses that (a) the three water-side × distance-band coefficients and (b) the 18
  class × distance-band coefficients are all zero (p values from χ² with q degrees of freedom and from the F(q, G − 1)
  approximation).
- **13. Matches gained by reading digraphs** (editor's major point 2, transliteration): for the ASCII names of the
  analysed settlements, also match the name with ae, oe and ue read as a, o and u (union with the original reading; the
  lexicon is not changed), and report the number and share of settlements (and of cells nearest them) that gain a lexicon
  match and that gain a water-side match, overall and by country (countries with at least 50 settlements). Also the number
  of matches lost if only the rewritten name were matched (the side effect of rewriting).
- **14. Recall of the etymology check by language** (editor's major point 2): recall of the lexicon's water-side reading by
  language, estimated as a domain with the stratum weights of the whole sample (as in v44), with percentile intervals from
  2,000 resamples within strata, beside the precision of item 11(f). The codes are the language model's codes of v44, not
  independent human codes.

## Addendum 2 (2026-10-27T06:05Z, before any result of items 1 and 2 had been seen): implementation record

Values seen when this addendum was written: those of addendum 1, plus the cluster results of items 3 (activation-grouped
split, corrected code), 4 (recurrence), 5 (verdict of the floodplain reading) and 6(b) (verdict without partly flooded
cells), and the local results of items 12–14. No fold of item 1 (rebuilt permutations) or item 2 (refitted plasmode) had
been seen. The computations are not changed.

- The permutations of item 1 (seeds 1–200, the seeds of B-5 item 1) are run in four ranges of 50 seeds (from 1, 51, 101
  and 151), as separate jobs per fold (`--perm-first`). The outputs are files per range, joined in the order of the seeds
  by the aggregation (v56).
- The run record of item 4 lists the AOIs with other events, so that the aggregation can compute R restricted to them.
- Cluster capacity: the permutation jobs used 3–4 GB of memory, and their request (20 Gi) limited how many could run at
  once. The 37 jobs that had not yet started (27 permutation and 10 plasmode jobs) were suspended (not deleted) and
  submitted again with the same arguments and ranges and lower requests (permutations 3 cores and 10 Gi, plasmode 8 cores
  and 14 Gi). The suspended jobs will not be resumed; they are listed in the research log.
- When the predictions of item 4 were fetched from the NAS, one file was truncated and two were missing; they were fetched
  again in parts with sha256 checks (every file was checked by reading the zip archives and the JSON).

## Addendum 3 (2026-10-27T06:15Z, before item 15 below was computed; the results of items 1 and 2 had still not been seen): test of the second-stage name coefficients

This answers the editor's minor point 12 (Pepe et al. 2013: test the coefficient of the new marker rather than a change
in performance) more directly than item 12 (the coefficients of the near-settlement models).

- **15. Joint test of the name-feature coefficients on the held-out predictions**: rebuild the table of the primary split
  with the corrected code (the confirmatory table and folds), and regress the EMS labels by logistic regression on the
  logit of the held-out M1c prediction of the corrected run (world_v2_lexfix) and the lexicon features of M3n (six class
  indicators of the nearest name, six distance-weighted class counts over the three nearest names and the water-side
  indicator: 13 columns), with AOI fixed effects. Wald test, with the cluster-robust covariance over location groups, of
  the hypothesis that all 13 coefficients are zero (χ² and F(q, G − 1)). The character-model score is fitted within each
  fold and is not included. The row order, y and AOI of the table are checked against the predictions of the corrected
  run before the computation.
