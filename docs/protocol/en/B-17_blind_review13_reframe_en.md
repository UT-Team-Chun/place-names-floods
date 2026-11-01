# B-17 Restructuring after the thirteenth round of first-time critiques: a study of the design's sensitivity, no claims of an upper limit, randomization inference and a check of coverage

Status: **planned (2026-10-31T07:42Z, commit cc992f3, fixed before the quantities of B and C below are computed)**.
The confirmatory verdict (equivalent-null) is unchanged. Everything here is post hoc (exploratory).
Protocol: [`preregistration_v2.md`](../preregistration_v2.md) (FROZEN 2026-10-25T20:58:37Z, commit e625924).

Values seen when this amendment was written: the confirmatory analysis and every result of amendments B-2 to B-16 (including
v77 to v79 and the reruns of addenda 3 and 4 of B-16), and the values that the four reports of round 13 read from the tables.
None of the quantities of B and C below has been computed. The summary of the distribution of R over the permutations of the
primary split (the 200 of B-7, corrected code: median, 2.5 to 97.5 %, the real names at the 28th percentile) has been seen in
earlier versions of the paper and in Table S6.

## Background

The thirteenth round of first-time critiques (NHESS version; main text 34 pp., Supplement 43 pp.; development repository at
5bf18bf, public repository at cf920a9; reports in [reviews/](../../research/reviews/), `2026-10-31_paper2_blind13_*.md`)
returned the paper before review (editor), with major revision from all three referees. The editor wrote that if the message
were still only that names did not detectably improve predictions, with neither an informative bound nor a reframing, the
paper would be declined. Common points:

1. **Claims of an upper limit**: putting the designated construction's limit (0.014 with z) in the abstract and concluding
   that the test "bounds this reading and assignment" contradicts the paper's own few-cluster limit of 0.141 (an odds ratio of
   about 4) and the six of ten constructions with no bound. The few-cluster rule is ad hoc and unvalidated (ν from realised
   contributions is noisy), and the cited sources do not support it (all four).
2. **Sensitivity of the design**: the verdict could be reached whether or not names carry moderate information. Put the
   design's sensitivity (ceiling, dilution of records, verdicts for injected effects, recall of the readings) before the
   verdict and in the abstract, and match the title to the reframing (editor M1 and M7, referee 1 M1, referee 3 M10).
3. **Results that cut against the headline**: the within-area gain of names on the GFM labels (planned odds ratio 0.55,
   shifted R 0.011, difference from reassigned names 0.028, joint Wald p < 0.001), the negative verdict of the centred
   comparison, the interval above zero on the split that keeps activations whole, and the positive verdicts after the language
   model's cut-off belong in the abstract and conclusions (editor M2 and M5, referee 1 M4, referee 3 M11).
4. **Placebo**: verdicts depend on one reassignment; use randomization inference over many reassignments (referee 2 C4 and
   C7, editor M5).
5. **Labels**: thresholds of 0.5 of the whole cell, channels in AOIs without a hydrography layer, flood traces; the label checks
   rescored fixed models. Refit the models (referee 1 M2 and M3).
6. **Baseline**: HAND from a surface model; a strong baseline with post-event information. Use a bare-earth baseline with
   pre-event information only (referee 1 M5).
7. **Readings**: unvalidated readings, the group carrying the variance read through "well" (chah), classes mixing opposite
   expectations, one coefficient for all languages (referee 3 M1 to M5 and M11, editor M3, referee 2 C8).
8. **Record**: the Supplement's timeline, list of amendments and planning status stop at amendment 15; the public README
   carries an old title (editor M6, referee 3 minor 1, referee 2 minor 4).
9. **Numbers and references**: the twelve items of the editor's 3.9 and the minor comments of referees 1 to 3.
10. **What cannot be checked**: the translations, evidence of the time of the freeze, deposit of the data, the statement on
    generative AI (authors; D).

## A. Rewriting (no computation)

1. **Restructuring**: write the paper as a study of the design's sensitivity whose main message is that this assignment
   design cannot detect the information names may carry around settlements and that the confirmatory verdict does not
   distinguish it, why, and what should be checked before settlement-level information enters gridded hazard models. Order of
   results: sample → what the design can detect (ceiling, dilution of known records, injected effects, recall of the readings)
   → the protocol's test (as a record, with the rule's operating characteristics) → within-area evidence after the result
   (randomization inference as the main test; the designated construction described without claiming a limit; the centred
   comparison and the split keeping activations whole under both readings; baseline and label checks refitted; the GFM
   discrepancy) → readings → labels near settlements. Move post hoc results out of the Methods. Cut the main text to at most
   70 % of its present length.
2. **Title** (proposal): "Machine-read settlement names did not detectably improve gridded predictions of Copernicus-mapped
   river flooding, in a test blind to settlement-level information" (20 words), naming the automated reading, the label
   source, the result and the design's limit.
3. **Abstract** (at most 200 words, counted after rendering): question → data → the protocol's record (unregistered, the rule's
   operating characteristics) → the design's sensitivity (ceiling, value kept by known records, the JRC map summarised around
   settlements also undetected) → inference after the result (no robust limit, the randomization p, the group carrying the
   variance read through "well") → the GFM discrepancy (marked post hoc) → conclusion (no answer about what names record;
   calibrate gridded tests with known records first). No "bound". Exploratory results in both directions on the same footing.
4. **Upper limits**: remove "bounds" from the text, discussion and conclusions. Keep the designated construction's values in
   Table 3 as descriptions, next to the coverage of the few-cluster rule (B1), whose result decides the wording (reporting
   rules). Give the odds ratios of the t(ν − 1) limits (3.97 to 4.18) and name the three constructions whose limit lies below
   the margin. Add D̂/SE(D̂) for each construction to Table 3.
5. **Positive control**: call the version of site HAND added to T0 (detected) the positive control; the version added to M1
   is redundant because M1 already contains each cell's HAND. Remove the M1 version from the paragraph on the design's
   sensitivity.
6. **Verdicts of injected effects** (editor 3.9 item 10, referee 2 C9 (d)): rewrite with the verdicts in the order of the
   decision table (B10).
7. **Exposure before the freeze**: state in Sect. 2.4 that the pilots showed the area-level relation between names and
   flooding that determined the estimand, and that the results of the earlier test in Japan had been seen.
8. **Record**: add amendments 16 and 17 to the Supplement's timeline (Table S3), content of amendments (S4) and planning status
   (S5); extend the census of estimates to amendments 16 and 17 (B14); replace "each number in the paper is marked" with "each
   group of estimates is assigned a planning status".
9. **Numbers and references**: the editor's 3.9 (31 % and 33 %; mixed runs for 0.043 and 0.430; the two definitions of
   language coverage; the flow of the etymology sample; the source of the margin for an effect that removes flooding; the source
   and denominator of 51 %; the two R within 500 m; the kind of interval for 1.02; the subset of −2.1 m; verdicts of injected
   effects; the missing amendment 16; "Sect. 6.2" in the list of Supplement tables) and every minor comment of referees 1 to 3,
   in particular: the sign convention of the text and Table 2; "a water-related name" → "a settlement with a water-related name"
   in the abstract; the truncated title and overlapping note in Fig. 4; that the label thresholds are fractions of the whole
   cell; the error of HydroRIVERS and the 100 m distance raster; the DEM behind the heights (FABDEM and HydroRIVERS); the GFM
   window and number of scenes; the units of Table S41; the sizes of the terciles in Table S46; the ratios of Table S9; the
   references for Satterthwaite and Benjamini–Yekutieli; citations of the Supplement's references.
10. **Limits of the readings**: state in the abstract and discussion that the readings are machine readings, that many
    languages are not read, that indigenous and minor names lie outside the test, and that the group carrying the variance is
    read through "well". Name dry-point siting (settlements named after fords, wells, springs or banks stand on dry ground beside
    water) as the hypothesis the data suggest (referee 3 M11). Move the edited lexicon and the element groups to the Supplement.
11. **Literature search**: give the method and dates in Text S1 and weaken the claim (no bibliographic databases).
12. **Statement on generative AI**: consolidate in one place; state that an agent ran the literature search, which quantities
    agents computed first, and that the Supplement was also drafted (the response to the policy is the authors'; D).
13. **References**: Stephens et al. (2012, 2014), Hawker et al. (2018), Scussolini et al. (2016), Pustejovsky and Tipton (2018),
    Satterthwaite (1946), Benjamini and Yekutieli (2001), Koch and Hercus (2009), Acheson et al. (2017) and others, each cited
    after its bibliographic details are checked.
14. **Public repository**: README and PROVENANCE with the new title and amendments 1 to 17.

## B. Local computation (script v80; random seed 20261106)

1. **Coverage** (referee 2 C4 and C6): with the saved cell predictions of the protocol's null simulation (power_v3_null, 100
   replicates) and of the simulation at the margin (power_v3_sesoi, 100 replicates), draw for each replicate 1,000 bootstrap
   resamples of the saved location groups (models fixed, seed 20261106 + replicate number) and compute (i) the protocol's 90 %
   interval (percentiles of Δ divided by D̂), (ii) the 90 % z Fieller interval (its upper end as the one-sided 95 % limit), (iii)
   the rule of B3 of amendment 16 (t(ν − 1), ν iterated at the upper limit of L+; no limit if ν − 1 < 1 or the set is
   unbounded). The true value is the mean point estimate of R over the 100 replicates of each scenario (likewise for Δ).
   Coverage = share of replicates whose interval contains the true value (two-sided 90 %) and share whose upper limit is at least
   the true value (one-sided 95 %; "no limit" counts as covering and its share is reported). Wilson 95 % intervals. The same
   for Δ alone (Brier units) with percentile and t(G* − 1) limits. The distribution of ν as well.
2. **Operating characteristics with the margin fixed** (referee 2 C1): with the bootstrap of B1 of B-12 (`b1_operating` in
   v69; 2,000 outer × 1,000 inner, seed 20261106), the share of outer resamples in which the rule with its margin fixed at the
   realised 0.05 D̂ (the inner 90 % percentile interval of Δ divided by the full-sample D̂) declares equivalence, for r0 ∈
   {−0.05, 0, +0.05} and both shift patterns, with Monte Carlo standard errors.
3. **Randomization inference** (referee 2 C4 and C7, editor M5): the position of the real names in the distribution over
   reassignments (composition rebuilt, second stage fitted). Sources: the 200 of B-7 on the primary split (corrected code) and
   the runs of C2 (primary split with the frozen reading; split keeping activations whole under both readings) and C3 (GFM
   labels). For the nested (M1c to M3n) and centred (M1s to M3d) Δ_j: one-sided p = (1 + #{Δ_j ≥ Δ̂})/(1 + N), two-sided
   p = (1 + #{|Δ_j − Δ̄| ≥ |Δ̂ − Δ̄|})/(1 + N), the placebo-adjusted R − mean(R_j) and the 2.5 to 97.5 % range of R_j. On the
   language-model cut-off split (the 200 redraws of C1): for each of the six comparisons, R − mean(R_p), p = (1 + #{R_p ≥
   R̂})/(1 + N) and the frequencies of verdicts (both readings).
4. **Second stage refitted on other label definitions** (referee 1 M2): join to the B1 dump of B-16 (primary split, frozen
   reading) the fractions of the training and evaluation cells (f_riverine, f_flood_any, f_trace, f_perm_water and f_observed of
   cells.parquet, by aoi, ix, iy) and GSW occurrence (cell_gsw); refit the second stage (M1c, M3n, M1s, M3d; base M0 and M1
   fixed) and score with these labels: thresholds 0.10 and 0.25 (y1 = riverine cover ≥ t, y0 = cover of all types < t, others
   dropped); flood traces removed (notrace); share of non-permanent-water area (nonwater); share of observed land (landobs:
   denominator f_observed − f_perm_water, dropped if not positive); cells with GSW occurrence of at least 50 % removed from
   training and evaluation (gswmask); the fraction as target (riverine cover over observed land, second stage fitted by a
   binomial log-likelihood weighted by the fraction, Brier against the fraction). For each definition: Δ, D (M0 and M1 scored on
   that definition), R and 90 % Fieller interval (4,000 draws of location groups, models fixed), the median change of within-AOI
   AUC (binary definitions) and the capped EMS estimate of Eq. (2) (binary definitions).
5. **The same comparison on GFM and EMS** (referee 1 M4): with the predictions of the confirmatory run (trained on EMS, models
   fixed), score (a) the EMS labels of the cells GFM can assess and (b) the EMS labels of every AOI, after shifting each AOI's
   predictions to its flooded share (as in the GFM comparison of B-15: a shift on the logit scale so that each AOI's mean
   prediction equals its flooded share, with the baseline gain from shifted predictions too). R, real minus reassigned names
   (M3np), 95 % intervals from 4,000 draws of location groups.
6. **The GFM gain and settlement attributes** (referee 1 M4): with the dump of C7 (GFM labels, corrected code, second-stage
   design), refit the second stage with settlement attributes added to both the reference and the names models: log10(1 +
   population) of the nearest settlement in GeoNames; class of its feature code (PPL; PPLA to PPLA4 and PPLC; PPLX; other); the
   mean WorldCover built-up and tree shares of the cells within 500 m of the nearest settlement whose nearest settlement it is;
   the first block that the lexicon assigns to the country of the nearest settlement (blocks with at least 2 % of the cells, and
   other). Δ (Brier), shifted R, difference from reassigned names (M3np with the same attributes), 4,000 draws of location groups.
7. **N-grams of the character model** (referee 3 minor 17): fit the character model (`nlearn_scores` settings) on all training
   cells of the primary split, with EMS and with GFM labels, and list the 20 one- to four-character n-grams of the names with
   the largest and smallest weights (coefficients of their hash buckets), with the number of names in which each occurs.
8. **Coefficients by language** (referee 2 C8, referee 3 M4): with the B1 dump of B-16 (frozen reading), add to the names model
   the 14 features and their interactions with the language area (the first block of the nearest settlement's country; blocks
   with at least 2 % of training cells, and other), fitted with an L2 penalty C ∈ {0.01, 0.1, 1} (standardised features). R and
   90 % Fieller interval (4,000 draws of location groups, models fixed).
9. **Decomposition of the dilution** (referee 2 C13): with the predictions of the existing runs for the records (per-cell and
   through-settlement versions of GSW from 2022, earlier EMS floods and the JRC 100-year map; models fixed), compute Δ on all
   ground and on the ground within 500 m of the nearest settlement, and the share kept by the through-settlement version on the
   500 m ground (change of variable) and on all ground (adding the change of support), with Fieller intervals (4,000 draws of
   location groups).
10. **Verdicts of injected effects** (editor 3.9 item 10): recount the verdict of each replicate of C4 of B-9 (every model
    refitted) in the order of the decision table (positive verdicts before equivalent-null).
11. **Correlation of the features** (editor M3): correlation matrix and variance inflation factors of the 14 name features on
    the training cells of the B1 dump.
12. **Flow of the etymology sample and definitions of language coverage** (editor 3.9 items 3 and 4, referee 3 minors 7, 19
    and 20): one flow from sampled names → names with a section on the name → names with a derivation → coded names →
    published names, from the records of v44, v57b and v69; the two quantities of language coverage (cells in countries without
    a block; cells in countries whose first language no assigned block reads) and the definition and source of "regions of
    another language".
13. **Subset of read names** (referee 3 M1 (e)): R scored only on cells whose nearest name matches any element of the lexicon
    (predictions of the point refit of B1 of B-16, models fixed), and R by tercile of the share of read names in each AOI.
14. **Census**: extend the census of estimates of B8 of B-15 to the interval estimates of amendments 16 and 17 (same families,
    BY q).
15. **Multiplicity**: the number of pre-specified verdicts, and the number of verdicts other than equivalent-null expected by
    chance from the protocol's null simulation.

## C. Cluster computation (utens; image tags recorded; the frozen reading via `--frozen-norway`)

1. **Placebo redrawn to 200**: language-model cut-off split, frozen reading, `--nread`, seeds 20261105 + k (k = 21 to 200),
   arms M0, M1, M3, M3p, M3r, M3rp, stacks M1c, M3n, M3np, M1s, M3s, M3ps, M1rc, M3rn, M3rnp, M3rs, M3rps
   (world_v2_b17plac021 to 200). The three comparisons that use the lexicon join the 20 frozen-reading redraws of B-16 (200 in
   all); the three that use only N-read join the 20 corrected-code redraws of B-16 (200 in all; the N-read comparisons are the
   same under both readings).
2. **More permutations** (`--perms-recomp`, 200, seeds 1 to 200, composition rebuilt): (a) primary split, frozen reading
   (world_v2_b17perm_frozen); (b) split keeping activations whole, frozen reading (world_v2_b17permact_frozen); (c) split keeping
   activations whole, corrected code (world_v2_b17permact).
3. **Permutations on the GFM labels**: `--labels gfm --perms-recomp`, 200, primary split, corrected code (world_v2_b17permgfm).
4. **A bare-earth baseline with pre-event information only** (referee 1 M5): arm M1b = M1 + the six FABDEM elevation
   statistics + heights above channels draining at least 100 and 1,000 km² measured on FABDEM (C7 of B-9) + log10 of the
   upstream area of the nearest channel (C1 of B-16); no JRC maps, GSW or WorldCover. Stacks M1bc, M3bn, M3bnp, M1bs, M3bd,
   M3bdp. Primary split and the split keeping activations whole, frozen reading (world_v2_b17bare, world_v2_b17bareact). R with
   the baseline gain of M1b, BS(M0) − BS(M1b), and with the main one. Verdicts from `_verdicts` of v56.
5. **Both stages refitted on other label definitions**: primary split, frozen reading, arms M0 and M1, stacks M1c, M3n, M3np,
   M1s, M3d, M3dp. Definitions: threshold 0.10 (world_v2_b17t010), 0.25 (world_v2_b17t025), share of observed land
   (world_v2_b17landobs), GSW occurrence of at least 50 % removed (`--exclude-gsw-perm 50`, world_v2_b17gswmask). Add t010,
   t025 and landobs to `riverine_labels` (the definitions of B4; tests added).
6. **Lexicon without wells and -ab**: a derived lexicon that drops chah, chahi, ab and aab (as words and as suffixes) from
   block fa of the frozen lexicon (docs/verification/v80_lexicon_no_chah_ab.yaml, made mechanically); primary split, frozen
   reading, arms M0 and M1, stacks M1c, M3n, M3np, M1s, M3d, M3dp (world_v2_b17nochah).
7. **Second-stage dump on the GFM labels**: `--labels gfm --dump-stage2`, primary split, corrected code (world_v2_b17gfmdump);
   used in B6.
8. **Scripts**: local computation in `backend/verification/r_okauchi/v80_b17_local.py`, aggregation of the cluster runs in
   `v81_b17_cluster.py`; output `docs/verification/v80_b17_local.json` and `v81_b17_cluster.json`.

## D. Decisions for the authors (not made in this amendment)

1. **Policy on generative AI** (a precondition of review): NHESS does not allow generative AI for text and interpretations.
2. **Toponymic check** (referee 3's condition 1): carry out the v69 check with human specialists, including the blocks and
   countries that carry the estimates (PK, BD, KH, MM, AU, DE, PT, GB), and report error rates per block. Without it, cut the
   toponymic interpretation (Sect. 3.5, most of Sect. 3.6, Table S26).
3. Checking the translations of the protocol and every amendment.
4. Verifiable evidence of the freeze (the development repository is not to be shared: a third-party time stamp, or an archive of
   the freeze commit and its ancestors only).
5. Deposit of the derived data (Zenodo DOIs) and access for reviewers; a DOI for the code; the terms for redistributing
   FABDEM-derived products.
6. Deposit of the conversation and prompt that produced the lexicon.
7. A statement of what each author checked, and CRediT roles consistent with it.
8. Whether to publish machine readings of indigenous names (CARE).
9. Format and journal.

## Reporting rules

- Everything is post hoc (exploratory), and the confirmatory verdict does not change. The protocol's verdict is reported as the
  record, with its operating characteristics.
- The abstract and conclusions do not use "bound". The wording of limits follows the coverage of B1: if the one-sided 95 % limit
  of the few-cluster rule exceeds the true value in at least 90 % of replicates in both the null and the margin simulations, its
  limit may be reported as "a limit from a rule checked by simulation". Otherwise no limit is claimed, the paper states that no
  limit with checked properties is available, and the values stay in a table as descriptions.
- The main post hoc test within areas is randomization inference (B3). Placebo verdicts are reported as frequencies over many
  draws.
- Every quantity of B and C is reported, whatever its direction. If exploratory results enter the abstract, results in both
  directions (for example the GFM discrepancy and the randomization p) enter on the same footing.
- Rules, units, numbers of replicates and strata are not changed after computation. Implementation errors are corrected in
  addenda.

## Addendum 1 (record of the implementation of B and C; 2026-10-31T08:11Z, commit 3bf202a; before the quantities of C1 to C7 and of B3, B6 and B14 were computed)

Quantities already computed when this addendum was written: B1, B2, B4, B5, B7 to B13 and B15 (v80), and, of B3, the position
of the real names among the permutations of B-7 (primary split, corrected code; a trial run of v81: nested R −0.0040,
one-sided p 0.72). Items 1 and 2 are discrepancies in the text found after computing; the computations do not change.

1. **Number of name features (B8, B11)**: the amendment says "14 features", but the names model has 15 (Table 1: six class
   indicators, six distance-weighted class counts over the three nearest names, the water-related indicator, the indicator that
   the lexicon reads a language of the country, and the character model's score). The coverage indicator is 1 for every
   analysed cell with the corrected code, so effectively 14 (\NNameFeatEff in the text). The B1 dump uses the frozen reading,
   in which the indicator is 0 for the Norwegian cells, so all 15 vary. B8 and B11 used the 15 features of the names model (the
   code was the same before computing). The tables say 15.
2. **Where B3 is computed**: B3 gathers cluster runs (C1 to C3 and B-7), so it is computed in v81, not v80, and written to the
   C1 (language-model cut-off split) and C2 (permutations) entries of docs/verification/v81_b17_cluster.json.
3. **Yardstick of the permutations**: the yardstick BS(M0) − BS(M1) that turns Δ_j into R_j is taken from a run on the same
   split as the permutation run (M0 and M1 do not use names, so every run on the same split gives the same values; the
   predictions of world_v2, world_v2_lexfix and world_v2_b13lexfix, and of world_v2_actgrp, world_v2_b7actfix and
   world_v2_b16strongact, were checked to differ by 0). Primary split: world_v2; split keeping activations whole:
   world_v2_actgrp; GFM: world_v2_gfm. As a check of the implementation, the location-group sums of the real M1c, M3n, M1s and
   M3d saved by each permutation run are compared with a run of the same reading and code (B-7: world_v2_lexfix; C2 (a):
   world_v2; C2 (b): world_v2_actgrp; C2 (c): world_v2_b7actfix; C3: world_v2_gfm) and the largest relative difference is
   reported (3 × 10⁻⁸ for B-7). Only permutation runs with all 10 folds × 4 chunks (50 each) are used.
4. **Real R in C1**: only the placebo arm changes between draws, so the real names' R is that of the confirmatory run
   (world_v2, llmcut), and the largest difference from each draw's R is reported as a check. For the three comparisons that use
   the lexicon, the verdicts of the 20 corrected-code draws (B-16 plac01 to 20) are also listed as a description.
5. **Submission of C7**: the first 10 jobs, submitted with every arm, were heavy and were suspended; they were resubmitted with
   `--arms M0,M1 --stacks M1c,M3n,M3np,M1s,M3d` (only the arms B6 needs; same tag world_v2_b17gfmdump).
6. **Arms of C4**: fitting M1b needs `--strong-baseline` (the flag of C1 of B-16). Submitted with `--arms M0,M1,M1b --stacks
   M1bc,M3bn,M3bnp,M1bs,M3bd,M3bdp`.

## Addendum 2 (descriptive checks for the remaining requests of the thirteenth round; 2026-10-31T08:19Z, commit 168452e; before the quantities below were computed)

Requests of the editor and referees that the body of this amendment did not take up and that can be counted locally. All are
exploratory, and the confirmatory verdict does not change. Script `backend/verification/r_okauchi/v82_b17_addendum2.py`, output
docs/verification/v82_b17_addendum2.json. Random seed 20261106.

1. **Leakage through the regional composition, shown directly** (editor M2): for each held-out AOI of the B1 dump of B-16
   (primary split, frozen reading), the AOI mean of the character model's score (the N-learn column of M1c's regional
   composition) and the AOI's flooded share (mean over cells). Among the other analysed AOIs of the same activation (the EMSR
   number in the AOI's name), those on the training side of its fold are its "training-side siblings"; their flooded share is
   the cell-weighted mean. Spearman correlations across AOIs (i) between the mean score and the AOI's own flooded share, for
   AOIs with and without training-side siblings, and (ii) between the mean score and the siblings' flooded share, for AOIs with
   siblings. 95 % intervals from 1,000 resamples of activations (seed 20261106).
2. **Number of GFM scenes** (referee 1, minor 3): for the 192 AOIs with GFM labels, the distribution of the number of
   Sentinel-1 scenes in the window (n_scenes of gfm_summary.json), and Spearman correlations of the number of scenes with the
   GFM flooded share (assessable cells) and with the share of cells mapped flooded by EMS that GFM also mapped flooded (per AOI,
   cells assessable on both, capped table).
3. **AOIs that lose their flooded cells to the cap** (referee 1, minor 11): the numbers of analysed AOIs with no flooded cell and
   with no dry cell after the cap.
4. **Matches of ab** (referee 3, minor 11): the settlements (nearest settlements of the analysed AOIs, frozen reading) matched by
   ab and aab of block fa of the frozen lexicon, counted by match as a word and as a suffix and by country; the words of the
   Pakistani names that contain a match, up to the 30 most frequent.
5. **Hedged derivations** (referee 3, minor 12): a coded derivation is hedged when its quote (v44 codes, field quote) contains
   any of these terms (case-insensitive, at word boundaries): perhaps, probably, possibly, may, might, likely, uncertain,
   unclear, disputed, believed, said to, supposedly, legend, legendary, folk etymology, tradition, traditionally, according to.
   Their number, and the lexicon's precision and recall without them (coders 1 and 2, with the weights and intervals of v44 and
   v57b).
6. **The two quantities of language coverage** (the rest of B12; already computed): "cells in countries whose first language is
   not among the languages read by the blocks applied to the country" (B3 of v69, 1.0 %) and "cells in countries whose first
   language has no block of its own" (complement of B9 (c) of v77; Urdu in Pakistan and Bengali in Bangladesh are read through
   the Persian and Hindi blocks; with the frozen reading Norway is included) in one table, with the definition and source of
   "regions of another language" (B14 of v66: a list of 13 regions of Natural Earth administrative units).
7. **Base models stopped early on location groups** (referee 1, minor 7; editor M4): C3 of B-13 (world_v2_b13gres, already
   computed) enters the text and the Supplement. No new computation.

## Addendum 3 (implementation of B6; 2026-10-31T09:27Z, commit 7704e09; before the quantities of B6 were computed)

1. **The C7 dump lacks the reassigned names**: `dump_stage2` saves the regional composition and the real name features only,
   not the placebo name features (`NAMES_P`) that B6 needs for the difference from reassigned names (M3np with the same
   attributes). `dump_stage2` is extended to save `xnp_tr` and `xnp_te` (the columns of M3np; a test is added), and C7 is rerun
   locally (world_v2_b17gfmdump2_local: primary split, GFM labels, corrected code, `--dump-stage2`, 10 folds). The cluster run
   of C7 (world_v2_b17gfmdump) is used only to check that the real arrays agree with the local dump.
2. **How the attributes are built** (making the wording of B6 concrete): population is GeoNames' population (zero kept as
   zero, log10(1 + population)); the feature class has four levels from the feature code (PPL; PPLA to PPLA4 and PPLC; PPLX;
   other; indicators with PPL as reference); built-up and tree cover are each cell's WorldCover lc_built and lc_tree
   (features.parquet) averaged, for each settlement, over the cells within 500 m whose nearest settlement it is (zero with an
   indicator for settlements without such cells); the language area is the country's first block (blocks with at least 2 % of
   the fold's training cells, the rest pooled). They are added in the same way to the second stages of the reference (M1c),
   names (M3n) and placebo (M3np) models, refitted with C = 1 on standardised features. As in the GFM comparison, each AOI's
   predictions are shifted to its flooded share, and R and the difference from reassigned names are given with 95 % intervals
   from 4,000 resamples of location groups (models fixed, seed 20261106).

## Addendum 4 (where the computation runs; 2026-10-31T10:00Z, commit d1dd75c; before the quantities of B6 were computed)

1. **C7 rerun and B6 on utens**: at the user's instruction to avoid computation on the user's Mac, the dump with the
   reassigned names, which item 1 of addendum 3 placed on the local machine, is rerun on utens (world_v2_b17gfmdump2: primary
   split, GFM labels, corrected code, `--dump-stage2`, with the arms and stacks of the cluster C7, `--arms M0,M1 --stacks
   M1c,M3n,M3np,M1s,M3d`, 10 folds; code at HEAD after 3a15fda, with the second stage and the dump as in 3a15fda). The refits
   of B6 are also computed in a utens job. The partial local dump (world_v2_b17gfmdump2_local, 8 of 10 folds) is not used.
2. **Checks**: B6c reports whether the real arrays (common keys) of the cluster C7 (world_v2_b17gfmdump, Miyabi-C) and of
   world_v2_b17gfmdump2 agree, and B6 reports whether the refits without attributes reproduce the saved predictions of M1c and
   M3n (largest difference). The definition of B6 (attributes, refits, scoring, resampling) is that of item 2 of addendum 3.
