# B-12 Restructuring and exploratory analyses in response to the eighth round of first-time critiques (NHESS version, 2026-10-29)

English translation, drafted with a language model; the Japanese file `docs/protocol/amendments/B-12_blind_review8_analyses.md`
is the record and governs where the two differ.

Status: **plan (2026-10-29T09:42Z, commit 4b1b327, before any quantity in B and C below was computed)**. The confirmatory
verdict is not changed. Everything below is **post hoc** and responds to an editor and three referees who read only the
PDFs of the manuscript and its supplement for the first time (language-model agents of a kind that does not load the
project's CLAUDE.md, told that they could consult the public repository at b6fb61e; reports in
`docs/research/reviews/2026-10-29_paper2_blind8_*.md`, prompts in the `_prompts.md` file in the same folder). Protocol:
`preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924). Code: with the correction of amendment B-4 (the
corrected code). Items C1 to C4 of B-9 are still being computed on the cluster, and this amendment does not change how they
are computed.

Values seen when this document was written: the confirmatory analysis and all results of B-2 to B-11 (except C1 to C4 of B-9,
and the items of B-10 C4 and B-11 B3 that use them). In particular, item **B7** of B-11 (an interval that refits the second
stage and recomputes the denominator in each replicate: with two-way weights over activations and countries, 90 % −0.055 to
0.036, **not** inside the margin; with location-group weights −0.020 to 0.008, inside the margin) and item **C1** (on the
country split, the real names lay at the 69th percentile of 100 reassignments with the composition rebuilt, one-sided
p = 0.31). Values estimated by the reports (referee 2's estimate that the fixed-margin rule would declare equivalence at
R = 0.05 with probability about 0.3 under an independent normal approximation; that the one estimate surviving BH is the
GFM odds ratio of the highland class, 2.45; the editor's p(1 − p) estimate of regional weights). None of the quantities in
B and C below has been computed.

## Main points (common to the four reports, or weighty)

1. **Policy on generative AI (editor)**: the NHESS policy (checked on its public page) says "Generative AI must not be used
   for text or interpretations" and allows only assistive use (sentence structure, readability, grammar). Sect. 2.7 and the
   acknowledgements state that a language model drafted the manuscript. The English translations have not been checked by
   the authors. → for the authors (D1). This amendment corrects the content of the manuscript but does not make the text
   the authors' own
2. **Title and abstract**: name the reading (matches of a lexicon compiled by a language model) and the assignment (nearest
   cells); state what the margin means (odds ratio), the ceiling (within 500 m) and the full uncertainty (intervals that
   include refitting and the denominator). The odds ratios near settlements were chosen after all versions had been seen
   and omitted the planned GFM estimate (0.55, 0.37 to 0.82), so remove them from the abstract or show all
3. **Uncertainty**: the interval with fixed models is narrow; B7 (two-way) is not inside the margin; achieved bounds range from
   0.011 to 0.046 across measures of the denominator; the operating characteristics of the fixed-margin rule are not shown;
   the Nadeau–Bengio interval was misapplied (a formula for repeated simple splits applied to the spread of 13 complete
   cross-validation estimates)
4. **Leak in the primary split and the reference model**: AOIs of one activation straddle folds; M1c beats M1s only on this
   split. Put the activation split beside the primary result, with both runs and the inflation of denominators where M0 is
   worse than a constant forecast
5. **Odds ratios near settlements**: not identified (label error varies with distance from settlements); post hoc choice;
   the estimate surviving BH (highland, GFM) is not reported; the pattern by class (floodplain terms in the predicted
   direction, water and crossing the other way) is not in the main text; the composition of the uncapped sample is not given
6. **Readings**: systematic problems of the lexicon (relative qualifiers as lowland/highland, floodplain micro-relief as
   highland, homographs, missing elements); country-level language assignment and the error in "13 %" (countries counted
   because codes differ although a block reads the language); the etymology check scores names, and p0 is not in the main
   text; N-read may use knowledge of places
7. **Labels**: the typing rule (any description without words such as river is "other", counted as riverine inland), the
   GFM time window (from the day before the earliest to the day after the latest EMS image), detail maps, built-up land,
   image timing → state precisely as limitations
8. **Length and structure**: main text about two thirds as long; Sect. 3 only pre-specified analyses; the Supplement only
   cited tables, with uncited, superseded and withdrawn tables in a separate PDF in the public repository
9. **Numbers and references** (listed in A12)

## A. Rewriting (no computation; corrections in the paper and Supplement)

1. **Title**: name the reading, the assignment and the within-area estimand, keeping "pre-specified (not registered)".
   Draft: "Water-related terms matched by a lexicon in settlement names, assigned to the nearest 250 m cells, did not improve
   within-area predictions of mapped river flood extent: a pre-specified test in 500 Copernicus EMS areas" (final wording
   decided when writing; it must contain the three elements and must not generalise to "landscape terms in settlement
   names")
2. **Abstract** (100 to 200 words): (i) what was tested (water, wetland, riverine, crossing and lowland terms matched by a
   lexicon compiled by a language model, given to all nearest 250 m cells; the matches recover about one fifth of the
   water-related derivations recorded on Wikipedia, and no toponymist has checked them), (ii) the result (R, the
   pre-specified 90 % interval, and the interval that refits the second stage and resamples activations and countries,
   −0.055 to 0.036, read as excluding an improvement above 0.036 of the baseline gain but not a loss of the size of the
   margin; this one-sided reading is post hoc), (iii) what the margin means (an injected odds ratio of 2.33 on all nearest
   ground; within 500 m not reached at 20; this ceiling was found after the result and the interpretation plan was
   abandoned), (iv) the dilution control (surface-water record, R 0.258 → 0.006, exploratory), (v) no odds ratios near
   settlements, (vi) the open question. Define EMS
3. **Sect. 3 only pre-specified analyses**: the H2 verdict; H1 (predicted positive; failed in the nested form); H3 (predicted
   monotone; failed); H4 (verdict on the 90 % interval; total increment positive-small, with the limits of that reading);
   splits (Table C1); sensitivity analyses (Table C2; "eight runs of seven analyses" not counting the defective first run of
   the sixth); **the protocol's two-way robustness interval for each split** (including the country split's 95 % interval
   −0.035 to 0.070, which exceeds the margin); H5 descriptive. Achieved bounds, refitting, CSI, other measures of the
   denominator, the activation split, other estimands, relaxed second stages, the JRC maps and the between-area regression
   move to a new exploratory section (A4)
4. **New exploratory section, "Uncertainty that the protocol's interval leaves out"**: second stage refitted (location
   groups −0.020 to 0.008; two-way −0.055 to 0.036); the range of achieved bounds (0.011 to 0.046 across denominators; 0.017
   refitted); the operating characteristics of the rules (B1); the activation split (both runs with the same weight, a row in
   Fig. 4; not-meaning is a 95 % rule, so give the 95 % interval; list every split on which M0 is worse than a constant
   forecast and show the inflated denominators (B17)); the leak (−0.74 from 12 pairs, with an interval (B10)); estimands
   without the composition (scoring within AOIs only; names centred within AOIs); relaxed second stages (distance-band
   interactions, the penalty on name coefficients only (B11), country indicators, continent interactions). Do not use the
   Nadeau–Bengio interval (report the range and SD of R over the 13 fold assignments as a description). For the 13 outlying
   null replicates, state that their own bootstrap standard errors were 3 to 10 times the median and unrelated to the
   denominator (v55), with the diagnosis of B5
5. **Sect. 2 (methods)**: the assignment was chosen for convenience, not on toponymic grounds (the protocol gives no reason;
   nor for three nearest settlements and the weight 1/(1 + d/1 km)). The second stage uses training-side base logits
   cross-fitted over five inner folds of location groups, and test cells use the base model fitted on all training cells.
   The regional composition in the estimand was chosen from two pilot-3 replicates without an effect (say two). The margin
   is a convention (5 % of the baseline gain, about the improvement from recalibrating the baseline model). The leak of the
   primary split (AOIs of one activation straddling folds) in Sect. 2.5. The typing rule exactly (descriptions without river,
   fluvial, flash, torrent, coast, storm surge or tidal are "other", and "other" counts as riverine in inland AOIs). The GFM
   time window exactly (from 00:00 of the day before the earliest EMS image to 24:00 of the day after the latest, with the
   distribution of B12). The depression index (cell mean elevation minus the minimum of the 3×3 neighbourhood of cell
   means); HAND resampled bilinearly to the 25 m sub-grid; the channel threshold of ASF HAND (checked in its documentation).
   Table 1: 15 name features, one constant, effectively 14. In the note to Table 1, state that positive verdicts also require
   the placebo condition
6. **Sect. 4 (assignment through settlements) as the centre of the exploratory work**: the 0.081 of the recurrence record is
   on the denominator of 227 AOIs; give the pooled-scale value (B7'). Rephrase as a statement about scheme, reading and labels
   jointly (the recurrence record, assigned through settlements, can exceed the margin where it exists). Lead with the bounds
   with the second stage refitted (1.57 nominal, 1.83 on the estimator's scale). The history of the reference (range in
   Table S19, 1.06 → 1.48 → 1.23, the withdrawn 1.16) briefly in the main text. Correct the scale error in the 1.27 sentence
   (on the nominal scale it exceeds the fixed bound 1.23 and lies inside the refitted bound 1.57). State that the bounds are
   not test inversions (the spread is assumed equal under the alternative; Table S37 shows it widens). Add the 1 km injection
   (C1) and restrict the ceiling statement to "within 500 m" (the ceiling within 1 km, 0.099, exceeds the margin)
7. **Near settlements (shorter, within a section on limits)**: retitle "Residual deficit of flood labels near settlements".
   Add the distance-band contrasts adjusted for the 10 m built-up share (B18). State that the odds ratios are not identified
   (label error, protection, siting on dry margins and image timing predict the same sign), and show all versions (cap or
   none × EMS and GFM × domain, with planning status) in a Supplement figure (specification curve, B19). In the main text,
   give the planned GFM estimate (0.55, 0.37 to 0.82, 43 cells in 16 groups) beside the uncapped one (0.82). Name the
   estimate surviving BH (highland class, GFM, AOIs whose EMS map was not drawn from Sentinel-1: 2.45, 1.51 to 3.98; EMS 1.96,
   1.24 to 3.08) and give the selection by GFM's exclusion of ground high above drainage as a candidate explanation (B6).
   The pattern by class (riverine 1.41, wetland 1.60, floodplain reading 1.45, island and water-meadow without au 2.76 in
   the predicted direction; water 0.77 and crossing 0.77 the other way) as post hoc hypotheses with the numbers of location
   groups and countries. The country and activation composition of the uncapped sample and the estimate without Pakistan
   and Bangladesh (B2)
8. **Readings (shorter)**: a table of the reading in the main text (share of cells with a match and with a water-side match by
   continent, top forms, precision, recall and p0 (0.18, 0.14 to 0.23), share read in another language (the corrected value
   of B3)). The systematic problems of the lexicon briefly (relative qualifiers, floodplain micro-relief, homographs,
   missing elements, one word in different classes in different blocks). N-read may contain knowledge of places (of 139
   settlements named after a water feature in Wikidata, N-read 66 and the lexicon 9). "Not of names in the toponymist's
   sense" → "not of etymologies, which toponymists establish from early forms". -näs is a water-side element. River names
   are 139 of 299 recorded namesakes and cannot be read by the lexicon. Forms stating a riverside position cover 0.3 % of
   cells
9. **Introduction and Discussion**: say in the introduction too that local knowledge is not tested. The place of the scheme
   (we found no hazard model that uses place names as predictors; the scheme is the simplest way to turn the recommendation
   to read names as warnings into a gridded feature, not a test of existing practice). A joint statement about the three.
   What a settlement-scale test needs (groups and directions fixed by toponymists by referent and setting, assignment to
   the referenced ground (same bank, same catchment), built-up footprints, many events (the GFM record), pre-regulation
   floodplains and wetlands, living and fossil terms, protection, principles for indigenous names (CARE))
10. **Sect. 2.7**: list, one line per amendment, the analytical and interpretive decisions that came from language-model
    critiques, and state that the authors take responsibility (to be confirmed by the authors; D1)
11. **Supplement**: only tables cited by the paper or the Supplement text. Uncited, superseded and withdrawn tables move to a
    separate PDF in the public repository (`supplement_archive.pdf`), as the Supplement's contents state. Fix the blank
    page 13
12. **Corrections of numbers and references** (all corrections of existing values):
    - Achieved bound of the confirmatory run in the abstract (0.011; done). The sign of ΔBS ("reduction in Brier score",
      negative = worse) consistent in the text, Table 2, Table C1 and the abstract (partly done). The interval of the run
      from 2016 onward labelled 95 % (done)
    - R on the GFM domain: give the ratio-bootstrap interval (−0.052 to 0.013) and the fixed-margin interval (−0.025 to 0.009,
      Table S64) with their kinds
    - "Unknown" of N-read 6.0 % (check the rounding rule). Denominators of 14 % (nearest settlements of analysed cells) and
      11 % (settlement points) for coordinates on whole minutes. Explain the difference between Tables S23 and S24 (883/882 and
      intervals)
    - Cross-references: 0.55 (Tables S25 and S64, not S11/Fig. S3); R on the GFM domain (Table S64); content of amendments
      Table S4, times Table S3; bounds under other references Tables S18 to S19 in the Table 4 caption; element groups also
      Table S23
    - Table S21 caption: sources of p1 and p0 (Table 5 and Sect. 6.2; p0 from v44). Table S25 column 6 counts the AOIs in each
      model. Table S26 counts AOIs once per country (with overlap). Table S63 Groups = groups holding flooded cells (52 for the
      capped row). Table S77 Sect. 5 rows for the uncapped estimate. Table S13 g00098 etc. are units of the activation split.
      Check the point estimate and upper limit of the fa row of Table S8
    - The garbled "English -mouth, -ey and 87 % wæsse" in Text S2. Equivalence for H4 on the 90 % interval (−0.016 to 0.006)
    - "Largest injected effect" (the estimator-scale calibration goes to nominal 3). The 1,541 relabelled cells refer to the
      capped table. "48 activations" are those with flooded cells near water-side names. Fig. S4 description (label-free (a),
      flooded shares (b)). One sentence on the flexible comparison of the language-model split (R 0.168, not-meaning).
      Explain "capacity to memorise places" in Fig. 1. Remove the time-split point in Fig. 4 (R undefined)
    - Table 2: remove the duplicated T4/flexible row; add fixed-margin intervals for flexible "lexicon only" and "learned score
      only" (B8). Table 5: intervals for the second coder (B9)
    - The approximate p of the census (from asymmetric intervals) conflicts with the 0.05 rule; remove the histogram. Report TOST
      p as "< 1/B (0.00025)"
    - One definition of the CSI of the HAND-only map (B15)
13. **References** (checked with Crossref before entry): Bernhofen et al. 2018, Skoulikaris 2022, Zhao et al. 2021, Garousi-Nejad
    et al. 2019, Smith 1956, Gelling 1984, Watts 2004, Niemeyer 2012, Nègre 1990–1998, Tent and Blair 2011, Carroll et al. 2020
    (CARE), Ahlers 2013. Versions, URLs and terms of GSW (1984–2021), HydroRIVERS v1.0, Natural Earth and the GFM product
    definition document

## B. Computations from existing artefacts (local, script v69, seed 20261101)

1. **Operating characteristics of the decision rules**: from the held-out predictions of the confirmatory run (frozen code)
   on the primary split, sums of squared errors by location group; a nested location-group bootstrap (2,000 outer × 1,000
   inner); with the true R set to −0.05, 0 and +0.05 (each outer replicate's ΔBS shifted to centre on R times that
   replicate's denominator), the share of replicates in which (i) the fixed-margin rule (the 90 % interval of ΔBS divided by
   the replicate's point denominator) and (ii) the ratio rule (the 90 % interval of R) declare equivalence. The share at the
   boundary is the size of the rule (nominal 5 %)
2. **Composition of the uncapped sample**: flooded cells within 500 m of water-side names (EMS, every evaluable cell) by country
   and activation, with shares and AOI counts; the odds ratio without Pakistan and Bangladesh (the formula of B-11 B2, pairs
   2,000, frequency-weighted form)
3. **Recount of language assignment**: whether the first GeoNames language of the country of each cell's nearest settlement is
   in the language list of the blocks applied to that country (the `languages` of the lexicon headers); share of cells not
   covered at country level; region-level (admin-1) flags counted separately
4. **Sequential inclusion counts**: exclusions at each step when the five criteria are applied in order (from v54)
5. **Outlying null replicates**: for the 13 replicates of power_v3_null, the location group with the largest contribution to
   ΔBS, with its country and flooded share
6. **Highland class**: record the estimate surviving BH (v41) and compare the share of cells within 500 m of highland names that
   GFM can assess with that near other names
7. **Recurrence record on the pooled scale**: the increment (ΔBS) of the recurrence record assigned through settlements on the
   227 AOIs, divided by the baseline gain of all cells
8. **Flexible lexicon-only and learned-score-only rows of Table 2**: fixed-margin 95 % intervals (location groups, 4,000)
9. **Second coder's intervals**: name-level bootstrap 95 % intervals for coder 2's precision and recall in Table 5
10. **Leak and reference model**: (a) on the primary split (corrected code), refit M1c without the AOI mean of the character
    model's score from the second-stage design (world_v2_b8stage2_local) and compare Brier scores with M1s and M1c;
    (b) a Fisher-z 95 % interval for −0.74 (12 pairs)
11. **Penalty on name coefficients only**: in the same design, scale up the base-logit column so that the penalty acts only on
    the name features (and composition), widen the penalty grid until the optimum is interior, and choose by five inner folds
    of location groups; R and the chosen penalty
12. **GFM time window**: distribution of the difference between each GFM image time and the EMS image day (from the GFM
    summary records)
13. **Flood typing**: number of AOIs counted riverine by the inland rule, number of AOIs containing typed non-riverine flooding,
    and the share of flood-trace polygons by sensor (from the product records)
14. **Cells neither flooded nor dry**: number of cells left out because most of their flooding was not riverine
15. **CSI of the HAND-only map**: one tie-breaking rule (random order) and all tables with that value
16. **Skill of the JRC maps themselves**: CSI of a map ranked by the JRC flooded fraction at RP100 (flagging as many cells as
    flooded), if the features are available locally
17. **Comparison with a constant forecast**: Brier scores of M0 and of a constant forecast (the training-side flooded share) on
    each split
18. **Built-up share at 10 m**: the distance-band contrasts (formula of Sect. 5.2) and the water-side odds ratio (capped table,
    pairs 2,000) with the ESA WorldCover built-up share (lc_built, already in the terrain feature tables) added
19. **Specification curve**: a figure of every version of the water-side odds ratio near settlements with planning status, label
    set, sample and domain (existing values only)
20. **p0**: p0 from v44 with its interval as a quantity of the main text
21. **Revised design of the audit sample** (before the toponymic check; uses no results): add a stratum of names not read as
    water-side (100 names from the AOIs that carry the near-settlement estimates, stratified by country), making a sample of
    250 names in which coders do not know which were read as water-side. In the guide, add the share coded twice (at least
    20 % per language) and how disagreements are resolved, a rule for site versus feature, how relative qualifiers are
    treated, and not to consult reports of the floods. The answers are not published

## C. Cluster computation (utens, corrected code)

1. **Injection within 1 km**: shape B with the distance set to 1,000 m (shape BK: cells whose nearest name is water-side and
   within 1,000 m of that settlement), the seven δ of B-6 (odds ratios 1.28 to 20.1) with the base models fixed, 20 replicates
   each (replicates 0 to 19; salts determined by shape, δ and replicate). Read the odds ratios at the margin and at the
   observed upper 90 % limits (fixed, fold spread, refitted) with the reference of amendment 8 (mean of the rebuilt-composition
   placebo)

## D. Not done in this amendment (for the authors, or a separate study)

1. **Making the text the authors' own** (NHESS policy), or choosing a journal whose policy fits: the authors
2. Checking the translations, archiving the development history with a third party, depositing the derived data (Zenodo,
   accessible to reviewers; the manifest is prepared in `data/zenodo/`) with a tagged release, the toponymic check (sample
   and guide of B21): the authors
3. Rebuilding the labels (detail maps, a 10 m built-up mask, per-image GFM, image time relative to the peak, protection), the
   baseline (HAND along flow paths, GFPLAIN250m, the JRC European maps, monthly pre-event GSW) and the lexicon (groups by
   referent and setting, language per name), and controls of N-read for knowledge of places: inputs to the design of the
   settlement-scale test (a separate study registered with a third party). Written as limitations in the paper
4. The confirmatory verdict is not changed

## Reporting rules

- B and C are all post hoc, decided after the confirmatory analysis and the results of B-2 to B-11 (except C1 to C4 of B-9) had
  been seen
- The two-way interval of B-11 B7 (not inside the margin) is reported in the abstract and in the new exploratory section beside
  the pre-specified interval. The one-sided reading (an upper bound on improvement) is labelled post hoc
- Removing the odds ratios near settlements from the abstract was decided in this amendment, after all versions had been seen
- The operating characteristics of B1 describe the size of the rules and do not change any verdict

## Addendum 1 (2026-10-29T10:02Z, commit 13de450, before the quantities below were computed)

Three points emerged while implementing B. Results of B1, B3 to B10, B12 to B17 and B19 to B21 had been seen; B2, B11 and B18
were running (results not seen). The allocation of B21 (countries of the names not read as water-side) showed that, in the
uncapped sample, much of the flooded ground within 500 m of water-side names lies in Cambodia and Myanmar.

1. **Correcting the "near settlements" flag of the census (implementation error)**: B8 of B-11 did not flag estimates whose path
   lacks le500m and similar tokens. Outputs of the Sect. 5 estimator (those carrying n_group_flooded or n_group_cells, as
   returned by v55.odds and similar) are near-settlement estimates whatever their path, so they are flagged by these keys
   too (the B6 variants of v57, the same-cell estimates of B5, the river-size strata and heights above the main river of B7,
   near_500m of B12, v61, the near estimates of C5 and C7 in v64, and so on). The size of the family, the number excluding
   the no-effect value and the BH q values are recounted with the corrected rule; Table S14 and the main text take the
   corrected values and say that they were corrected
2. **Benjamini–Yekutieli adjustment**: in the near-settlement family, report the BY q values, valid under arbitrary dependence,
   beside BH (referee 2, M7; the estimates are strongly dependent variants of the same quantities)
3. **Adding matched forms to the audit sample**: to the 250 names of B21, add 50 names from the form × country strata that supply
   flooded cells within 500 m of water-side names in the uncapped sample and are absent from the B11 sample of B-11 (Cambodia,
   Myanmar and others), by the rule of B11 (up to three names per stratum, from the heaviest strata; weights from the uncapped
   sample's flooded and water-side cells), making 300 names. The answers are not published

## Addendum 2 (2026-10-29T10:30Z, commit dd007b0, before the quantities below were recomputed)

Three points emerged after the results of B19 and B18 had been seen, and before the results of B2, B11 and C1 were seen.

1. **Correcting the label flag of B19 (implementation error)**: B19 decided whether an estimate used GFM labels also from the
   file name, so it marked as GFM the EMS-label estimates in v41, v45 and v50 (EMS labels on the cells that GFM can assess).
   The flag is set from the path: GFM when a path segment starts with odds_gfm, equals gfm or equals odds_v1 (the V1 variant
   of the GFM labels in v50), and for the gfm and gfm_not_s1 rows of B2 of B-11; EMS otherwise. EMS estimates also record
   their domain: the cells that GFM can assess (v41, v45, v50 and the ems_same_cells and ems_not_s1 rows of B2) or all cells
   (capped: v39, v55; uncapped: the ems_all row of B2). The ranges in the text and the supplement figure use the corrected
   flag
2. **B19 reads the corrected census**: B19 is "every water-side odds ratio within 500 m in the census", so it reads the census
   corrected in addendum 1 (v69_estimates.csv)
3. **Numbering of the supplement archive (implementing A10)**: the supplement keeps the tables that the paper and the
   supplement text cite and the tables that those tables cite (so that every reference resolves). Tables moved to the archive
   are numbered A1, A2, ...; references in the archive document use S numbers (supplement) and A numbers

## Addendum 3 (2026-10-29T10:32Z, commit 1881566, before B19 was recomputed)

This emerged after B19 had been computed as corrected in addendum 2, and before the results of B2, B11 and C1 were seen.
B19 took only rows whose path contains wetside and le500m, so it dropped the versions of the Sect. 5 estimator whose
near-settlement flag addendum 1 corrected (their paths lack le500m). This is the error of addendum 1. The values of the
dropped rows had each been reported in earlier amendments, and were seen when the rows were listed.

1. **Rows that B19 takes**: rows of the corrected census whose estimand is the odds ratio of flooding within 500 m of
   water-side names relative to other settlements of the same AOI (one pooled value). (a) Rows whose path names the
   water-side class within 500 m (wetside and le500m). (b) Versions flagged by the keys of the Sect. 5 estimator that
   estimate the water-side class: wetside_adjusted_cci of odds_ems and odds_gfm in v55, partial in v55 (versions of the label
   definition), the B6 versions of v57, odds_adjusted_h_main and odds_same_cells_unadjusted of B7 in v57, the restricted
   samples of v61 (without Pakistan; without Pakistan and Bangladesh), and C5_near and C7_near of v64 (definition of
   permanent water; bare-earth elevation)
2. **Rows not taken**: values by stratum (severity in v55, river size in B7 of v57), the rows leaving out one activation at a
   time (B6 of v63, which has its own table), and other classes (the floodplain reading: v55, B5 of v57, B2 of v66;
   river names: near_500m of B12 in v57)
3. **Duplicates**: rows with the same estimate and 95 % interval to three decimals are counted once (for example uncapped of
   B6 in v57 and ems_all of B2 in B-11). Versions without the cap (uncapped in the path, and B2) are marked
