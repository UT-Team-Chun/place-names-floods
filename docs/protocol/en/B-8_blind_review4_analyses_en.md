# B-8 Corrections of reporting and post hoc analyses after the fourth blind review (NHESS version) (2026-10-27)

English translation drafted with a language model; the Japanese file `docs/protocol/amendments/B-8_blind_review4_analyses.md`
is the record and governs where the two differ.

Status: **plan (2026-10-27, before any quantity of parts B and C was computed)**. The confirmatory verdict is not changed.
Everything below is **post hoc** and responds to an editor and three referees (language-model agents that read only the
PDFs of the manuscript and of its Supplement; reports in `docs/research/reviews/2026-10-27_paper2_blind4_*.md`, commit
6361fa7). The paper calls these "post hoc analyses fixed in this document, after the confirmatory results were known and
before they were computed". The word "registered" is used only for deposit with a third party.
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924). Code: the code with the correction of
B-4 (corrected code).

Values seen when this document was written: all results of the confirmatory analysis and of B-2 to B-7 (the values in
the manuscript version of commit 38afd03), and values that the reports point out or approximate. In particular the
reports give: bounds with the mean of the B-7 reassignments as the no-effect reference of about 1.25 (all ground) and
about 6 (500 m) and a detectable odds ratio of about 1.55 (referee 2, by interpolating Table S29); a 500 m odds-ratio
interval of about 0.6–1.6 after a simple correction for misclassification (referee 3); and a smallest margin at which
equivalence holds of about 0.011 (referee 2). None of the quantities of parts B and C has been computed yet.

## Main points of the reports (all major revision; the editor would return the manuscript before review)

1. The pre-specified test (R pooled over areas) cannot detect the effect that motivates the claim (local to the site of
   a settlement and the ground around it). This was found only after the results (the ceiling of B-6), yet the
   Introduction writes the area-wide estimand as if it had been the design intention. The interpretation plan said it
   would read equivalent-null as a "settled null" on the place-specific information of names
2. The readings of the names (the lexicon) were not validated independently (only the model that built the lexicon
   coded the Wikipedia excerpts). The lexicon departs from its stated policies (two-letter affixes, river names in the
   Tupi block, words of regional languages). The language of the Wikipedia editions was inappropriate for some countries
3. The manuscript reads like an audit trail and is too long. Its provenance (development history, critiques, data)
   cannot be checked by third parties
4. Near settlements: on the same cells both EMS and GFM give odds ratios below one, yet the paper says the label sets
   disagree in sign. The estimates by element class, which are closest to the claim, are not in the main text
5. The bounds and the detectable effect use as no-effect reference the mean of reassignments withdrawn in B-7. The
   500 m plasmode curves carry large Monte Carlo error
6. The reference model M1c is unstable across splits (activation split, refit on GFM). The intervals hold the fitted
   models fixed
7. The abstract and Sect. 5.4 claim more than the evidence. The abstract breaks the reporting rules fixed in B-2
   addendum 3, B-3 and B-6 item 4

## A. Corrections of reporting (no computation; fixed in the paper and the Supplement)

1. **Places that broke our own reporting rules** (referees 1 and 2). The abstract (i) gave the settlement-point positive
   control and the 500 m odds ratio, which were not planned, in the same voice as the confirmatory result without saying
   so, and did not give beside them the planned all-cell positive control (water-side names, pooled, do not mark low
   cells) (rule of B-2 addendum 3); (ii) said nothing about the designated replication on GFM labels (C2) (rule of B-3);
   (iii) gave the GFM odds ratio without the EMS value on the same cells and the exclusion of built-up land and seasonal
   water in the same sentence (rule of B-6 item 4). The abstract and the Conclusions are rewritten to follow the rules.
   The GFM odds ratio is taken out of the abstract (written by the rule it becomes long, and it would conflict with the
   policy of not interpreting a single interval without corroboration)
2. **Error in Sect. 4.2**: it said the verdict was equivalent-null "in the label and sample variants of Table C5", but
   the run with only AOIs whose flood type is recorded was negative (R −0.009, 90 % −0.019 to −0.001), and two rows
   have no verdict. Corrected, and written beside the negative verdict of the 2016-onward sensitivity analysis
3. **Error "the two label sets do not agree on the sign near settlements"** (all): on the same cells EMS gives 0.79 and
   GFM 0.55, both below one. The change from 1.01 to 0.80 is a change of sample (the AOIs that GFM can assess). Text and
   Conclusions corrected (written from the test of B5)
4. **Reading of the recurrent-flooding control** (referee 1): assignment through settlements keeps about a seventh of the
   information (0.012/0.088 overall, 0.081/0.587 within the 227 AOIs with other events), and within those AOIs the point
   estimate assigned through settlements exceeds the margin. The sentence that a name "would have to carry more ... than
   such a summary does" is replaced by: names would need to carry about as much information as a record of past mapped
   floods around every settlement, and far more than recurrent surface water around it (R 0.006). Because "other events"
   include later events, the control tests what the design can register, not a predictor available in advance
5. **A planned diagnostic that was not reported** (referee 2): the SD of the standardised deviations (R − mean R)/SE of the
   synthetic-label replicates (δ = 0), planned in B-2 addendum 3 item 3 (computed, 0.89), and the median bootstrap SE
   (0.0021, equal to the robust SD of 0.0021) are written in the text, and they correct the statement that the SD between
   replicates is 2.1 times the bootstrap SE (zero is covered less often because name features slightly worsen held-out
   Brier scores under the null, mean R −0.0008, not because the intervals are too narrow)
6. **Corrections of the description of the methods** (referee 3): (a) positions: the lexicon has entries only for whole
   words, beginnings and ends; no entry matches inside a word (four letters or more). (b) Unlike its header (no short
   forms), the frozen lexicon contains nine two-letter affixes (German -au, Swedish and Norwegian -as, Norwegian -oy,
   Danish -aa, Turkish and Kazakh -su, Persian -ab, Indonesian ci-). The English suffixes -burn, -bourne, -borne, -pool,
   -well and -wells also read as water names derived from transferred names or personal names (Auburn, Melbourne,
   Liverpool, Cromwell). (c) Unlike "no names of places or rivers", three Tupi entries (parana, iguacu, paranapanema)
   are also the names of large rivers. (d) Unlike "regional languages are not read", some blocks contain words of
   regional languages (Indonesian ci-, krueng, way; Punjabi and Sindhi words in the Persian block). The list of unread
   regions is extended (Catalonia, the Basque Country, Galicia, Ireland, Scottish Gaelic, Sámi, Aboriginal Australia,
   Māori, Alsace-Moselle, French Flanders). (e) "Blocks of related languages" become lexical groupings (not genealogy).
   (f) The denominator of the agreement between the language model and the lexicon is the names of each lexicon class,
   not all names. (g) Wikipedia editions were chosen by "the language of the country's first block", so Persian for
   Pakistan, Ukrainian for Russia and Hindi for Bangladesh; names of these countries were hardly checked. (h) Finnish ala
   and yla mostly mark downstream and upstream along a river
7. **Times** (editor): times written inside the files (e.g. addendum 1 of B-7 at 05:20Z, the interpretation plan at
   06:10 JST) are the writer's notes and differ by minutes from the commits that recorded them (05:16Z, 06:03 JST).
   Table A2 gives commit times and lists every addendum (addenda 2 and 3 of B-7 added). This difference is stated in the
   appendix
8. **Interpretation plan and the limit of the design** (editor, referee 2): the interpretation plan said it would read
   equivalent-null as "the place-specific information of names is below the SESOI (a settled null)". That a test pooled
   over areas cannot detect an effect confined to the ground near settlements was found after the results (the ceiling of
   B-6). The paper narrows this reading to an area-wide effect and says why and when (Introduction, Sect. 2)
9. **Construction of the confirmatory placebo** (referee 2): like the reassignments withdrawn in B-7, the confirmatory
   placebo (M3np) combines reassigned features with the regional composition of the real names. It does not affect the
   primary verdict (the placebo condition is used only when the 95 % interval of R lies above zero), and the paper says
   so. The "80th percentile" placebo is the placebo of the corrected code (R −0.006), and the paper says so
10. **Three intervals for one estimate** (all): the pairs bootstrap of the 500 m EMS odds ratio is unified on the 2,000
    replicates of B-7 (0.77–1.31); the others (1,000 and 200 replicates) go to the Supplement with their counts
11. **Other statements**: the definition of 2.9 % (all cells GFM cannot assess; 2.0 % for cells mostly under the exclusion
    mask); 6.5 % against 6 %; "was right for fewer" becomes a share; the interval level of the HAND control (Table 2
    against Table S19); the reference of the footnote of Table S7; the two-way bootstrap described as multiplied weights
    (the pigeonhole bootstrap of Owen 2007) with its replicate counts (1,000 and 4,000); the 53 units of the country
    split explained in the main text; the dependence of the H4 TOST on an inflated yardstick shown in the table and the
    figure; when the p values were computed stated in Sect. 2.7; the number of post hoc runs whose label differed (since
    "no post hoc analysis changes a verdict" is true by definition); the bound of H3; the GFM HAND threshold cited from
    Sect. 4.4.3 of PDD v4.0 (≥ 10 m); a rewording of "GFM does not keep frequent local flooding"; literature on the
    under-detection of flooding in built-up areas; a source for floods affecting the most people; why overview maps were
    preferred to detail maps; the definition of "flooded share near each name" for N-learn; a source for the training
    cut-off of Gemma; overlaps in Table S23 (-bach also counts as -ach); the language-model prompt in the Supplement; and
    that the three events of the Japanese test were development data for this protocol whose results had been seen

## B. Local computations (from existing artefacts; script v57)

1. **No-effect reference of the bounds and the detectable effect, and Monte Carlo error** (referee 2, editor)
   - The primary no-effect reference ΔBS_null is the mean ΔBS of the 200 reassignments with the composition recomputed
     (B-7 item 1). Beside it: zero, the mean of the reassignments centred within AOIs (B-7), and the mean of the B-5
     reassignments (the values used so far)
   - Curves: the fixed-base plasmode of B-6 (forms A, B, CA, CB) and the refitted one of B-7 (form A). After the
     replicates of C1 are added, the 50-replicate curves are primary and the 10-replicate values go to the Supplement. The
     reading is that of B-6 (bound: the odds ratio at which ΔBS_null + g(δ) reaches the upper 90 % limit of the observed
     increment, interpolating linearly on the δ grid; 80 % power: Φ((ΔBS_null + g(δ))/SE − 1.96) ≥ 0.8)
   - Monte Carlo intervals: replicates are resampled with replacement within each δ (2,000 times, seed 20261028), the mean
     curve is rebuilt, and the odds ratios are read the same way; 2.5th to 97.5th percentiles, and the share of resamples
     that do not reach the target within the grid
   - **Rule**: the paper's bounds and detectable effects use the mean of the rebuilt B-7 reassignments as reference, and
     a table shows how they depend on the reference. They are not called "conservative". A bound whose Monte Carlo 95 %
     interval reaches the top of the grid is written as "not determined within the simulated range"
2. **Achieved equivalence bound** (referee 2): max(|lower|, |upper|) of the 90 % interval of the primary comparison (on
   the R and Brier scales), for the frozen and the corrected code, and its maximum over the combinations of five
   yardsticks and two bootstraps (Table S2)
3. **Location groups needed for inconclusive verdicts** (rule of the interpretation plan; referee 2): for the inconclusive
   comparisons on the primary split (total increment M1s → M3s, flexible refit M1 → M3) and for R within 500 m (Table S6,
   not planned), assuming that interval half-widths scale with the inverse square root of the number of location groups
   and that the point estimate stays, the number needed for (i) the 90 % interval to lie within ±0.05 and (ii) the lower
   95 % limit to exceed zero (and its ratio to the present 195)
4. **Permutation test** (referee 2): the one-sided p value of the real increment in the distribution of the 200 rebuilt
   reassignments of B-7, (1 + #{ΔBS_q ≥ ΔBS_real})/(1 + 200), and the label that the decision table would give with the
   placebo condition replaced by the 95th percentile of this distribution (descriptive)
5. **GFM** (all)
   - (a) **Label difference on the same cells**: on the cells GFM can assess (those of 0.79 and 0.55), the difference of
     the log odds ratios of water-side × within 500 m (GFM − EMS), by a pairs cluster bootstrap over location groups
     (2,000 replicates, seed 20261028, both fitted on the same resample); its 95 % interval and the upper 95 % limits of both
   - (b) **Small-sample correction**: the jackknife variance from estimates leaving out one location group at a time
     (CV3, (G − 1)/G Σ(β_(g) − β̄)²) with a t(G − 1) 95 % interval, for the GFM odds ratio and the EMS odds ratio on the
     same cells
   - (c) The values on the 68 AOIs whose EMS map was not drawn from Sentinel-1 (planned in B-3; in the Supplement) move
     to the main text
   - (d) The EMS odds ratio of the **floodplain reading** on the cells GFM can assess (the cells of the GFM 0.31)
   - (e) **Diagnosis of the failed refit**: for M0, M1, M1s and M1c refitted on GFM labels and for M1 trained on EMS
     labels, on the same cells and GFM labels: (i) the median (quartiles) of the within-AOI AUC; (ii) the Brier score split
     into a between-AOI part (weighted mean squared difference between AOI-mean prediction and AOI flooded share) and a
     within-AOI part; (iii) the distribution of AOI flooded shares on the training and evaluation sides. **Reading rule**:
     if the median within-AOI AUC of M1 refitted on GFM lies within 0.05 of that of M1 trained on EMS, write that the
     refitted model learned the terrain ordering within AOIs but the AOI-level flooded share did not transfer between
     location groups; otherwise write that the relation with terrain could not be learned from the GFM labels
   - (f) The designated comparison ΔBS(M1c → M3n) = +5.7×10⁻⁴ (Table S20) and the Brier score of M1c (0.1155; M1s
     0.0778) go to the main text
6. **Variants of the near-settlement odds ratio** (EMS, within 500 m, the model of B-2's v39: AOI fixed effects + logit of
   the held-out M1 prediction + water-side × distance band; intervals clustered by location group and by a pairs cluster
   bootstrap with 2,000 replicates). All are descriptive and not planned in the protocol
   - (a) Only AOIs whose flood type is recorded
   - (b) Without the 215 AOIs that lack a hydrography layer
   - (c) Without Oceania and the Americas (settler societies), and without Australia alone
   - (d) **Cleaned reading**: water-side re-read with a lexicon without the nine two-letter affixes and the six English
     suffixes of A6(b) (-burn, -bourne, -borne, -pool, -well, -wells)
   - (e) **Uncapped sample**: every evaluable cell of the 500 AOIs (the cells and features tables of each AOI); the
     terrain prediction is the held-out M1 of the uncapped sensitivity run (world_v2_s2a_nocap); names read with the
     same lexicon (corrected code)
   - (f) **Correction for misclassification** (probabilistic bias analysis): with the rare-outcome approximation, the
     odds ratio for truly water-related names OR = (OR*(1 − p0) − (1 − p1))/(p1 − OR* p0); p1 (share of names read as
     water-side with a water-related derivation) and p0 (share of other names with a water-related derivation) are drawn
     from Beta(count + 1, rest + 1) using the counts of the Wikipedia check of B-5, and OR* is drawn from the pairs
     bootstrap distribution of the main estimate (not (e)) (20,000 draws, seed 20261028); median and 2.5th to 97.5th
     percentiles. The paper states beside it that p1 and p0 were measured with the codes of the model that built the
     lexicon
   - (g) By class (existing values) with Holm-adjusted intervals, as the main description in the text. The floodplain
     reading (B-7) is written as "a reading fixed after the class estimates were seen", beside the pooled value (1.01)
   - (h) **River size**: settlements are split by the upstream area (UPLAND_SKM) of their nearest HydroRIVERS reach (< 100
     km², 100–1000, ≥ 1000), with the water-side odds ratio in each stratum; and a version that adds as covariate the
     height above the nearest reach with upstream area ≥ 1000 km² (B7(a) below)
   - (i) Fractional outcome in every distance band: a quasi-binomial (logit) regression with the riverine cover of the
     250 m cell as outcome, giving the water-side ratio within 500 m, 500 m–1 km and beyond 1 km
   - **Rule**: the paper's conclusion about near-settlement flooding (unresolved) changes only if the 95 % interval of (e)
     (pairs) excludes one and (a) and (d) point the same way; even then it is written as a post hoc association without
     corroboration
7. **Further positive controls** (referee 1)
   - (a) **Height above the main river**: for each settlement, the nearest point (within the AOI grid) of HydroRIVERS
     reaches with upstream area ≥ 100 km² and ≥ 1000 km²; height = mean elevation of the settlement's cell − minimum
     elevation of the cell of that point. Within-AOI difference water-side − other (fixed effects, cluster bootstrap over
     location groups, as for the HAND control), and by class
   - (b) Within-AOI difference of the (log) upstream area of the nearest reach
8. **Calibration** (referee 2, editor): for the primary, country, time, language-model cut-off and activation (frozen and
   corrected) splits and the refit on GFM, for the held-out predictions of M1, M1s, M1c and M3n: (i) calibration intercept
   and slope (logistic regression on logit p), (ii) Murphy's decomposition of the Brier score (reliability, resolution,
   uncertainty; ten equal-count bins), (iii) median within-AOI AUC
9. **Skill of the terrain baseline in flood-mapping terms** (referee 1): for the held-out M1 of the primary split and for
   the ranking by HAND alone (cell minimum): median (quartiles) over AOIs of CSI, hit rate and false-alarm ratio when as
   many cells as were observed flooded in each AOI are predicted flooded, and of the within-AOI AUC
10. **Post hoc intervals that exclude the null by chance** (referee 2): the number of post hoc interval estimates in the
    status list of the Supplement (Table S33), how many exclude the no-effect value (0 or 1), and 0.05 times that number
11. **Audit of the lexicon (descriptive)** (referee 3): for each form of A6(b)–(d) (two-letter affixes, the six English
    suffixes, the three Tupi words, words of regional languages) and for hi char: the number of analysed settlements
    matched and of cells nearest them, by country; the most frequent names read through -ab; a table of matched forms by
    country (Supplement)
12. **Language-independent hydronym reading** (referee 3): GeoNames features of class H whose code starts with STM
    (rivers, streams) or LK (lakes), within 5 km of the settlement. Their asciiname is split into lower-case letter tokens;
    a core is a token of four letters or more that is among neither (i) the 200 most frequent tokens of all class-H names
    (generic words such as river, creek, rio) nor (ii) the 200 most frequent tokens of settlement names in that country.
    A settlement name "contains a hydronym" when one of its tokens equals a core, or starts or ends with it leaving two
    letters or more. Reported: the share of settlements with a match, its overlap with water-side, recall among the 139
    settlements recorded on Wikidata as named after a water feature (lexicon: 9/139), precision in the sample of the
    Wikipedia check, the positive control (within-AOI difference of settlement-point HAND), and the odds ratio within
    500 m (EMS, given terrain). All descriptive and not planned in the protocol
13. **A second coder (a language model of another developer)** (editor, referee 3): the 462 Wikipedia excerpts of B-5
    (settlements with a section on the name) are coded by a model of another developer that took no part in the lexicon
    or the first coding (the default model of OpenAI's Codex), with the rules of the guide for human coders
    (v51_audit_guide.md), without seeing the lexicon reading, the stratum or the first codes. Two codes are added to the
    guide: whether the name was transferred from a water feature or is commemorative (yes, no, unknown), and whether what
    the name refers to is at this site (yes, no, unknown). Excerpts are given 20 at a time with a fixed prompt; output is
    JSON. Reported: agreement with the first codes (Cohen's κ by class and for water-related or not), precision and recall
    under the second codes (broad and strict), and the share of transferred names. **Rule**: if precision or recall under
    the second codes differs from the first by 10 points or more, write that the measured rates depend on the coder. The
    paper states that this does not replace coding by human experts

## C. Cluster computations (utens)

1. **More plasmode replicates** (fixed base, forms A, B, CA, CB, the δ grid of B-6, replicates 0–49; replicates 0–9 use
   the seeds of B-6 and must reproduce B-6). Run world_v2_b8plas50
2. **Interval with the second stage refitted inside the bootstrap** (referee 2): the primary comparison of the corrected
   code (M1c → M3n). In each fold the design of the second stage is saved (cross-fitted base logits on the training side,
   composition and name features, the same on the evaluation side, location groups; world_eval --dump-stage2, run
   world_v2_b8stage2); for each of 500 bootstrap replicates (multinomial weights over location groups, seed 20261028 +
   replicate), the second stages of M1c and M3n are refitted in each fold with the training-side weights and Brier scores
   are summed with the evaluation-side weights. Base predictions are fixed. The 90 % and 95 % percentile intervals of
   R_b = ΔBS_b / the point estimate of the yardstick (the fixed-margin form of the decision) are set beside the
   fixed-model intervals. **Rule**: if the refitted 90 % interval lies within ±0.05, write "equivalent-null also with the
   variability of refitting"; if not, write that the fixed-model interval understated the uncertainty (the confirmatory
   label is not changed)
3. **Operating characteristic of the decision at the margin** (referee 2): 40 replicates (100–139) injecting an odds
   ratio of 2.34 (the margin translation of B-6) in form A with every model refitted (the injection rule of B-7 item 2).
   For each replicate the Brier sums of each evaluation-side location group are saved, and the replicate is judged with
   the same cluster bootstrap (4,000) and decision table as the primary comparison. Reported: the share of replicates
   whose 90 % interval lies within ±0.05 (Wilson interval) and the distribution of labels. **Rule**: if the lower limit of
   the Wilson 95 % interval of that share exceeds 0.05, write that the fixed-model interval is anticonservative at the
   margin. Run world_v2_b8margin
4. **Recurrent surface water assigned to each cell directly** (referee 1): JRC GSW v1.4 occurrence, seasonality and maximum
   extent averaged over each 250 m cell of the AOI grid (30 m pixels read at a 5 × 5 set of points in the cell). Features
   = [occurrence, seasonality, maximum extent]; second stage M1cW → M3oW (composition = AOI means of the same features).
   Set beside the settlement-assigned GSW control of B-6 (R 0.006). **Reading rule**: if R assigned to cells exceeds the
   margin and R assigned through settlements does not, write "for GSW too, assignment is the bottleneck"; if both lie
   below the margin, write "recurrent surface water carries little information about the mapped extent once terrain is
   known, however it is assigned". Run world_v2_b8gswcell
5. **Primary comparison with the cleaned lexicon** (referee 3): with the lexicon without the forms of B6(d) (a derived file
   made mechanically from the frozen lexicon; the frozen file is not changed), M1c → M3n with the composition recomputed.
   Run world_v2_b8lexclean

## D. Not done here (author decisions, or a separate study)

- Independent coding by toponymists (two coders per language, κ), extended to English (including Australia and North
  America), German, Urdu/Punjabi, Bengali, Romanian, Swedish and Portuguese, using expert sources (Key to English
  Place-Names, NSW Geographical Names Register and others). The authors will commission it
- A time-stamped archive of the development history (Software Heritage, Zenodo), release of the critiques and of the
  agents' computations, the data deposit (DOI), the authors' check of the translation of the protocol, and whether a
  summary of the Japanese test goes into the Supplement. The authors decide
- A bare-earth DEM (FABDEM) and MERIT Hydro HAND, building footprints (GHSL, WSF), protection standards (FLOPROS),
  administrative units (LAU), detail-scale products, the age of names (Wikidata P1249), and readings of each name in its
  own language. These belong to the separate settlement-scale study (to be registered with a third party in advance) and
  are stated as limitations

## Reporting rules

- All quantities above are reported in the Supplement. The main text carries only those that carry the argument (B1,
  B5(a)(e), the main variants of B6, B7(a), C2, C4 and summaries of B12 and B13)
- The confirmatory verdict is not changed. Where results disagree with earlier readings, this is written as it is
- The order of computation, thresholds and group definitions are not changed after results are seen. Any change is added
  to this document, and the results before the change are also reported

## Addendum 1 (2026-10-27, before any quantity of parts B and C was computed): implementation record

- **Code**: `topsem.world.b8` (with tests); in `world_eval.py` the options `--lexicon`, `--dump-stage2`,
  `--plasref-deltas`, `--plasref-reps`, `--gsw-cell` and `--plasmode-rep-first`; `world_boot_refit.py` (C2);
  `v57_lexicon_clean.py` (the derived cleaned lexicon `docs/verification/v57_lexicon_clean.yaml`; 15 entries removed: nine
  two-letter affixes and six English suffixes; the mapping of countries to languages is that of the frozen lexicon);
  `v58_gsw_cells.py` (per-cell GSW for C4)
- **C1**: replicates 0–9 are the values of the B-6 run (world_v2_b6plasmode); replicates 10–49 run in ranges of ten
  (world_v2_b8plas50, `--plasmode-rep-first`). The salt of a replicate is fixed by its number, so splitting into ranges
  gives the same values. As a check that B-6 is reproduced, fold 0 replicates 0–9 are rerun in the same run; if they do
  not match B-6, everything is rerun
- **C2**: the second-stage design is saved for the primary split with the corrected code (world_v2_b8stage2; arm M0;
  stacks M1c and M3n). The bootstrap is world_v2_b8bootref, 50 jobs of ten replicates (0–499). In each fold the
  unweighted refit must reproduce the saved predictions (the maximum absolute difference is recorded), and each weighted
  refit starts lbfgs from the unweighted coefficients (the point of convergence is the same)
- **C3**: the injected δ is the δ of the B-6 margin translation itself (0.8484278831, odds ratio 2.336). Replicates
  100–139 run as four jobs of ten per fold (world_v2_b8margin)
- **C4**: per-cell GSW is retrieved in the cluster script stage (v58); the evaluation run is world_v2_b8gswcell (arms M0
  and M1; stacks M1c, M3n, M1cW, M3oW)
- **C5**: world_v2_b8lexclean (`--lexicon docs/verification/v57_lexicon_clean.yaml`; arm M0; stacks M1c, M3n, M3np)
- **Resources**: evaluation and design saves 4 CPU and 12 Gi; plasmode ranges 2 CPU and 10 Gi; refitted plasmode 4 CPU
  and 14 Gi; bootstrap 2 CPU and 10 Gi; GSW retrieval 2 CPU and 12 Gi. Nodes taurus and virgo are avoided

## Addendum 2 (2026-10-27, before any result of B13, B5 or B6 was seen): implementation record

- **B13**: `v57b_second_coder.py`. The prompt is the class and status rules of the guide plus transferred (a name taken
  from another place or from a feature not at the settlement, or commemorative), at_site (whether what the name refers to
  is at this site) and spring_or_well_only (whether the only water element is a spring or well, for scoring the strict
  reading); the PROMPT in the file is the record. Coding runs 20 items at a time in an empty temporary folder, and a batch
  is retried up to three times if the output ids do not match the input. As a format check, the first three items were
  coded once to see the shape of the output (after which the wording of transferred was made explicit as "a feature not
  at the settlement"); the final coding codes everything again, these three items included
- **Implementation of B6 and B7**: from the nine regional HydroRIVERS files, every vertex of every reach that touches the
  extent of an analysed AOI (± 0.2°); nearest points via a KD tree on unit-sphere coordinates. For the class adjustment
  (B6(g)) Holm-adjusted p values are given with Bonferroni 95 % intervals, because Holm's procedure does not define
  simultaneous intervals. In the uncapped sample (B6(e)) the readings are the water-side flag of the nearest name in the
  uncapped sensitivity run (frozen code), and the distance to the nearest settlement is joined from each AOI's feature
  table by coordinates (within 1e-5°). Every pairs bootstrap uses 2,000 replicates
- **B5(e)**: the predictions refitted on GFM labels (world_v2_gfm) and the EMS-trained M1 (frozen run) are joined by AOI
  and coordinates

## Addendum 3 (about 2026-10-27T10:20Z, before any result of C1–C5 was seen): failure of the cluster NAS and change of the computing environment

- **Failure**: the NAS of utens (the GlusterFS FUSE mount `/mnt/nas`) disconnected on every node ("Transport endpoint is
  not connected", input/output errors, files not visible). The second-stage design saves of C2 (world_v2_b8stage2)
  completed in all ten folds, but some files on the NAS cannot be read. The jobs of C1, C3, C4 and C5 failed reading their
  inputs (nothing was computed). Most failures were due to the NAS; C5 also failed because the derived lexicon was not in
  the image (.dockerignore excludes docs/verification and the Dockerfile copies only the v44 JSON). The remaining jobs were
  suspended. Remounting the NAS is an operation on shared infrastructure and is left to the authors
- **Change**: C1–C5 run on the local Mac (arm64, 16 cores) with the same code (commit 53b1fb5 or later; `TOPSEM_DATA_DIR`
  is the local data/). Because the last digits of floating-point results may differ from the confirmatory runs on the
  cluster (amd64), local runs carry the suffix `_local`, and the maximum absolute difference of predictions from the
  corrected run (world_v2_lexfix, cluster) is recorded
- **Reduced scope** (local computing time): C1 runs only forms B and CB (500 m, the curves with large Monte Carlo error)
  for replicates 10–49 (forms A and CA keep their ten replicates). C3 is reduced to the 20 replicates 100–119 and runs last
  (if it does not finish in time, the replicates completed are reported). The C2 bootstrap keeps 500 replicates

## Addendum 4 (about 2026-10-27T10:45Z, before any result of C1–C5 was seen and before the 50-replicate curves of forms B and CB were read): changes to the plan and the tolerance for the reproduction check of C1

- **B10 is not counted**: the status list of the Supplement (Table S33) lists groups of estimates and does not define a
  unit for counting single intervals. Deciding after the results what counts as one interval, and whether its no-effect
  value is 0 or 1, would introduce choices, so the B10 count is not computed. The paper states only that the post hoc
  analyses produced several hundred intervals, of which about one in twenty can exclude the no-effect value by chance
  alone, and repeats the existing policy of not interpreting a single interval without corroboration from another
  analysis or another set of labels
- **Tolerance for the reproduction check of C1**: floating-point rounding may differ between the local machine (arm64)
  and the cluster (amd64) (addendum 3). For replicates 0–9 of fold 0 (forms B and CB), the check passes if, for every form
  and δ, the absolute difference between the mean over the ten replicates of the locally rerun ΔBS(δ, replicate) and the
  mean of the B-6 values is less than 0.1 times the standard error of the ten B-6 values of fold 0 (between-replicate SD
  / √10); the 50-replicate curves then combine replicates 0–9 of B-6 with local replicates 10–49. Otherwise replicates
  0–9 (forms B and CB) of every fold are rerun locally and the curves of forms B and CB use only local replicates 0–49
  (the ten-replicate B-6 values are listed in the Supplement)
- **Implementation of B1** (v59, written before the results were seen): Monte Carlo resamples draw replicates with
  replacement independently for each form and δ, from a single generator (seed 20261028) in the order reference, form,
  quantity. The refitted curve of amendment 7 (form A, three replicates) is read the same way, and the paper says that
  intervals from resampling three replicates are coarse. The interval of the odds ratio detectable with 80 % power comes
  from the same resamples
- **Verdicts of C3**: the injected runs have no placebo arm, so the decision table is applied without the placebo
  condition (not-meaning cannot occur). The report is the share of replicates whose 90 % interval lies within ±0.05
  (Wilson 95 % interval) and the distribution of labels, as planned
- **"Exceeds the margin" in the reading rule of C4**: a point estimate of R of at least 0.05 (the point-estimate condition
  of positive-meaningful in the decision table). The settlement-level value is the GSW control at the settlement points
  of amendment 6 (v52, R 0.006)
- **Yardstick of C2**: computed from the M0 and M1 predictions saved locally and listed beside the value of the corrected
  cluster run (world_v2_lexfix)

## Addendum 5 (about 2026-10-27T11:15Z, before any result of C1 or C3 was seen): the NAS is back, and C1 and C3 return to the cluster with their original scope

- **NAS restored**: at about 2026-10-27T11:12Z, read-only check pods (on cancer and gemini, with the same hostPath as the
  jobs) could read the NAS files (we did not remount it)
- **C1 and C3 on the cluster as originally planned**: locally a single C1 job (one fold, one range of ten replicates)
  takes more than an hour, so the reduction of scope in addendum 3 is withdrawn and the scope of the plan and addendum 1
  is restored. C1 runs forms A, B, CA and CB for replicates 10–49 (world_v2_b8plas50, cluster amd64, the environment of
  B-6), with the reproduction check on replicates 0–9 of fold 0. C3 runs replicates 100–139 (world_v2_b8margin). The
  suspended jobs are resumed, and the jobs that failed on NAS errors (C1: fold 0 replicates 0, 10 and 20, fold 5
  replicates 20) are recreated with the same image and arguments. The local C1 and C3 jobs in the queue are cancelled
  before any output (the local rerun of replicates 0–9 of fold 0 is allowed to finish and is reported as a check of
  reproducibility between the local machine and the cluster)
- **B1**: all four forms use 50-replicate curves (B-6 replicates 0–9 and cluster replicates 10–49). The reproduction rule
  of addendum 4 (its tolerance) applies to the cluster rerun of replicates 0–9 of fold 0. C3 reports the share among 40
  replicates with a Wilson interval
- **C2, C4 and C5** stay local (their designs are saved or their runs are done or running)
- If the NAS fails again and ranges remain unfinished, those ranges fall back to the local plan of addendum 3 (C1: forms
  B and CB, replicates 10–49; C3: replicates 100–119), and the paper reports which replicates ran where

## Addendum 6 (about 2026-10-27T11:20Z, before the estimates below were computed): an editor's request missing from the plan (the Pakistan stratum), and where the flooded cells near water-side names lie

- **Omission**: the editor's minor point 17 (the Pakistan stratum ur in Table S8 has 31.5 % flooded cells and R −0.051;
  state which activations dominate it and show the estimates near settlements with them left out) was not in the plan of
  amendment 8
- **Description (v60, computed)**: with the table and AOI selection of the odds ratios near settlements, the cells within
  500 m of settlements with water-side names that flooded on the EMS labels number 361, in 27 countries and 48
  activations; Pakistan holds 90 and Bangladesh 87, about half, and the largest activation (EMSR097, Bangladesh 2014)
  holds 87. These are counts, not estimates; none of the estimates below had been seen
- **Computations (v61, all descriptive and not planned in the protocol)**: the main odds ratio of the water-side class
  within 500 m (the same model as the main estimate of B6, with cluster-robust intervals and a pairs bootstrap of 2,000
  replicates, seed 20261028), (a) without the AOIs in Pakistan and (b) without those in Pakistan and Bangladesh; and
  (c) a description of the Pakistan stratum: evaluable cells, flooded cells, and cells and flooded cells near water-side
  names, by activation
- **Rule**: reported as descriptions in the Supplement, with (b) in one sentence of the main text. The rule that would
  change the conclusion near settlements (unresolved) stays that of B6 (the interval of the uncapped sample excludes one,
  with the typed-AOI and cleaned-reading estimates in the same direction). If an estimate without these countries
  excludes one, it is still described as an uncorroborated estimate after a choice of regions

## Addendum 7 (about 2026-10-27T11:26Z, before the quantity below was computed): referee 3's minor point 10 (precision of the strict reading by language), missing from the plan

- Referee 3's minor point 10 (the Slovene precision of 94 % may reflect spodnji and dol, which Wikipedia confirms
  trivially; report the precision of the strict reading by language) was not in the plan of amendment 8
- **Computation (v62, descriptive)**: with the existing codes of the first coder (v44), the precision of the strict
  reading (amendment 6, item 5: water without springs and wells, wetland, riverine) for each language with at least 30
  coded names (Wilson 95 % intervals), listed in the Supplement beside the precision of the broad reading by language

## Addendum 8 (about 2026-10-28T15:47Z, before the results of C1 were fetched or aggregated): one C1 job failed to save, and is rerun on the cluster

- Of the 41 C1 jobs, the one for fold 1, replicates 20–29, failed when saving its output after computing it
  (2026-10-27T14:30Z; the NAS disconnected briefly: `Transport endpoint is not connected`). The other 40 C1 jobs and the
  40 C3 jobs completed afterwards, and the NAS is back
- Addendum 5 sends ranges left unfinished by a further NAS failure back to the local plan, but the local machine does not
  reproduce the plasmode values within the tolerance (the local C1 check: up to 1.4 standard errors), so this job is rerun
  on the cluster with the same image and arguments. The outputs of C1 and C3 have not yet been fetched or aggregated (the
  reproduction check on replicates 0–9 of fold 0, and C3, may be aggregated without waiting for the rerun: neither
  includes the range rerun)

## Addendum 9 (about 2026-10-28T16:03Z, before the manuscript and tables were changed; no new computation): the refitted-model bound and margin translation are biased low by the coarse grid — removed from the tables and text, and reported as a ratio of gains

- **How it was found**: in a read-through of the manuscript we noticed that the upper bound "with every model refitted" in
  Table 3, 1.16 (B1's reading of the amendment-7 refit, form A, three replicates, read the same way), is smaller than the
  bound with the base models fixed, 1.26. Refitting reduces the gain of the same injection to 0.91–0.95 of the
  fixed-model gain (amendment 7, item 2), so the bound should be larger, not smaller
- **Cause**: the grid of the amendment-7 refit is δ = 0.5, 1.0, 1.5, with no point below 0.5. The bound (δ about 0.15) is
  read on the straight line from 0 to 0.5, but the gain curve is convex (on the amendment-6 grid the gain at δ = 0.25,
  4.0 × 10⁻⁵, is 0.58 of the value of the straight line from 0 to 0.5, 6.8 × 10⁻⁵), so the line overstates the gain and
  understates the bound. The amendment-6 grid with the base models fixed has δ = 0.25, near the bound (δ = 0.23), so the
  bias there is small. For the same reason the odds ratio at the margin with refitted models (2.34, on the line between
  δ = 0.5 and 1.0) is biased low; this is why it equals the fixed-model value although the gain is 0.95 of it
- **Rule (no new computation)**: the refitted-model bounds (B1's refit_A, all references) and the refitted-model odds ratio
  at the margin are removed from Table 3, Table S25, the supplementary table of the refit and the text. Instead we report
  that the refitted gain was 0.91–0.95 of the fixed-model gain for the same injections (three effect sizes, three
  replicates each), so translations and bounds with refitted models would be slightly larger than those with the base
  models fixed. The range of margin translations (2.3–4.1) is set by the fixed-model and synthetic-label values (it is the
  same with or without the refitted reading). The injection of C3 (odds ratio 2.34) is the fixed-model translation of
  amendment 6 and is unaffected
- The research log's account of amendment 7 ("the margin corresponds to an odds ratio of 2.34 also with refitting") and of
  amendment 8 are annotated with this bias
