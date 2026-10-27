# B-6 Post hoc analyses in response to the second blind review (NHESS version, 2026-10-27) — English translation

> English translation (2026-10-27) of `docs/protocol/amendments/B-6_blind_review2_analyses.md`. The Japanese file is
> the record; where the two differ, the Japanese file governs. In the paper this amendment is called amendment 6.

Status: **planned (2026-10-27, before any of the analyses below was computed)**. The confirmatory verdicts are not
changed. Everything below is **post hoc**, added in response to an editor and three referees (language-model agents)
who read only the PDF of the manuscript for the first time (reports in
`docs/research/reviews/2026-10-27_paper2_blind2_*.md`). The paper describes these as "post hoc analyses fixed in this
document after the confirmatory results were seen and before they were computed", distinguished from registration
with a third party. The word "registered" is used only for registration with a third party.
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924).
Code: the code with the B-4 correction (the Norwegian lexicon). Comparisons are set beside the B-4 rerun
(world_v2_lexfix) and, where needed, beside the run with the frozen code (world_v2).

Values seen when this document was written: all results of the confirmatory analysis and of B-2 to B-5 (those in the
paper version of commit 4ac0124), and values that the referees estimated in their reports (for example, that a bound
anchored on the permutation null would be an odds ratio of about 1.6–2.5, and that the ceiling of a local effect is
R ≈ 0.035). None of the quantities below had been computed.

## Analyses (the numbers are used in reporting)

### 1. Ceiling of a local effect (analytical, from existing predictions)

The upper limit of what the primary estimand can register of an effect confined to the ground near settlements. If
the names of a set S of cells were perfectly right (the predictions of the cells of S equal their labels), the
numerator of the primary comparison could grow at most by ΔBS_max(S) = (1/N) Σ_{i∈S} (p_i − y_i)² (p is the held-out
prediction of M1c, N the number of cells), and R_max(S) = ΔBS_max(S) / yardstick. S is (a) the cells whose nearest
settlement's name is water-side and that lie within 500 m of the settlement, (b) the same within 250 m, (c) within
1 km, and (d) all cells whose nearest name is water-side. Both for the confirmatory run (world_v2) and the corrected run
(world_v2_lexfix). **Interpretation rule**: if R_max(a) is below the margin of 0.05, write that "an effect confined to
the ground near settlements (500 m), however strong, could not be judged meaningful by the confirmatory test".

### 2. Bounds recomputed by injecting effects into the observed labels (plasmode)

The three referees and the editor pointed out that the test inversion of B-5 item 4 (the 5th percentile of the world
of synthetic labels) does not match the real-data null (the permutation distribution) in location or spread. We
therefore inject effects into the observed labels, refit only the second stage, and measure the increment.

- **Injection**: in each fold of the primary split, the base models (M0, M1) are fitted once as in the confirmatory
  analysis (the same mechanism as B-5 item 1; the cross-fitted logit of M1 is held fixed), and the cells of the target
  whose label is 0 are turned to 1 with probability q = (p′ − p)/(1 − p), where p is the M1 prediction of the cell
  (cross-fitted on the training side, held-out on the evaluation side) and p′ = p e^δ / (1 − p + p e^δ). The injection
  is applied to both the training and the evaluation side (a world in which names really have an effect). The character
  model (N-learn) is refitted on the injected training labels.
- **Target forms**: (A) all cells whose nearest settlement's name is water-side (the lexicon's reading); (B) those of
  them within 500 m of the settlement; (C) forms with misreading: for each settlement a latent indicator L of a "truly
  water-related name" is drawn, with L = 1 with probability π₁ for settlements that the lexicon read as water-side and
  π₀ otherwise (π₁ = the true rate in the water-side stratum of v44, π₀ = the true rate in the other stratum, both
  unweighted). The targets are all cells nearest a settlement with L = 1 (C-A) and those within 500 m (C-B). The name
  features remain the lexicon's reading.
- **Grid**: δ ∈ {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0} (odds ratios 1.28–20), 10 repetitions for each δ and form (random
  seeds determined by form, δ and repetition). The Brier scores of M1c and M3n are recorded in each fold, and ΔBS(δ) is
  computed over all folds.
- **Gain**: gain(δ) = mean over repetitions of ΔBS(δ) − ΔBS(0) (ΔBS(0) is the increment of the same run without
  injection = the confirmatory value of the corrected run), with the standard error between repetitions.
- **Bound (primary)**: anchored on the real-data null. The null increment ΔBS_null = the mean ΔBS of the 200
  permutations of B-5 item 1; U90 = the upper limit of the 90% interval of the observed increment (one-sided 95%).
  δ_A = the largest δ with ΔBS_null + gain(δ) ≤ U90 (linear interpolation in the log odds ratio between grid points;
  "> 20" beyond the end of the grid). Reported as the odds ratio e^δ_A.
- **Bounds (secondary)**: (i) δ_B anchored on zero (gain(δ) ≤ U90); (ii) the same computation with U95, the upper limit
  of the 95% interval.
- **Odds ratio of the margin**: δ_M with gain(δ) = margin (0.05 × yardstick; what size of effect the margin of the
  primary comparison corresponds to), and δ_M′ with the null anchor added (ΔBS_null + gain(δ) = margin).
- **Reporting rule**: the bounds and the translation of the margin in the main text are replaced by the values of this
  computation (the values of B-5 item 4 stay in the appendix, with a note that the spread of the synthetic world is
  smaller than the uncertainty of the real data). The values of form C are written as bounds for "truly water-related
  names", not for "the lexicon's class". The approximation of holding the base models fixed (the effect of names does
  not enter the base models) is stated. For forms whose bound exceeds the grid, write that "this design cannot bound an
  effect of that form from above".

### 3. Controls that carry flood information (positive controls of what the design can register)

The oracle of B-5 item 5 (the terrain of the settlement site) carries only terrain information and is not an upper
bound for the case in which names carry flood information that terrain lacks (as the referees pointed out). We
therefore give each settlement a feature that carries flood information, assign it to cells as names are assigned
(from the nearest settlement to the cell), and compute the increment of the same form as M1c → M3n on the primary split
(the regional-composition term = the AOI mean of that feature).

- **(a) Ceiling of the event's labels (an upper bound that leaks by construction)**: the flooded share of the
  EMS-evaluable cells within 500 m of the nearest settlement's point (the labels of that event, leaving out the cell
  itself). Settlements with no evaluable cell within 500 m get the AOI median. Because the evaluation-side feature is
  built from evaluation-side labels, it leaks as a prediction; it is written as the upper bound "if names knew the
  flooding around the settlement in that event perfectly".
- **(b) Memory of recurrent water (no leakage)**: the mean, over the 30 m pixels within 500 m of the settlement point,
  of the occurrence (share of time with water, %) and seasonality (number of months with water) of JRC Global Surface
  Water v1.4 (1984–2021, Pekel et al. 2016), and the same means excluding permanent water (occurrence ≥ 90 and
  seasonality = 12). Retrieved per AOI by range reads.
- The comparison is, for (a) and (b) each, an M1c-type reference → the model with the feature added; R and the label
  of the decision table of §6.
- **Interpretation rule**: if the R of (a) does not reach the margin, write that "this design cannot register even
  information that knows the event's flooding around settlements perfectly as an increment that reaches the margin",
  and that the confirmatory null says nothing about the settlement-scale claim. If it does reach the margin, write that
  "the design can register settlement-scale flood information; the null for names shows that names do not carry that
  information at this level". (b) is written in the same way as the increment of "information on recurrent water".
- **(c) Positive control with recurrent water (no flood labels)**: the occurrence and seasonality of (b) (within 500 m
  of the settlement point), compared in the same form as the positive control of B-2 (the HAND of the settlement
  site): the within-AOI difference between settlements with water-side names and other settlements (cluster bootstrap
  over location groups, 4,000 replicates). Also by lexicon class and for the strict reading of item 5 below.

### 4. Diagnostics of the GFM labels, and a variant that restores seasonal water

- **Layer components**: rebuild the GFM labels with the same rules and keep, for each cell, the shares of exclusion
  mask, permanent water and seasonal water (reference water = 2) (the exclusion layer is 0/1 and does not distinguish
  causes, so causes are not separated). For the EMS-evaluable cells, the share that GFM cannot assess, by component,
  overall, within 500 m of settlements with water-side names and within 500 m of other settlements. The EMS flooded
  share of the cells that GFM cannot assess. Among excluded cells, the share with HAND (cell mean, this study) of 10 m
  or more and the shares with high CCI built-up and tree cover (approximations of GFM's exclusion causes).
- **Variant V1 (seasonal water assessable)**: assessable = observed water extent (ensemble_water_extent) 0/1, exclusion
  = 0 and reference water ≠ 1 (not permanent water); flooded = assessable and observed water extent = 1 (water observed
  on seasonal water counts as flooded). Images, windows and thresholds (v ≥ 0.5, f ≥ 0.5) as in B-3. The AOIs are the
  192 adopted in B-3. On V1, compute B (odds ratios of the water-side group within 500 m and over all cells, adjusted
  for terrain) and C1 (rescoring) of B-3.
- **Separating the selection of AOIs from the selection of cells**: the EMS odds ratio (water-side group within 500 m,
  adjusted for terrain) on all EMS-evaluable cells of the 192 AOIs (whether or not GFM can assess them).
- **Number of contributing units**: the number of AOIs and places with flooded cells in the water-side group within
  500 m, and the number of flooded cells (EMS and GFM).
- **Robustness of the intervals**: for the odds ratios of the water-side group within 500 m (GFM, EMS on the same cells,
  EMS on all AOIs), the interval of a pairs cluster bootstrap over location groups (1,000 replicates, refitting the
  regression) and the range when leaving out one place at a time.
- **Diagnosis of C2 (refitting)**: constant predictions (the training-side and the evaluation-side flooded share), the
  Brier scores of M0, M1, M1s, M1c, M3n and M3d, and the increment M1s → M3d (a comparison without the
  regional-composition terms). From the existing output of world_v2_gfm.
- **Agreement**: CSI (agreement on flooding / (either flooded)), P(GFM = 1 | EMS = 0) and P(EMS = 1 | GFM = 0) (overall
  and by post-event image type; planned in B-3), and the rescoring restricted to AOIs whose EMS post-event image is not
  Sentinel-1 (computed in v41) are reported.
- **Reporting rule**: if the abstract gives the GFM odds ratio, the EMS value on the same cells and the exclusion of
  built-up land and seasonal water are written in the same sentence. If the odds ratio of the water-side group moves
  towards 1 under V1, write that "the GFM result depends on masking".

### 5. Documentation of the name readings, and a strict reading

- **Lexicon table**: for each block, the languages it covers (YAML notes), the number of entries per class, the
  countries mapped to it, and the matching rules (word, prefix, suffix, infix); the forms excluded because of
  collisions (those in the YAML notes).
- **Match statistics**: the share of analysed cells whose nearest name matched any class or a water-side class, by
  block and by continent; the 20 water-side forms with the most matches (form, position, number of cells and of
  settlements).
- **Precision by match type**: among the v44 codes (with a derivation), for names that the lexicon read as water-side,
  the precision by the position of the matched form (word, prefix, suffix, infix) (unweighted, Wilson intervals).
- **Strict definitions**:
  - Strict reading of the lexicon = the classes water (without the forms for springs, wells and watering places),
    wetland and riverine. The forms for springs, wells and watering places are those in the table below (fixed in this
    document before computation, like the classes; decided from the meaning of the words only and applied to no name).
  - Strict recorded derivation = the coded classes water (except those that denote only springs or wells), wetland and
    riverine. Among water codes, those whose quoted evidence contains a word for springs or wells (the words of the table
    below and spring, well, fountain, source, Quelle, Brunnen, fontaine, puits, pozzo, fonte, fuente, manantial,
    studnia, studna, studenec, izvir, izvor, kút, forrás, bunar, ain, bir) and no word for rivers, streams, lakes, ponds
    or marshes (river, stream, brook, creek, lake, pond, marsh, Fluss, Bach, See, rivière, ruisseau, lac, étang, rzeka,
    potok, jezioro, reka, jezero, folyó, patak, tó) are taken as "springs or wells only".
  - Precision and recall of the strict reading (with the same stratum weights as v44), and precision and recall by
    recorded class (unweighted, with counts).
- **Examples of names**: among the v44 codes, for each water-side class of the lexicon, the two lowest-numbered names
  whose recorded derivation matches the class and the lowest-numbered one whose derivation does not (name, country,
  reading, recorded derivation, URL). The settlement names of the example AOI of Fig. 3 and their readings.
- **Primary comparison with the strict reading** (cluster): the primary comparison (the M1c → M3n form, primary split)
  with the name features restricted to the classes of the strict reading (water bodies, wetland, riverine, and their
  union).
- **Preparation of a human check**: a random sample of 100 of the v44 codes, stratified by stratum and language (seed
  20261028), and a guide that states the coding rules. The coding by a person (preferably a toponymist) and the
  agreement (κ) are requested by the author and added when available. Not done in this version (stated as a limitation
  in the paper).

### 6. Sensitivity runs (cluster, primary split, corrected code)

- (a) **Feature codes of the gazetteer**: remove PPLX (section of populated place), PPLH (historical), PPLQ (abandoned)
  and PPLW (destroyed) from the candidates for the nearest settlement, and merge duplicates with the same asciiname
  within 100 m.
- (b) **A stricter rung without elevation**: the within-area increment of names on T0− (position, the remaining land
  cover, settlement layout), i.e. T0 of the ladder without the two HydroRIVERS distances and without the CCI water and
  wetland fractions.
- (c) **Excluding AOIs whose products lack a permanent-water layer**: exclude the AOIs whose EMS product has no
  hydrography layer and whose permanent-water share was set to 0.
- (d) **Excluding AOIs without a flood type**: exclude the AOIs treated as riverine because their flooding has no type
  and they lie more than 10 km from the coast.
- (e) **The language-model reading (N-read) on the primary split**: M1rc → M3rn (the same arms as on llmcut). State in
  the same paragraph that a language model may remember events at the place.
- (f) **Permutation distribution on the activation-grouped split**: 200 within-AOI permutations of names (the same
  mechanism as B-5 item 1) on the same split as world_v2_actgrp of B-2. Re-read the not-meaning verdict against the
  permutation distribution (the rank of the real increment) rather than a single permutation.

### 7. Further summaries from existing predictions (local)

- **AOIs weighted equally**: R built from the mean over AOIs of the differences in Brier score (each AOI weighted
  equally), with the interval of a cluster bootstrap over location groups.
- **Excluding settlements with rounded coordinates**: the near-settlement odds ratios (EMS and GFM, the water-side group
  and by class), excluding the cells nearest settlements whose GeoNames coordinates are whole arc-minutes (up to
  ±0.9 km).
- **Reporting the planned power and coverage results**: B-2 item 6 (coverage of the intervals on synthetic labels with
  δ = 0) and the recomputation of its addendum 3 are written in the main text (already computed, null_coverage of v37),
  together with the distribution of verdicts at δ = 0 (equivalent-null 91/100, negative 7/100, etc.).

## Reporting rules

- All quantities above are reported (appendix or supplement). The main text gives summaries of 1, 2 and 3, the
  diagnostics of 4 and the results of V1, and the precision and examples of 5.
- The confirmatory verdicts are not changed. Where the results above disagree with the confirmatory analysis or with the
  readings of B-3 to B-5, this is written as it is.
- Scope of claims: the abstract and conclusions are limited to "within mapped flood extents, the settlement names as
  read here add nothing to 250 m terrain-based predictions of flood extent". The settlement-scale claim is written,
  according to the results of 1 and 3, as "this design cannot answer it" or "it answered it", and is not asserted
  untested. No claim about risk communication is written, because this study did not examine it.
- The order of computation, thresholds and group definitions are not changed after the results are seen. Any change is
  added to this document, and the results before the change are reported too.

## Addendum 1 (2026-10-27, before computing items 8–15 below): implementation record and descriptions requested in the reports

B-6 quantities already computed when this addendum was written: 1 (ceiling), 7 (AOI weighting, rounded coordinates,
coverage), 3(c) (positive control with recurrent water), all of 4 except V1 and the layer components, the
documentation and the precision of the strict reading in 5, and the primary comparison with the strict reading (cluster)
in 6. The other cluster runs were in progress, and their results had not been seen.

### Implementation record (the computations are not changed)

- Operational definitions fixed during implementation (each before the results were seen): in 4, an "excluded cell" is
  a cell that GFM cannot assess and whose exclusion share is at least 0.5; "high CCI built-up or tree cover" is a share
  of at least 0.5 (means also given). In 5, "a word appears" means that the quotation, lower-cased and split into
  letter sequences, contains the word of the table or the word with s or es added. In 5, the sample for the human check
  is allocated in proportion to the counts of stratum × language (largest-remainder rounding).
- In the pairs bootstrap of item 4, places drawn more than once are refitted as distinct places and distinct AOIs.
- Failures and resubmissions of cluster runs: (i) the retrieval of GSW for 3(b) stopped because the container has no
  CLAUDE.md and could not find the repository root (fixed in commit 4dc4ecf); (ii) the plasmode of item 2 stopped
  because the input for π₁ and π₀ (the v44 JSON) was not in the image (fixed in fd07605). Both stopped before computing
  and were resubmitted. (iii) The T0− arms of 6(b) were not run in the first run (world_v2_b6t0x), because the check
  whether to run an arm also required the N-learn feature, which is built inside each fold, to be in the table (fixed in
  f516c6a). The new run is world_v2_b6t0xb. The T0 and primary-comparison values of the first run are kept.

### Further descriptions (items 8–15; the confirmatory verdicts are not changed)

8. **Minimum detectable odds ratio** (continuation of item 2): for each form, the δ at which power(δ) = Φ((ΔBS_null +
   gain(δ)) / SE − 1.96) equals 0.8 (linear interpolation between grid points). SE is the bootstrap standard error of
   the confirmatory increment of the corrected run (width of the 95% interval / 3.92), and ΔBS_null the permutation mean
   of item 2. A normal approximation of "a positive verdict requires the lower limit of the 95% interval to exceed 0";
   the placebo condition is not included.
9. **Reference points for the baseline** (primary split, held-out predictions of the corrected run): (a) the CSI of a map
   that, in each AOI, flags as flooded the share of cells equal to the observed flooded share, in increasing order of
   HAND (cell mean) (an optimistic reference, because the threshold uses the labels); (b) the CSI of the predictions of
   M1 and M3n thresholded by the same rule; (c) the Brier score of a constant prediction (the training-side flooded share
   of each fold).
10. **Spatial correlation of residuals**: for the held-out residuals of M1 (y − p), the correlation between pairs of cells
    within AOIs by distance band (0–0.5, 0.5–1, 1–2, 2–5, 5–10, 10–20 km). Up to 1,000 cells drawn at random per AOI
    (seed 20261028). A description of the 5 km buffer.
11. **Regional languages that are not read**: the number and share of cells whose nearest settlement lies in Wales (GB),
    Catalonia, the Basque Country, Galicia, Valencia, the Balearic Islands or Navarre (ES), or Brittany or Corsica (FR)
    (by GeoNames admin1; the codes are checked against the admin1 of well-known towns). The number and share of cells in
    IE, NZ and AU (countries where names of indigenous or minority-language origin are read with the English lexicon).
12. **Riverside constructions**: the share of cells and settlements whose nearest name contains the word sur, upon, am,
    an (followed by der, dem or den), nad, nod, aan or on (only when followed by the) (a literal count, regardless of
    language).
13. **Description of the language-model reading**: the shares of the classes of the nearest settlements' readings
    (including unknown and parse failures), whether the confidence was used in decisions, and by continent the shares of
    cells whose nearest name the language model and the lexicon read as water-side.
14. **Inclusion flow of AOIs**: 999 AOIs → 500, the change in the number of activations, counts by criterion, and the
    characteristics of the excluded AOIs (flood type, flooded share, number of evaluable cells, post-event image type,
    continent).
15. **The p values required by the protocol (not reported so far)**: the TOST p values of H2 and H4 (inverting the
    percentile bootstrap: p = max(P*(ΔBS* ≥ m), P*(ΔBS* ≤ −m)), m the fixed margin), the two predictions of H3 (the
    one-sided bootstrap p of Δ0 − Δ4 and the exact one-sided p of Spearman's ρ), the one-sided bootstrap p of an
    improvement for H4 (P*(ΔBS* ≤ 0)), and the Holm adjustment of H3 and H4 (H3 represented by the p of Δ0 − Δ4). For
    both the frozen run (world_v2) and the corrected run (world_v2_lexfix).

## Table: forms taken as springs, wells and watering places (removed from the class "water" in the strict reading of item 5)

| Block | Forms |
|---|---|
| en | well, wells, spring, springs, soak, waterhole |
| de | born, brunn, brunnen, quell, quelle |
| fr | fontaine, fontaines, source, sources, puits |
| it | fontana, fontane, fonte, pozzo, pozzi, sorgente |
| es | fuente, fuentes, manantial, pozo, pozos, ojo |
| pt | fonte, fontes, nascente, poco, cacimba |
| el | krini, pigi, pigadi, pigadia, vris, vrises, vrisi, vrys, vryses, vrysi, kefalovryso |
| sh | bunar, izvor, izvori, studenac, vrelo, vrela, vrutak |
| sl | izvir, studenec, studenci |
| bg | bunar, cheshma, izvor, izvori, kladenets, kladenci |
| sq | burim, burimi, pus, puse |
| ro | fantana, fantanele, izvor, izvorul, izvoare |
| hu | forras, kut |
| cs | pramen, prameny, studna, studnice, studnicka |
| pl | studnia, studzianka, zrodla, zrodlo |
| uk | dzherelo, istok, kolodez, kolodiazi, krinitsa, krynytsia |
| lt | saltinis, saltiniai |
| lv | avots, avoti |
| et | allikas |
| fi | lahde |
| sv | kalla, kallan |
| no | kilde |
| da | kilde |
| tr | bulak, kaynak, kuyu, pinar |
| ar | ain, ayn, aioun, ayoun, bir, biar, abar, hassi |
| fa | chah, chahi, chashma, chashmeh |
| vi | gieng |
| id | sumur |
| tl | bukal |
| sw | chemchemi, kisima |
| bantu | chitsime |
| af | fontein, put |
| mg | loharano |
| qu | pukyu, puquio |
| nah | ameyal |
| ka | tskaro, tskarosi |
| kk | bulak, bulaq |
| so | ceel |
| am | minch, minchi |

Blocks not in the table (nl, hi, my, th, km, zh, ng, tupi) are taken to have no water form that denotes only a spring or a well.

## Addendum 2 (2026-10-27, before the results of items 2, 3(a), 3(b), 6(b)–(d) and 6(f) were seen): supplement to the implementation record

- The fix recorded in addendum 1 (ii) (fd07605, re-including only the v44 JSON in .dockerignore) did not take effect
  in the image build, and the second submission (image 1793050943290310000) stopped before computing for the same
  reason. Commit 71cd71d puts the whole of docs/verification into the image (about 4 MB), and the third submission
  (image 1793052265534416000) was made. The computations are unchanged. The jobs of the first two submissions are
  kept, not deleted; the three that were waiting to retry were suspended.
- Because the cluster lacked capacity, the jobs of T0− (world_v2_b6t0xb) and the folds of the permutation run
  (world_v2_b6actperm) that had not started were suspended for a time so that the plasmode could run first. The
  suspended jobs had not started, no computation was lost, and they run with the same inputs after resumption.

## Addendum 3 (2026-10-27, when no plasmode result existed): correction of the diagnosis in addenda 1 and 2

- The diagnosis in addendum 1 (ii) and addendum 2 (the re-inclusion in .dockerignore did not take effect) was wrong.
  The Dockerfile copies only docs/protocol and docs/data, so the v44 JSON does not reach the image even when
  .dockerignore puts it into the build context. The third submission (71cd71d, image 1793052265534416000) also stopped
  before computing for this reason (its ten jobs are suspended and kept).
- Commit b787ab1 adds a COPY of the v44 JSON to the Dockerfile. The fourth submission was made after a small pod run
  with the new image (1793059748252581000) confirmed that the file is in it (same arguments; only the memory request
  changed from 12 Gi to 16 Gi). The computations are unchanged.
