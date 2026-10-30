# B-16 Rewriting after the twelfth round of first-time critiques, and checks of the design's sensitivity, the labels, the baseline and the readings (English translation)

<!-- Translation of docs/protocol/amendments/B-16_blind_review12_sensitivity.md, drafted with a language model and not
     yet checked by the authors. The Japanese file is the record and governs where the two differ. -->

Status: **plan (2026-10-29T23:25Z, commit 57fd30e, fixed before the quantities of parts B and C below were computed)**.
The confirmatory verdict (equivalent-null) is not changed. Everything here is post hoc (exploratory).
Protocol: [`preregistration_v2.md`](../preregistration_v2.md) (FROZEN 2026-10-25T20:58:37Z, commit e625924).

Values seen when this was written: the confirmatory analysis and every result of amendments B-2 to B-15 (including B1 to
B11 of v76), and the values that the twelfth-round reports read from the tables. None of the quantities of parts B and C
below has been computed. The sums of squared errors by location group of the 200 reassignments of B-7, which B9 (a) uses,
already exist, and the summary of their distribution (the 95 % range of R over the reassignments) has been seen in Fig. S1b.

## Background

The twelfth round of first-time critiques (NHESS version, manuscript 30 pages and Supplement 38 pages, development
repository version b60b1d2, public repository 33f3c59; reports in [reviews/](../../research/reviews/),
`2026-10-30_paper2_blind12_*.md`) returned: the editor, "return before review"; the three referees, major revision. The
common points:

1. **Title and abstract**: "did not improve" and "declared equivalence" cannot be the headline given the design's
   sensitivity (effects near settlements cannot reach the margin, the reading's recall is 18 %, the positive control gets
   the same verdict, an injected odds ratio of 1.5 is called equivalent in 40 of 40 replicates) and the few clusters. Write
   "no detectable improvement; improvements of the margin's size and effects near settlements cannot be excluded". The
   title names only the water-related flag, but the test used 14 features (every class and the character model). The
   abstract's exploratory comparison (names below earlier records and the JRC maps) does not survive the corrections for
   multiplicity and few clusters, and the abstract cites only exploratory results that disfavour names. Spell out EMS (all
   four).
2. **Summary of sensitivity**: put the ceiling, the misclassification, the positive control, the detectability of injected
   effects (Table S22) and R on the ground near settlements (0.060) beside the verdict (referee 1 M1, referee 2 M7, referee 3
   M11, editor M4).
3. **Estimand**: the protocol names the pooled ΔBS(M1c → M3n) as the primary estimand and uses R to express the margin. The
   paper's "the primary estimand is R" differs from the protocol. The margin is derived from the sample (5 % of the point
   estimate of D), and D depends on an extrapolating M0 (referee 2 M1, M2).
4. **Few clusters and the choice of post hoc construction**: choose one construction that respects the dependence, for a
   reason independent of the results, and treat the others as sensitivity analyses. t(ν) versus t(ν − 1); ν of Δ (1.7)
   versus ν of L+ (2.9). Randomization inference by reassigning names (stratified by the relief of settlement sites);
   refitted estimates without one and without the top five groups (referee 2 M5).
5. **The reference model**: the composition term helps only through the leak of the primary split and worsens the
   predictions on every other split. Make the comparison of names centred within AOIs the main within-area analysis
   (referee 2 M6, editor M6).
6. **Placebo**: a necessary condition of the verdict rests on one reassignment (referee 2 M3).
7. **Inference near settlements**: permuting labels is the wrong null. With names reassigned, the odds ratios centre below
   one (GFM 0.86, uncapped EMS 0.90), contradicting "the estimator itself does not lean below one". The assumptions needed
   for a "bound"; the dependence on the sensor (optical over radar ratio 0.50) (referee 2 M10, referee 1 M3, editor M5).
8. **Labels**: detail maps (71 AOIs), the delay of the images, the threshold of binarisation and a fractional target, a
   uniform permanent-water mask, strata of sensor and schema generation, a table of agreement between the GFM and EMS labels
   (referee 1 M2, M4).
9. **Terrain and hydrology of the baseline**: a surface model, HAND above small channels, no upstream area, no protection
   (referee 1 M5).
10. **Readings**: the toponymists' check (authors), the language of each name and unread names coded as zero, homographs,
    romanisation, knowledge of places in the language-model reading, agreement between the coders, calling the check "an
    agreement between language models and Wikipedia", the edited lexicon as a variant, element groups as hypothesis
    generation, a hydronym reading (names that contain the name of a river) (referee 3, editor M3).
11. **Errors in numbers and references**: many (every minor item of the four reports).
12. **What cannot be verified**: the policy on generative AI, the translations, evidence of the freeze time, the data
    deposit (authors' decisions; part D).

## A. Rewriting (no computation; changes to the paper and Supplement)

1. **Title**: "Settlement names assigned to the nearest grid cells did not detectably improve within-area predictions of
   mapped river flooding" (17 words). It refers to every feature of the names (every class of the lexicon and the
   character model) and says that no detectable improvement was found.
2. **Abstract** (at most 200 words, counted after typesetting): EMS spelled out, no coined phrases. In order: what was
   tested (the lexicon was compiled by a language model and has not been checked by toponymists, stated at that point); the
   protocol's verdict (names slightly worsened the Brier score; the rule declares equivalence in a third of samples at the
   margin); the upper limit of the designated post hoc construction and the few clusters (B1); the design's sensitivity
   (0.7 % of cells near water-related settlements, the ceiling); the dilution (records kept −2 % to 22 % of their value); the
   conclusion (about this reading and assignment; what names record remains open). The exploratory "names below records"
   is removed from the abstract.
3. **Estimand and margin**: describe the primary estimand, as in the protocol, as the pooled Δ = BS(M1c) − BS(M3n), with R
   the ratio that expresses the margin as a share of the baseline gain. State that the margin is derived from the sample
   (5 % of the point estimate) and give its odds-ratio translation in both directions (effects that add and remove flooding)
   and over the 95 % interval of D (1.83 to 2.89, Table S21). Call "about the improvement that recalibration alone brings"
   a post hoc description. State that the analogy with non-inferiority margins set at a lower limit differs from trials whose
   margin comes from external data.
4. **A paragraph on sensitivity** (right after the verdict in Sect. 3.2): (i) 0.7 % of cells within 500 m, ceiling 0.027;
   (ii) the lexicon's recall of 18 % and 6 of the 168 settlements named after a river, stream or lake; (iii) the positive
   control (site HAND added to M1) is also equivalent-null; (iv) with every model refitted, an injected odds ratio of 1.5 is
   called equivalent in 40 of 40 replicates, and 2.0 is positive in 40 but equivalent in 35 (Table S22); (v) R on the ground
   within 500 m of settlements, 0.060 (90 % −0.036 to 0.162). State that the verdict can come out the same whether or not
   names carry moderate information.
5. **Post hoc upper limits**: the construction of B1 (super-clusters of the primary split, top level only, second stage
   refitted) is the designated post hoc result, and the other constructions of Table 3 are sensitivity analyses (the reason
   for the choice is the dependence structure, given by referee 2 before the result). The few-cluster reference follows the
   rule fixed in B3. Remove "bounds the improvement at a few per cent" and write "once the concentration of the variance is
   acknowledged, no upper limit below the margin is robust" (following the reporting rules of part B). Bring the path of the
   analysis (the ratio test adopted after the amendment-11 interval was not inside the margin; the one-sided limit and the
   range chosen after every interval had been seen) into the main text.
6. **The comparison with names centred within AOIs** (M1s → M3d) becomes the main post hoc within-area analysis in
   Sect. 3.3, with M1c → M3n the protocol's record. State in Sect. 2.3 that the composition term helps only through the leak
   of the primary split.
7. **Near settlements**: permuting labels checks only the estimator; inference by reassigning names (values of Table S25:
   two-sided p = 0.03 for the planned GFM estimate, reassignment median 0.86, uncapped EMS reassignment median 0.90) is the
   main inference; correct "the estimator itself does not lean below one"; the sensor dependence (Table S25) and the uncapped
   GFM estimate (0.82) in the main text; the upper limit without two activations (1.41) in the "bound" sentence; for many
   water-related elements the expected direction at the site is negative (referee 3 M13 (d)).
8. **Readings** (referee 3): the sentence on qualifiers in Sect. 3.5 (the relation to floodplains is not tested); the
   sentence on the character model (drop "not the language or region"); state the failed predictions (nested H1, the H3
   ladder); cite H5 as Table S9; bring into the main text the assignment of languages by country, unread names coded as
   zero, the whole-word blocks, the romanisation mismatches, homographs, transferred names and hydronyms; call the check "a
   check of agreement between language models and Wikipedia"; the edited lexicon is a variant; element groups are
   hypothesis generation; local interpretations of names and names of events (Isoda et al. 2019); the language model's
   prompt forbids knowledge of places, with no guarantee that this holds.
9. **Methods details** (referee 1 minor 9–12, 17, 18, 21; referee 2 minor 12, 13; referee 3 minor 9, 10): the penalty of
   the second stage (C = 1) and the result with the penalty chosen by inner cross-validation (B-12 B11, R 0.004); early
   stopping on a random tenth of the training cells; Not Analysed polygons removed from the observed area; traces counted as
   flooded; the exact inland rule; the cap shared in proportion to each AOI's cells; the resolution of HydroRIVERS; the river
   network, absence of defences and variables of the JRC maps; the GFM version; SAR shadow and layover and emergent
   vegetation; the source of the native spelling given to the language model.
10. **Supplement**: remove the author names (NHESS); no findings in the Supplement that are not in the paper (the sensor
    dependence of Table S25 and the island group of Table S27 either mentioned in the paper or moved to the archive); capped
    versus uncapped in Tables 5 and S41; one definition of "R over a constant forecast" in Tables S9 and S10; 2.34 versus
    2.33 in Table S19; the interval named in the caption of Fig. S7; the test behind the Holm p of Table S26; the ordering of
    forms in Tables 5 and S39; the TOST range of Text S1 including the pigeonhole and t versions; Supplement references
    (Carter et al. 2017; Webb 2023).
11. **References**: FABDEM (Hawker et al. 2022), Wikidata and Wikipedia, Tent and Blair (2011), Coates (2006), Gelling
    (1984), Gotway and Young (2002) (change of support), Imbens and Kolesár (2016), MacKinnon and Webb (2017), Davison and
    Hinkley (1997), Gneiting and Raftery (2007), each cited only after its bibliographic details have been checked.
12. **Literature search**: search J-STAGE, CiNii and the onomastic journals (Names, Onoma, Nomina, Beiträge zur
    Namenforschung and others) on the web, and record the queries and dates in Text S1; state that bibliographic databases
    (Scopus, Web of Science) were not available, and weaken the claim in the Introduction.
13. **Numbers and references**: every minor item of the twelfth round (four reports). In particular: write "within 500 m"
    for "effects confined near settlements", "exceeded that of reassigned names", the attribution of ν, 0.036 → 0.005 (done),
    a table for the odds-ratio range, the definitions of the country counts (done), Norway's 68 % (done), the band indicators
    of Eq. (2) (done), the EmsSameOR clash of Table S33 (done), the R of the 2022+ rows of Table 4 is unstable (give Brier
    units too), and "not meaningful" in prose for "not-meaning".

## B. Local computations (script v77; seed 20261105)

1. **The designated post hoc construction** (proposed by referee 2, M5): save the second-stage design of the primary split
   again, with the AOI, country and activation identifiers (`world_eval.py --split primary --arms M0,M1 --stacks
   M1c,M3n,M1s,M3d --dump-stage2`, local, frozen reading). Taking as units the connected components of the primary split's
   location groups linked by a shared activation (super-clusters), draw them at the top level only, 500 times (seed 20261105
   + b), and refit the second stage with weights in each draw (character model fixed). Comparisons: M1c → M3n (the
   protocol's record) and M1s → M3d (centred). Point estimate, 90 % Fieller interval, one-sided 95 % upper limit with z and
   with the few-cluster reference of B3. This construction is the post hoc result in the paper.
2. **Refits without one and without the top five groups**: in the saved design of B1, remove each of the five location groups
   with the largest contributions to the variance of Δ, and all five, from training and evaluation, refit the second stage,
   and report R and Δ.
3. **Rule for the few-cluster reference** (one rule for every construction): compute the effective degrees of freedom from the
   unit contributions to L+(θ) = Δ − θD at θ = the reported upper limit θ_U (iterating θ_U and ν to convergence). The main
   reference is t(ν − 1) (no upper limit when ν − 1 < 1), with z and t(ν) beside it. Also give ν of Δ alone (G*). Recompute
   for every construction with saved resamples (Table 3, B1), and correct the error of B-15 B3, which used the ν of M1c → M3n
   for the centred comparison. Reporting rule: if the main reference's upper limit of any construction exceeds 0.05 or is
   unbounded, the paper says "once the concentration of the variance is acknowledged, no upper limit below the margin is
   robust".
4. **Injected odds ratios of the upper limits**: translate every upper limit of Table 3 and B1 (z and t) into odds ratios with
   the mapping of Table S10 (four no-effect references, nominal scale); "unbounded" where unbounded.
5. **Inference near settlements**: (a) for the planned GFM estimate and the EMS estimate on the same cells, the p value and
   interval of a restricted (null-imposed) wild cluster bootstrap (Webb six-point weights, 9,999 draws), beside the
   unrestricted version of B-15 B10; (b) the same inference for the direct estimate on all ground (1.02 in Table S33);
   (c) the upper 95 % limit of the uncapped EMS estimate without the two activations (EMSR052, EMSR130) in the "bound"
   sentence. Inference by reassigning names (Table S25) uses existing values (no new computation).
6. **Checks of the labels** (models fixed, primary split, corrected code): (a) strata of AOI attributes: sensor of the image
   (radar, optical, unknown), schema generation (0–3), presence of a hydrography layer, AOIs without type and footprint, delay
   of the image (from the event time of the activation to the first image, terciles); in each stratum R (with the stratum's
   yardstick), Δ, and the Eq. (2) estimate on the capped EMS table; (b) scoring with other targets: the Brier score against
   the riverine flooded fraction (f_riverine; fractional target), and labels re-binarised at thresholds 0.10 and 0.25 (same
   evaluable cells); (c) scoring without cells whose mean GSW water occurrence (1984–2021) is at least 50 %; (d) R by the
   distance from the cell to its nearest settlement (within 500 m, 0.5–1 km, beyond 1 km), with the band's yardstick and the
   primary yardstick.
7. **Agreement between the GFM and EMS labels**: on cells both can assess, the 2 × 2 table, agreement, κ, and with EMS as the
   reference the GFM hit rate, miss rate, false-alarm ratio and CSI, by distance band and by the sensor of the EMS map.
8. **Bias analysis of the odds ratios near settlements**: in the 500 m band, give the ground nearest water-related names and
   the other ground each a sensitivity (0.6, 0.8, 1.0) and a specificity (0.99, 0.995, 1.0), correct the 2 × 2 counts, and
   report the range of corrected odds ratios for the planned GFM estimate and the uncapped EMS estimate (the AOI adjustment
   approximated by correcting the counts within the band).
9. **Checks of the readings**: (a) from the sums of squared errors by location group of the 200 reassignments of B-7 (with
   the composition rebuilt), the 95 % interval of R − Rp and the verdict of the nested comparison with each reassignment as
   the placebo, and the share of reassignments that change the verdict; (b) the contributions to Δ by language block and by
   country, and a table of the top ten location groups (contributions to Δ and D, country, cells, flooded share, sensor,
   generation, hydrography layer, most frequent forms); (c) subsets scored with the models fixed: cells whose nearest
   settlement is not on whole arc-minutes; PPL settlements only; cells whose nearest settlement's GeoNames modification date
   (allCountries.zip) is before the event; cells in countries whose first language a block reads and outside regions of
   another language; (d) κ between the two coders by class, and the flow of the check's sample (drawn → with a section →
   with a derivation → coded); (e) a hydronym reading: an indicator that the nearest settlement's name contains the name of a
   river or lake within 5 km (the definition of B-10), added to M1c in the second stage (B1's saved design, refitted), with R,
   Δ and a centred version.
10. **Checks**: the definition of "R over a constant forecast" in Tables S9 and S10; the difference between Tables 5 and S41
    (capped or not); the ordering of forms in Tables 5 and S39; what the Holm p of Table S26 refers to; the list of the 11
    Belgian -au names of Table S28 (iv); how the counts of the check (353, 360, 245) correspond.

## C. Cluster computations (utens; image tags recorded)

1. **A strong baseline** (referee 1 M5): M1+ adds to the base models FABDEM elevation statistics (replacing GLO-30), the
   height above the nearest large river (HydroRIVERS reaches with an upstream area of at least 1,000 km²; with FABDEM, the
   cell's mean elevation minus the elevation at the nearest point of that river), the upstream area of the nearest reach
   (log), the GSW water occurrence, seasonality and maximum extent (1984–2021; including the event year for events up to
   2021), the JRC hazard fractions and depth, and the WorldCover 2021 built-up and tree shares (a map after many events).
   Arms M0, M1+; stacks M1+c, M3+n, M3+np, M1+s, M3+d. Primary split and the split that keeps activations whole. Δ, R (with
   the M1+ yardstick), the centred comparison, the placebo.
2. **A names model that can act near settlements** (referee 1 M1 (b)): stacks M3nb, with the name features multiplied by
   indicators of the cell's distance band to its nearest settlement (within 500 m, 0.5–1 km, beyond 1 km), and the centred
   M3db. Primary split. R on all cells and on the ground within 500 m.
3. **Stratified randomization** (referee 2 M5): within each AOI, reassign the names within terciles of the settlement sites'
   HAND (the mean within 177 m of B-14 C2), rebuild the features and the composition, and refit the second stage (200 times).
   The distribution of R of the nested and centred comparisons, and the one-sided p (share of reassignments with R at least
   the real R).
4. **Redrawing the placebo** (referee 2 M3): redraw, 20 times, the within-AOI reassignment of the protocol for the flexible
   refit on the language-model split (M1 → M3 with placebo M3p), the total increment (M1s → M3s with M3ps) and the nested
   comparison (M3np), and count how often each verdict changes.
5. **Detail maps** (referee 1 M2 (i)): retrieve the detail maps of the 71 analysed AOIs that have them (100 products) and grid
   them on the AOI's 250 m grid by the same rules (locally is fine). Within the detail maps' extent, the hit rate and
   false-alarm ratio of the overview labels against the detail labels by distance band and by whether the nearest name is
   water-related; the names increment (models fixed) and the Eq. (2) estimate on those cells with the detail labels.

## D. Authors' decisions (not done in this amendment)

1. **The policy on generative AI** (editor M1, a precondition for review): NHESS does not allow generative AI for text or
   interpretations. The authors either write the text as their own or choose a journal that accepts declared drafting.
2. Checking the English translations of the protocol and every amendment.
3. Verifiable evidence of the freeze (because the development repository is not to be shared, the authors decide on a
   third-party time stamp or an archive of the freeze commit and its ancestors only).
4. Deposit of the code and derived data (DOI) and access for referees.
5. The toponymists' check (the 300 names of v69, including the blocks and countries that carry the results) and a reanalysis
   with the checked readings.
6. Deposit of the conversation and prompt that produced the lexicon.
7. Publication of machine readings of indigenous names (CARE).
8. Testing whether the language-model reading uses knowledge of places (the same names with other countries, invented
   names); needs GPU computing.
9. Article type and journal.

## Reporting rules

- Everything is post hoc (exploratory), and the confirmatory verdict is not changed. The protocol's verdict is reported as
  the record, with its operating characteristics.
- The post hoc result in the paper is the construction of B1 (designated), and the other constructions are sensitivity
  analyses. The few-cluster reference follows the rule of B3.
- Every quantity of parts B and C is reported whatever its direction. The abstract carries no exploratory result, or, if it
  does, carries results in both directions under the same standard.
- No rule, unit, number of draws or stratification is changed after it has been computed. Implementation errors are corrected
  in addenda.

## Addendum 1 (record of the implementation of parts B and C; 2026-10-29T23:40Z, commit f709cf1, before their quantities were computed)

1. **C2 overlaps B-10**: the model of names × distance bands was already run in C1 of B-10 (run world_v2_b10stk, primary
   split; M1cd = M1c + band indicators, M3nd = M1cd + the 14 name features × band indicators; R −0.0075, 90 % −0.017 to
   0.0013, equivalent-null; v67). The uncentred version of C2 is the same model, so it is not rerun on the cluster. Instead,
   (a) from the saved predictions of world_v2_b10stk, Δ = BS(M1cd) − BS(M3nd) and R (with the band's yardstick
   BS(M0) − BS(M1) and with the overall yardstick) are computed in each distance band (within 500 m, 0.5–1 km, beyond 1 km)
   with 4,000 draws of location groups (models fixed, seed 20261105). The distance is taken from P_d1 of features.parquet
   after joining the cells to cells.parquet by AOI and coordinates rounded to 10⁻⁶ degrees (the computation stops unless every
   cell is joined). (b) The centred version (M1sd = M1s + band indicators; M3db = M1sd + name features centred within AOIs ×
   band indicators) and M1cd → M3nd are fitted locally from the dump of B1 and included in the super-cluster draws of B1.
   R on all cells and on the ground within 500 m (band and overall yardsticks) is reported.
2. **Dump for B1**: `world_eval.py --tag world_v2_b16_stage2_local --split primary --dump-stage2` (local; the confirmatory
   seeds, table seed and caps). P_d1, P_gid1, ix and iy of the training and evaluation cells are added to the dump (extending
   `dump_stage2`). The frozen reading is reproduced by setting the lexicon features of Norwegian names to 0 and rebuilding the
   regional composition and the centred features (B9 of B-15, `_designs(frozen=True)` of v76). The corrected lexicon's
   reading is reported as point estimates only. As a check of the implementation, the differences between BS(M1c) and
   BS(M3n) of the refit with the frozen reading and those of the confirmatory run (world_v2) are reported.
3. **Super-clusters**: the location groups of the primary split (5 km buffer) of the analysed AOIs are joined when they hold
   AOIs of the same activation (`merge_by_activation`). Each cell belongs to the super-cluster of its AOI. Their number and
   the largest super-cluster's share of cells are reported. In each of 500 draws (seed 20261105 + b, multinomial weights) the
   training and evaluation cells are weighted and the second stage is refitted with the weights (`stage2_weighted`, started
   from the unweighted fit; the predictions of the base models M0 and M1 are fixed). Arms: M1c, M3n, M1s, M3d (main); M1cd,
   M3nd, M1sd, M3db (C2); M1cH, M3nH, M1sH, M3dH (B9 (e)). Points are the unweighted refits. Intervals: the 90 % Fieller
   interval from the covariance of (Δ, D) over the draws (its upper end is the one-sided 95 % upper limit), and the percentile
   and normal-approximation limits of L± (`_levels_detail` of v73).
4. **B2**: the five location groups with the largest contributions to the variance of Δ, γ_g = ((e_g − Δ̂ n_g)/N)² (e_g the
   sum over the group's cells of the cell differences for M1c → M3n, from the unweighted refit with the frozen reading of
   B1), are removed one at a time and all five together from training and evaluation, and the second stage is refitted
   (unweighted). Δ, D (on the remaining evaluation cells) and R for M1c → M3n and M1s → M3d. The base predictions are fixed
   (the base models trained with the removed groups are not refitted).
5. **Details of the rule of B3**: from the contributions x_g(θ) = Σ_{i∈g}(e_i − θ d_i) of the construction's resampling
   units, ν(θ) = (Σγ_g)²/Σγ_g² with γ_g = ((x_g − x̂ n_g)/N)² (`_nu` of v76). θ_0 is the Fieller upper limit with z, and
   θ_{k+1} the Fieller upper limit with the critical value t_{0.95}(ν(θ_k) − 1), iterated until the change is below 10⁻⁶ or
   for at most 100 iterations. If ν − 1 < 1, or the Fieller interval is unbounded, the result is "no upper limit". z and
   t(ν) are reported with the same ν. G* is ν at θ = 0 (Δ alone). Units: super-clusters for B1; location groups for (a) and
   (b) of Table 3; countries for (c) and (d); the smaller of location groups and activations for the pigeonhole (as in v76);
   super-clusters for B3 of B-15. For the centred comparison ν is computed from its own cell differences.
6. **B4**: the mapping to odds ratios is `_or_conversions` of v73 (the curves of B2 of B-13 and the four no-effect
   references). "No upper limit" stays "no upper limit".
7. **Inference of B5**: a restricted wild cluster score bootstrap (in the form of Kline and Santos 2012). With Eq. (2)
   (`design` of v39, banded, the logit of M1 as a covariate, AOI fixed effects), the model is fitted under the null
   δ_≤500 = δ₀ (that column as an offset of δ₀); the efficient score of column j, adjusted for the other coefficients with the
   information matrix after the AOI effects are profiled out, is summed by location group (s̃_g), and
   W = (Σs̃_g)²/Σs̃_g², W*_b = (Σv_gb s̃_g)²/Σ(v_gb s̃_g)² (Webb six-point weights, 9,999 draws, seed 20261105, the same v for
   every δ₀), p = (1 + #{W* ≥ W})/(1 + 9,999). The 95 % interval is the set of δ₀ on a grid of log odds ratios (−1.6 to 1.6,
   step 0.01) with p > 0.05 (if not contiguous, its outer range, with a note). Targets: (a) the planned GFM estimate and the
   EMS estimate on the same cells (v41 table, the same rows as B10 of v76); (b) the direct estimate on all ground (1.02 in
   Table S33: the unbanded water-related indicator on the capped EMS table, the same rows as ORWetAllAdj). (c) cites the
   existing value (1.41).
8. **B6**: predictions of the corrected run world_v2_lexfix (primary split, 10 folds, models fixed). The yardstick is
   BS(M0) − BS(M1) for each stratum and target. (a) Sensor from gfm_summary (s1 and sar as radar, optical as optical, none as
   unknown); generation from meta.json; the hydrography layer from layers.hydro of meta; AOIs without type and extent by the
   rule of B-6 6(d); the imaging delay in days from the event time (event_time of the inventory, else activation_time) to the
   first post-event image (Post-event dates of the v38 source table, else image_times of the inventory), split into terciles
   of the AOIs with a known delay (AOIs as units, unweighted), with AOIs of unknown delay as a separate stratum. In each
   stratum, R, Δ and the capped EMS estimate of Eq. (2) (cluster-robust 95 % interval, `fit` of v39). (b) The fractional
   target is the mean of (p − f_riverine)² on the same cells; the labels with thresholds t = 0.10 and 0.25 are
   y_t = 1[f_riverine ≥ t] on the same cells. (c) Cells with occ (%) of cell_gsw at least 50 are removed (cells without GSW
   are kept). (d) R in each band (band and overall yardsticks). Intervals from 4,000 draws of location groups (seed 20261105).
9. **B7**: cells with both GFM and EMS labels in the v41 table. With EMS as the reference: hit rate = P(GFM=1 | EMS=1), miss
   rate = 1 − hit rate, false-alarm ratio = P(EMS=0 | GFM=1), CSI = a/(a + b + c), and Cohen's κ; overall, by band and by
   EMS sensor.
10. **B8**: within 500 m, sensitivities Se ∈ {0.6, 0.8, 1.0} and specificities Sp ∈ {0.99, 0.995, 1.0} are given to the ground
    of water-related names and to other ground separately (9 × 9 = 81 combinations), and the flooded counts are corrected as
    a = (a_obs − (1 − Sp)N)/(Se + Sp − 1). Combinations giving a corrected count of 0 or less, or N or more, are dropped as
    impossible. Corrected odds ratio = estimated (AOI-adjusted) odds ratio × corrected crude odds ratio / observed crude odds
    ratio. The range over all combinations and over the 9 with the same Se and Sp on both grounds. Targets: the planned GFM
    estimate (v41 table) and the uncapped EMS estimate (the table of NearUncapOR).
11. **B9**: (a) Rp_j = [SSE_ref,j − SSE_names,j]/(N D̂) (the increment from the rebuilt composition M1 + c(q_j) to
    M1 + c(q_j) + q_j of reassignment j); for the centred comparison SSE(M1s) − SSE_cent,j. The 95 % interval of R − Rp_j from
    4,000 draws of location groups (models fixed, seed 20261105, the same draws for every j); the verdict from `decide` of
    world_aggregate with the interval of R − Rp_j as the placebo interval. The text states that, unlike the protocol's placebo,
    these rebuild the composition. (b) The language block is the one the lexicon assigns to the country of the nearest
    settlement. Contributions are Σe_i/N (they add up to Δ). The top ten location groups are ordered by γ_g as in B2. (c) (i)
    Cells whose nearest settlement does not have both latitude and longitude on whole arc-minutes (×60 within 10⁻⁴ of an
    integer); (ii) nearest settlement with feature code PPL; (iii) nearest settlement whose GeoNames modification date
    (allCountries.zip) is before the AOI's event date; (iv) countries where the block the lexicon assigns reads the first
    language in the Languages column of countryInfo.txt (ISO 639-1, the part before "-"; no restriction by language region,
    since no data on language regions are at hand). Δ and R (subset and overall yardsticks) for each subset, 4,000 draws of
    location groups, scored with the predictions of the unweighted refit with the frozen reading of B1. (d) Cohen's κ by
    class on the part of the audit coded by a second coder, if such codes exist; otherwise, a statement that none exist. The
    flow of the sample from the records of v69. (e) H is `hydronym` of v57_hydronym.parquet (the definition of B12 of B-8: the
    settlement name contains the core of the name of a GeoNames water feature within 5 km), joined by P_gid1.
    M1cH = M1 + the composition of names + c(H), M3nH = M1cH + H; centred, M1sH = M1s → M3dH = M1s + (H − c(H)).
12. **C1**: arm M1p (strong baseline) = T0 + TER_4 (GLO-30, including HAND; not replaced) + the six elevation statistics of
    FABDEM (columns fb_*, features_fabdem.parquet of C7 of B-9) + h_above_100 and h_above_1000 (same file) + log10 of the
    upstream area of the nearest channel (new: UPLAND_SKM of the nearest vertex among all HydroRIVERS channel vertices within
    the AOI ± 0.2°, read as by `river_vertices` of `world_fabdem.py`; features_upland.parquet) + hz_f10, hz_f100, hz_f500 and
    hz_d100 (features_jrc.parquet of C5 of B-10) + o_gswc_occ, o_gswc_seas and o_gswc_max (cell_gsw.parquet) + lc_built and
    lc_tree (features.parquet). Missing values are left as they are (gradient boosting). Stacks M1pc = M1p + regional
    composition, M3pn, M3pnp, M1ps, M3pd, M3pdp. Splits: primary and `--group-by-activation`. R uses the yardstick of M1p,
    BS(M0) − BS(M1p). Verdicts from `_verdicts` of v56 (4,000 draws of location groups, models fixed, the protocol's rules).
13. **C3**: in each AOI, the gids of the P names referenced in it are split into strata by the within-AOI terciles (quantiles
    1/3 and 2/3; ties to the lower stratum) of the settlement-site HAND (computed as o_site_hand of C2 of B-14) and permuted
    within strata. Settlements without a site HAND form a separate stratum. The random numbers are made as in
    `placebo_perms`, with the stratum number added (seeds j = 1 to 200). The features (lexicon and N-learn) and the regional
    composition are rebuilt and the second stage refitted (as in `run_perms_recomp` of B-7).
    p = (1 + #{R_j ≥ R̂})/(1 + 200), and the share #{R_j ≥ R̂}/200.
14. **C4**: the verdicts that depend on the protocol's placebo are two: the flexible comparison on the language-model cut-off
    split (not-meaning) and the total increment of N-read (positive-small, with the placebo condition met) (confirmatory
    values, B1 of v76). The redraws run the language-model cut-off split with N-read and every default arm and stack (M3p,
    M3rp, M3ps, M3rps, M3np and M3rnp change with each reassignment). Seeds via `--placebo-seed` are 20261105 + k
    (k = 1 to 20). Verdicts from the protocol's aggregation (world_aggregate); for each comparison, the number of redraws whose
    verdict differs from the protocol's.
15. **C5**: analysed AOIs for which the old distribution holds a DETAIL product with the same AOI number (MONIT excluded, the
    latest version of each DETAIL number) and whose labels come from a product other than DETAIL. The layers of all DETAIL
    products are treated together as one product and rasterized onto the AOI's grid (meta.json) with the same rules (flood,
    trace, extent, observed and permanent-water fractions). Cells compared are those assessable in both. Hit rate, miss rate,
    false-alarm ratio and CSI of the overview labels against the detail labels, by band and by presence of a water-related
    name. With the predictions of the confirmatory run (models fixed), Δ, D and R on the cells within the detail maps' extent
    scored with the detail labels and with the overview labels, and the capped EMS estimate of Eq. (2).
16. **Scripts**: local computation in `backend/verification/r_okauchi/v77_b16_local.py`, cluster computation in
    `v78_b16_cluster.py`, detail maps in `v79_b16_detail.py`. Output `docs/verification/v77_b16_local.json` and others.

## Addendum 2 (correction of the implementation of C5; 2026-10-30T00:15Z, commit 50c0c75, before its comparisons were computed)

After the rasterization of C5 (`--build` of v79) and before the comparisons (`--score`), we found that the area-of-interest
layer of the DETAIL products of generation 0, and of some of generation 1, outlines the whole AOI, not the frame of the
detail map (AOI01 of EMSR026: the layer covers about 3.8 square degrees, while the flood polygons lie within about
0.15° × 0.11° around a town). Taking that layer as the extent of the detail map would count ground outside the frame as
observed and not flooded on the detail map. AOIs whose DETAIL products' area-of-interest layers together cover at least 0.9
of the AOI's cells (the layer outlines the whole AOI, so the extent of the detailed delineation is unknown) are therefore
removed from the comparisons of C5, and the removed AOIs are listed.

Only the record of the rasterization was seen (for each AOI, the number of products, their generations, the number of cells
inside the area-of-interest layer and the number of cells evaluable on both maps); no quantity of the comparisons had been
computed. The share of cells covered by the layer splits into values of 0.572 or less and 0.965 or more (the latter: 21 AOIs
of generation 0 and 3 of generation 1), and the threshold 0.9 lies between them.

## Addendum 3 (correction of the implementation of C4; 2026-10-30T01:20Z, commit 71573dc, before the redraws with the frozen reading were computed, after the counts of the redraws with the corrected code were seen)

The 20 redraws of C4 (world_v2_b16plac01 to 20, seeds 20261106 to 20261125) ran on the current code, which reads Norwegian
names with the lexicon (the correction of amendment 4). Their held-out predictions for every arm other than the placebo
therefore equal those of the corrected-code run (world_v2_lexfix, protocol seed), not those of the confirmatory run
(world_v2, frozen reading): in all 20 redraws the largest difference over the 37 non-placebo arms and stacks is 0, against
up to 0.67 from the confirmatory run. C4 assumed that only the placebo arms change between redraws, but for the three
comparisons that use the lexicon (nested M1c to M3n, total M1s to M3s, flexible M1 to M3) the difference from the
protocol's verdict mixes the correction of the reading with the redraw. The three comparisons of H4, which use only N-read
(M1rc to M3rn, M1s to M3rs, M1 to M3r), are not affected: their arms and their original placebo arms (M3rnp, M3rps, M3rp)
have the same predictions in the confirmatory and the corrected-code runs (largest difference 0).

Values seen when this addendum was written:

- C4 of v78 (number of redraws whose verdict differs from the protocol's verdict): nested 0/20, total 0/20, flexible
  20/20 (all inconclusive), H4 nested 0/20, H4 total 18/20 (not-meaning 18, positive-small 2), H4 flexible 0/20.
- Verdicts of the corrected-code run on this split (placebo with the protocol's seed): nested equivalent-null, total
  equivalent-null, flexible inconclusive (R 0.110, 95% −0.0003 to 0.250), H4 nested equivalent-null, H4 total
  positive-small, H4 flexible inconclusive.

Corrections:

1. **Redraws with the corrected code** (the existing 20): count the redraws whose verdict differs from the verdict of the
   corrected-code run (world_v2_lexfix, placebo with the protocol's seed). For the three comparisons of H4 these counts
   equal those against the protocol's verdict.
2. **Redraws with the frozen reading** (new runs): add `--frozen-norway` to `world_eval.py`, which replaces the key `NO`
   of the lexicon's country-to-language mapping with `FALSE`, as the YAML 1.1 reading at the freeze did (Norway then has no
   language); nothing else changes. (a) Check: with this flag and the protocol's seed (20261025), run the language-model
   cut-off split with arms M0, M1, M3, M3p and stacks M1c, M3n, M3np, M1s, M3s, M3ps only (world_v2_b16placf00, without
   N-read), and confirm that their held-out predictions equal those of the confirmatory run (world_v2; largest difference
   at most 10⁻⁹; otherwise the results of (b) are not used and this is reported). (b) With the same flag, arms and stacks,
   run seeds 20261105 + k (k = 1 to 20; world_v2_b16placf01 to 20). (a) and (b) may be submitted together, but the verdicts
   of (b) are not computed until (a) has passed. The verdicts of the three comparisons that use the lexicon come from
   `_verdicts` of v56 (the protocol's aggregation), and we count the redraws whose verdict differs from the protocol's.
3. **Reporting**: the paper's result of C4 is 2. (redraws with the frozen reading) for the three comparisons that use the
   lexicon, and the existing redraws of 1. for the three comparisons of H4 (same reference as the protocol's verdict). The
   counts of 1. with the corrected code, and those of C4 in v78 (which mix the correction of the reading), are listed in a
   table of the Supplement, the latter marked as mixed.

## Addendum 4 (correction of the implementation of C1; 2026-10-30T03:13Z, commit 569f5e5, before the centred comparison was computed, after the nested comparison and the redraws with the frozen reading were seen)

The runs of C1 (world_v2_b16strong and world_v2_b16strongact) did not run the stacks M3pd and M3pdp (names centred within
AOIs over the strong baseline model, and their placebo), because `b6_stack_available`, which decides whether a stack is
run, treated the centred features `d_*`, which are built inside each fold, as missing columns of the table (the run record
says "stacks not run (features missing): ['M3pd', 'M3pdp']"). The nested comparison (M1pc to M3pn, placebo M3pnp) and M1ps
were run.

Values seen when this addendum was written:

- The nested comparison of C1 (baseline gain of M1p): primary split R 0.000 (90% −0.0024 to 0.0027), equivalent-null;
  split keeping activations whole R 0.0039 (90% 0.0011 to 0.0069), not-meaning; the Brier score of each arm.
- The redraws of C4 with the frozen reading (addendum 3, item 2): the check run equals the confirmatory run (difference 0
  over the 10 arms and stacks); over the 20 redraws, the nested and the total comparisons were equivalent-null each time,
  and the flexible comparison was positive-meaningful 14 times and not-meaning 6 times (the reassigned names' R from
  −0.189 to 0.166, median 0.026; 0.137 for the protocol's draw).

Corrections:

1. In `b6_stack_available`, a feature that begins with `d_` can be run when its name without `d_` is a column of the table
   or a feature built inside each fold. Nothing else changes (a test is added).
2. On both splits, rerun arms M0, M1 and M1p and stacks M1ps, M3pd and M3pdp (world_v2_b16strongd and
   world_v2_b16strongactd; seeds, table and folds as in C1). First check that the held-out predictions of M0, M1, M1p and
   M1ps equal those of the C1 runs (largest difference at most 10⁻⁹; otherwise the centred comparison is not used and this
   is reported). If they agree, compute the centred comparison (M1ps to M3pd, placebo M3pdp) with the aggregation of C1
   (`_verdicts` of v56, 4,000 draws of location groups, models fixed, baseline gain of M1p).
3. Reporting follows the rules of amendment 16 (both splits, whatever the direction).
