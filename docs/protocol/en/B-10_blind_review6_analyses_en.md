# B-10 Corrections to the reporting and exploratory analyses in response to the sixth round of first-time critiques (NHESS version, 2026-10-29)

English translation, drafted with a language model; the Japanese file `docs/protocol/amendments/B-10_blind_review6_analyses.md`
is the record and governs where the two differ.

Status: **plan (2026-10-29T00:36Z, commit 5d48bb3, before any quantity in B and C below was computed)**. The confirmatory verdict
is not changed. Everything below is **post hoc** and responds to an editor and three referees who read only the PDFs
of the manuscript and its supplement for the first time (language-model agents, told that they could consult the public
repository; reports in `docs/research/reviews/2026-10-29_paper2_blind6_*.md`, commits d1f426c and 583aa42). The paper
describes these analyses as "fixed in this document after the confirmatory result had been seen and before they were
computed". Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924). Code: with the correction of
amendment B-4 (corrected code). Items C1 to C4 of amendment B-9 (low-δ grid, injections that remove flooding, effects that
differ between settlements, power with refitted simulations) were running on the cluster when this document was written;
this amendment does not change how they are computed or reported.

Values seen when this document was written: all results of the confirmatory analysis, of B-2 to B-8 and of items C5 to C7
of B-9 (as printed in the version of the paper at commit 28cc841), and values that the reports point out or estimate. In
particular: the direct terrain-conditional odds ratio 1.02 (95 % 0.83–1.25) and 1.01 (0.77–1.31) within 500 m; the
residual correlations of B-7 item 11(g) (−0.74 below 50 km, −0.21 at 50–100 km, −0.10 overall over 359 pairs, against
−0.003); the subset of B-7 item 4 (R 0.081, 0.005–0.172, 227 AOIs); the GFM-refitted increment of the names (+5.7×10⁻⁴)
and its Brier scores (M1c 0.1155, M1s 0.0778); the difference between the total increment of the language-model readings
and its placebo (0.007–0.045); real minus placebo in the country split (+0.015); and that the CSI of the HAND-only reference
moves between 0.03 and 0.25 with its definition (referee 1's estimate). None of the quantities in B and C below has been
computed.

## Main points of the critiques (common to all four)

1. The confirmatory test bounds an area-wide effect (on the Voronoi ground of the nearest settlement); an effect at the
   scale of settlements cannot be registered at any strength. The title, abstract and conclusions still lead with the
   confirmatory verdict. The abstract should say that the settlement-scale claim is untested
2. The odds-ratio bounds in the abstract (1.23–1.57) depend on a simple model of the labels (the terrain prediction is the
   true probability, flips are independent between cells, the effect is the same in every location group), and are read
   on a straight line from zero below the first grid point (1.28), which reads the bound too low for a convex curve. Put
   the direct estimate (1.02, 95 % 0.83–1.25) in the abstract
3. Quantities that were planned but not reported: B-7 item 11(g) (correlation of held-out residuals between AOIs of the
   same activation), the subset of B-7 item 4 (the recurrence control in the 227 AOIs with other events), and B-8 item
   B5(f) (the GFM-refitted increment of the names and the Brier scores of M1c and M1s)
4. Errors: the GFM calibration slope 0.16 belongs to M1c (the terrain model's is 0.38); the note to Table S28 says "2.49
   standard errors" (2.80) and "20 assignments" (13); the HAND of the "island" group was computed over all cells, not at
   the settlements, and the group includes German au, aue, auen (water meadows); the achieved bounds 0.011 and 0.017 come
   from different versions of the code (0.013 and 0.017 in the same version); the range of the H5 strata; "one exception";
   values that appear only in the text; values in Table S58 that are not in the paper; and others (the 21 items of
   Sect. 4 of the editor's report)
5. Readings of the names: no toponymist has checked them; folding to ASCII (å, -näs); forms with several origins (-ach
   from -acum, English bank as "slope", doab); mixed classes (springs and wells, the relative modifier "lower"); regional
   and indigenous languages read with the national language (Welsh, Gaelic, Irish); validation is between language models
6. Labels and reference: if misclassification differs between name groups, the odds ratio need not move toward one;
   sensor and product generation; the definition of permanent water; the reference is weak (median CSI 0.13); the
   HAND-only reference is full of ties and is not defined (0.03–0.25 depending on the definition); use HAND along the
   channel network or a modelled flood-hazard map as a reference
7. Statistics: AOIs of the same event fall on both sides of a fold (M1c collapses in the activation-merged split); the
   intervals hold the model fixed; the second-stage penalty does not regularise the AOI-level features; there is no
   interaction between names and distance
8. Length: halve the main text, move most of the post hoc detail to the supplement, add a list of model symbols

## A. Corrections to the reporting (no computation; in the paper and supplement)

1. **Framing**: in the title, abstract, introduction and Sect. 2, before any result, state what was tested (the area-wide
   effect of names assigned to the nearest settlement) and what cannot be detected (an effect at the scale of
   settlements). Remove the verdict label from the abstract. Say in the abstract that the settlement claim is untested,
   and that gazetteer names are not local knowledge (practice, oral tradition, field names). Limit the conclusions to
   "readings compiled by a language model" (including in the title)
2. **Numbers in the abstract**: remove the odds-ratio bounds from the abstract; give the direct estimates (terrain-
   conditional odds ratios over all ground and within 500 m) and say in Brier and CSI units what lying inside the margin
   means. Next to the confirmatory verdict, say that the same verdict appears for records of past flooding and for surface
   water assigned through the settlements
3. **Margin**: say that it has no flood-risk justification (it was set by convention). Describe "about what calibration of
   the terrain model alone gains" as a post hoc observation
4. **Planned but unreported quantities**: report B-7 item 11(g), the subset of B-7 item 4 and B-8 item B5(f) in the text
   or supplement. Put in the main text the total increment of the language-model readings (the only positive verdict, in
   the cut-off split; 95 % interval of the difference from its placebo 0.007–0.045) and real minus placebo in the country
   split (+0.015)
5. **Corrections of errors**: correct all 21 items of Sect. 4 of the editor's report and every factual error among the
   referees' minor points (attribution of the slope 0.16, the note to Table S28, the code version of the bounds, the range
   of H5, how the sensitivity analyses are counted, the code version of the water-side share, the margin as OR 2.33 versus
   2.34, the definitions of near cells and location groups, the count of Swedish -as, the count of etymologies, the terms
   of the Murphy decomposition, cross-references, the number of units in the country split, numbers of repetitions, the
   definition of the HAND-only reference, the contents of Table S58, values that appear only in the text)
6. **"Toward one"**: rewrite to say that misclassification that differs between name groups (missed flooding inside
   buildings concentrates at settlements; water-side settlements differ in height above the river) can bias in either
   direction
7. **Rationale for the composition term**: the placebo hardly changes the set of names in an AOI, so it cannot test whether
   the composition carries meaning. Say in the abstract and methods that the estimand excludes between-area information by
   design. Item 3.1(d) of the editor's report is answered by C3
8. **The "island" group**: rename it "islands and water meadows", report the value at the settlements (cells within 177 m)
   and without au, aue, auen (B1), and correct the claim
9. **Form**: halve the main text (about 6,500 words). Turn Sect. 3.4 into a table; move Fig. 4b and most robustness numbers
   to the supplement; move Appendix A to the supplement, keeping a short list of deviations. Add a list of model symbols
   as an appendix. Label the interval type of each row of Fig. 7. Note f of Table 3 becomes "50 each". Keep "water-side"
   as the protocol's name and define it at first use
10. **Literature**: HAND-based inundation mapping and comparisons of global models with observations (Johnson et al.
    2019, Dottori et al. 2016, Sampson et al. 2015, Trigg et al. 2016, Nardi et al. 2019), the Sendai Framework, Kelman
    et al. 2012, Brázdil et al. 2006, Di Baldassarre et al. 2015, and intervals for a fixed model (Nadeau and Bengio 2003,
    Bates et al. 2024). Check each reference against Crossref before adding it
11. **Public repository**: bring README and PROVENANCE up to B-10 and v66, correct the number of critique rounds, say that
    the gateways package is an unused template

## B. Computations from existing outputs (local, script v66; random seed 20261030)

1. **"Islands and water meadows"**: the within-AOI difference in HAND at the settlements (cells within 177 m of the
   settlement point), and for the group without au, aue, auen (same definition) the differences over all cells and at the
   settlements and the odds ratio within 500 m (same formula as B-5 item 8)
2. **Floodplain reading** (the riverine and wet classes of B-7 item 5): the odds ratio leaving out one country at a time
   (each country with flooded group cells within 500 m), leaving out activation EMSR097, the numbers of group cells and of
   flooded cells within 500 m by form (composition), and the odds ratio re-read without the following forms: English
   bank, banks; German au, aue, auen; forms for embankments and reclaimed land (levee, dijk, dam, digue, dige, hraz,
   argine, argini, dig, nasip, dambis, pylimas, tanggul, penger, polder, bonifica). Same formula as the main version of
   B-8 item 6 (AOI fixed effects + logit of the terrain prediction + class × distance band; location-group robust
   interval and pairs cluster bootstrap with 2,000 resamples)
3. **Power of the odds ratio near settlements**: from the standard error of the log odds ratio in the pairs bootstrap of
   the main estimate (water-side × 500 m), the smallest odds ratio detectable with two-sided 5 % and 80 % power
   (exp(2.80 × SE), in both directions)
4. **Intervals for the main 500 m estimate**: 95 % intervals from the location-group jackknife (CV3, t distribution) and
   from the wild cluster bootstrap (Rademacher, 2,000 draws, score version with the null imposed)
5. **Range of the misclassification correction**: with the formula of B-8 item 6(f), OR = (OR*(1 − p0) − (1 − p1))/(p1 −
   OR* p0), on the 16 points of the grid p1 = 0.6, 0.7, 0.8, 0.9 (share of water-related names among names read as
   water-side) and p0 = 0.05, 0.10, 0.20, 0.30 (share of water-related names among the other names). OR* is the point
   and the pairs 95 % limits of the main estimate (within 500 m). Points where the denominator is zero or negative are
   reported as "not defined"
6. **Validation of the readings by element cluster and by country**: precision and recall of the two coders by
   dictionary element cluster and by country (clusters and countries with at least 10 coded names)
7. **Sensor and product generation**: R of the main comparison (corrected code; the same predictions scored on the subset,
   the yardstick also on the subset) and the odds ratio near settlements (pairs), by the imagery of the EMS map (Sentinel-1,
   other radar, optical, not recorded) and by product generation (0–3)
8. **Main comparison without Africa**: R from the corrected code's predictions without the 40 African AOIs (yardstick on
   the same sample)
9. **Activation-merged split**: with the corrected code's run world_v2_b7actfix, R under the 5 yardsticks × 2 bootstraps,
   and a table of the Brier scores of M0, M1, M1s, M1c and M3n in each split (main, country, time, cut-off,
   activation-merged)
10. **Within-AOI evaluation with GFM**: on the held-out predictions of the run that refits every model to GFM (B-3), the
    within-AOI evaluation of B-9 item B4 (each AOI's flood rate given)
11. **The Sect. 5 estimator on injected labels**: fit the formula of Sect. 5 (AOI fixed effects and the logit of the
    terrain prediction) to labels from the fixed-base injections (forms A and B, δ = ln 1.25, ln 1.5, ln 2, ln 3, 20
    each) and set the mean odds ratios within 500 m and over all ground beside the nominal odds ratio
12. **Aligning the HAND-only reference**: rank the HAND-only map (i) by minimum HAND (ties broken by mean HAND, then in
    random order) and (ii) by mean HAND, and report CSI, hit rate and false-alarm ratio both as the median over AOIs and
    pooled over all cells (the terrain model and the names model in the same two ways). Align Tables S15, S20 and S46
    with these definitions
13. **Between-area relation within language areas**: in AOIs of countries read with the German, English and French
    clusters, regress the AOI flood rate (logit) on the composition of the names (share water-side) and the AOI mean of
    the terrain model's prediction (activation fixed effects, standard errors clustered by activation), and describe the
    composition coefficient
14. **Cells read in another language**: (i) at country level, countries where the first language in the languages column
    of GeoNames countryInfo is not among the dictionary languages applied to that country (the group codes bantu and sh
    count their member languages as the same); (ii) within countries, first-level administrative units of Natural Earth
    1:10 m (including the region and geonunit fields): Wales, Highland and Eilean Siar in Scotland, Ireland and Northern
    Ireland, Catalonia, Valencia, the Balearic Islands, the Basque Country, Navarre and Galicia, Brittany, Alsace,
    Friesland, Finnmark and Troms in Norway, Norrbotten in Sweden and Lapland in Finland. For each, the number of analysed
    cells, how many were read as water-side, and how many flooded within 500 m. Record the list of matched units
15. **Distribution of the language-model (Gemma) reading classes**: number and share of names by N-read class (nearest
    names of analysed cells; unweighted and weighted by cells)

## C. Cluster computations (utens, corrected code, main split)

1. **Interactions of names with distance bands, and region-specific coefficients**: (a) reference M1cd = M1c + indicators
   of distance band (within 500 m, 500 m–1 km, beyond 1 km), names model M3nd = M1cd + name features × band indicators (no
   main effects of the name features). (b) Reference M1cr = M1c + continent indicators, names model M3nr = M1cr + name
   features × continent indicators. Both with R = [BS(reference) − BS(names)] / [BS(M0) − BS(M1)] (same yardstick,
   intervals and verdict rules as the main comparison; they do not replace the confirmatory verdict). The continent is
   taken from the country of the cell's nearest settlement
2. **Nested cross-validation of the second-stage penalty**: choose the L2 strength C of the second stage of M1c and M3n
   by nested 5-fold cross-validation (Brier) over training location groups, from 13 points between 10⁻⁴ and 10² (M1ct,
   M3nt). Report R and the chosen C
3. **Reference with country indicators**: a reference with indicators of country (that of the cell's nearest settlement)
   in the second stage (M1k), the same plus the composition (M1kc), and the same plus names (M3kn). Report R of the names
   (M1kc → M3kn) and the increment of the composition after the country indicators (M1k → M1kc, divided by the
   yardstick). Indicators of countries absent from the training side are all zero and have no effect on the evaluation
   side (recorded)
4. **Spatially clustered injection**: inject on the same target as form A (all ground nearest a water-side name), with
   flips correlated within the ground of a settlement (from a standard normal z_s shared by the settlement and e_i per
   cell, u_i = Φ(√ρ z_s + √(1−ρ) e_i), ρ = 0.5; the marginal flip probability is the same as in form A). δ = ln 1.1,
   ln 1.25, ln 1.5, ln 2, 20 each, fixed base. Report the gain curve and the bound (B-8 criterion) beside form A
5. **A modelled flood-hazard map as reference**: from the inundation depths of the JRC global river flood hazard maps
   v2.1.2 (CEMS-GloFAS, 3 arcsec, ODbL) for return periods of 10, 100 and 500 years, compute for each cell the inundated
   fraction (three return periods) and the mean depth at 100 years, and add them to the terrain model (M1h; the location
   model M0 is unchanged). R = [BS(M1hc) − BS(M3hn)] / [BS(M0) − BS(M1h)]. **Condition to run**: the map is available for
   at least 95 % of the analysed cells (sea and tiles without a map set to 0, and recorded). **Caveat**: the return-period
   flows are estimated from a reanalysis that includes the event years, so the map is not a quantity available
   beforehand; it is reported as a control of how strong the reference is
6. **Reporting rules**: C1–C3 and C5 are given verdicts with the same decision table as the main comparison, but do not
   replace the confirmatory verdict. For C4 the bound is read as in B-8 item B1, written ">" outside the grid

## D. Not done in this amendment (authors' decisions, or a separate study)

1. Checking by a toponymist (the sample of 100 names and the high-volume element clusters), and rebuilding the checking
   guide on authoritative place-name dictionaries (the guide is rewritten under A; the authors commission the check)
2. Third-party preservation of the development history (including the freeze), disclosure of the freeze record to the
   editor, the design document and pilot outputs, a summary of the Japanese tests, checking the English translation of the
   protocol, and a statement of what each author verified (authors' decisions)
3. Deposit of derived data and code with a DOI (done by the authors)
4. A settlement-scale test (outcomes per settlement over many events, building footprints, an element set fixed by a
   toponymist, third-party registration), labels from detailed maps, stratification by FLOPROS, HAND along the channel
   network (MERIT Hydro requires registration): a separate study
5. Publication of the instructions and records from building the dictionary (the authors check whether records exist)

## Implementation and resources

- B is computed by `v66_b10_local.py`. For C, `world_eval.py` gains the stacks M1cd, M3nd, M1cr, M3nr, M1ct, M3nt, M1k,
  M1kc, M3kn, M1h, M1hc and M3hn, `--plasmode-forms AC` (with `--plasmode-rho`) and `--hazard jrc`; `world_jrc.py` fetches
  the hazard tiles and builds per-AOI features. Aggregation by `v67_b10_aggregate.py`
- C1–C3 and C5 are 10-fold evaluations; C4 is split into Jobs by ranges of fold × repetition (as in B-9 item C3)
- Random seed 20261030

## Addendum 1 (2026-10-29T01:06Z, commit 38e8dc9, before any result of C and after some results of B): implementation record and one addition

Seen: B1, B6, B12, B13 and B15 (first run of v66; for example, B1 gives a difference in HAND at the sites of "islands and
water meadows" settlements of −3.2 m, and −5.8 m without au, aue, auen). No quantity of C, and none of B2–B5, B7–B11 or
B14, has been seen. This addendum records the implementation and adds an arithmetic translation that uses no data (item 3
below); it changes no rule for computing or reporting any quantity.

1. **Cluster (C) implementation** (code eacf687, image 1793234945410716000)
   - C1: the name features are the 15 columns of the main comparison (the 6 dictionary classes × nearest name and count
     over three names = 12 columns, water-side, any match, N-learn). M1cd has 3 band indicators; M3nd has 15 × 3
     interaction columns. Continents are the GeoNames countryInfo continent codes (AF, AS, EU, NA, OC, SA). Cells
     without a country of the nearest settlement have all indicators zero
   - C2: the inner folds are the 5 location-group folds already used to cross-fit the base model. The second-stage inputs
     are standardised over the whole training side before the inner split. C is chosen per outer fold, with the Brier
     score pooled over all inner held-out cells. The chosen C and the grid of Brier scores are kept in the fold record
     (b10_tuned)
   - C3: the country levels are the countries of the nearest settlements of all cells in the table (no labels used). The
     fold record keeps the countries present on the training side
   - C4: the settlement normal z_s depends only on the repetition ("B-10/AC/settlement/{rep}", the same for every δ). The
     cell uniform e_i uses the form-AC, δ and repetition salt (built as in B-6). Cells without a settlement are
     independent. Repetitions 0–9 and 10–19 run as two Jobs
   - C5: from the JRC server (https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard/, v2.1.2), the
     RP10, RP100 and RP500 depth and permanent-water files of the 61 10° tiles covering the analysed AOIs (244 files,
     3.58 GB) were fetched on 2026-10-29 (no fetch errors). Depth nodata (−9999) is dry land or sea (not distinguishable).
     Permanent-water pixels (value 1 in the permanent-water tile) are removed from numerator and denominator. Rasters are
     resampled by nearest neighbour to the 25 m subgrid and aggregated per 250 m cell. A cell's cover is the share of its
     pixels inside the 100-year tile; the run condition is "cover ≥ 0.5 for at least 95 % of the analysed cells". Missing
     values (outside tiles, all permanent water) are set to 0 and their share recorded
   - Jobs: comparisons (b10stk, 4 CPU, 16 GiB), hazard map (b10haz, 4 CPU, 14 GiB, resumed after the run condition is
     checked), form AC (b10plasac, 2 CPU, 8 GiB). They are resumed after all B-9 Jobs have been resumed, with at most 40
     B-9 and B-10 Jobs running together
2. **Local (B) implementation** (v66)
   - B2: the excluded forms are removed within the riverine class only, at any matching position (word, prefix, suffix,
     infix). The wet class is unchanged. A re-reading through the same code path that keeps all forms is reported beside
     it to show that it reproduces lex_floodplain_n1. Countries are those of the cell's nearest settlement
   - B4: the grid of the wild cluster score bootstrap runs over log OR −0.80 to 0.80 in steps of 0.02. The restricted fit
     is IRLS with an offset (AOI fixed effects); the efficient score is the weighted residual of the regressor after
     projection on the other regressors and the AOI fixed effects. The interval is the range of accepted points (whether
     it is contiguous is recorded)
   - B7: the imagery class is gfm_summary's sensor (s1 if any post-event image is Sentinel-1, sar if only other radar,
     optical, none if not recorded)
   - B11: p for injection is the corrected code's held-out M1 prediction (M1 uses no names, so it equals the frozen code's;
     the largest difference is recorded). The salts are those of B-6 forms A and B. Injection is on held-out cells (each
     cell once)
   - B12: ties in minimum HAND are broken by mean HAND and then at random (seed 20261030); ties in mean HAND at random
   - B13: the AOI flood rate enters as the empirical logit log((k + 0.5)/(n − k + 0.5)). The composition is the share of
     analysed cells whose nearest name is water-side. AOIs alone in their activation drop out with the fixed effects.
     Standard errors are CR1 clustered by activation, intervals use t(G − 1)
   - B14: the group code sh counts sh, sr, hr, bs and cnr; bantu counts ny, bem, sn, nd, zu, xh, rw, rn, lg, ln, kg, ts, tn,
     st, ss, ve, nr, ki, yao, tum, loz, ng, hz, sw and mg. Units are matched on the name, region and geonunit fields of
     Natural Earth 1:10 m first-level units (zip fetched from naciscdn on 2026-10-29, sha256 efc59726…2ba05); the matched
     units are recorded
3. **Addition (an arithmetic translation, referee 3, 3.1(c))**: a table showing how an effect confined to the ground a
   name refers to is diluted in the odds ratio for "all ground nearest the name", under stated assumptions. Let φ be the
   share of the nearest ground that the name refers to and OR_n the odds ratio there. With the rare-outcome
   approximation, the risk ratio on the ground of truly water-related names is RR = 1 + φ(OR_n − 1), and with misreading
   the observed odds ratio is OR* = (p1 RR + 1 − p1)/(p0 RR + 1 − p0) (p1 = 0.76 and p0 = 0.18 from Table 4; also on the
   grid of B5). φ = 0.05, 0.1, 0.25, 0.5; OR_n = 1.5, 2, 3, 5. It uses no data, so its values cannot depend on any result

## Addendum 2 (2026-10-29T03:20Z, commit e682643, before any result of part C, after all results of B1–B15)

Seen: every result of B1–B15 (v66). No quantity of part C has been seen (still running on the cluster).

1. **A placebo for B10**: the within-AOI R of B10 (the run world_v2_gfm, every model refitted on the GFM labels) was
   0.011 (90 % 0.004–0.019, 95 % 0.003–0.021), so its 95 % interval lies above zero. For such an R the decision table of
   the protocol requires the difference from the placebo (M3np, names reassigned among the settlements of each AOI; the
   not-meaning row), which the plan of B10 did not state. The held-out predictions of M3np from the same run are scored in
   the same way (the logit of each AOI shifted so that the mean prediction equals the observed flooded fraction), and the
   placebo's R and the 90 % and 95 % intervals of R − R_placebo are computed with the same bootstrap draws as B10
   (location groups, 4,000 replicates, seed 20261030), and the row of the decision table is applied. For comparison the
   same quantities are given for the within-AOI evaluation of the corrected EMS run (world_v2_lexfix; item B4 of
   amendment 9, whose R is negative) with the draws of that item (seed 20261029). This addition was decided after the
   result of B10 had been seen, and the paper says so
2. **A description added to B1**: the numbers of location groups, AOIs and countries holding the flooded cells of each
   group within 500 m were added to the odds ratios of B1 (descriptive; no estimate changes; added after the results of
   B1 had been seen)
