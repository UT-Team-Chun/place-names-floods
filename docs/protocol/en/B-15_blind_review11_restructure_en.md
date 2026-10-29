# B-15 Restructuring after the eleventh round of first-time critiques, and checks of the uncertainty of the post hoc upper limit

<!-- English translation drafted with a language model; not yet checked by the authors. The Japanese file
     (docs/protocol/amendments/B-15_blind_review11_restructure.md) is the record and governs where the two differ. -->

Status: **plan (2026-10-29T21:03Z, commit e47fdf6, fixed before any quantity of part B below was computed)**. The confirmatory verdict (equivalent-null) is
unchanged. Everything here is post hoc (exploratory). Protocol: [`preregistration_v2.md`](../preregistration_v2.md)
(FROZEN 2026-10-25T20:58:37Z, commit e625924).

Values seen when this was written: the confirmatory analysis and every result of amendments B-2 to B-14 (including item B11
of B-12 and v75 of addenda 5 and 6 of B-14). None of the quantities of part B below has been computed. The values that
enter the table of item B1 (the metrics of the confirmatory run) have been seen; in particular, the comparison with
flexible models refitted on floods after the language model's training cut-off (R 0.168, 95% 0.039 to 0.328, placebo
0.137, R − Rp 0.031, verdict not-meaning) was read from the metrics after the eleventh critiques. It was in neither the
text nor any table.

## Background

The eleventh round of first-time critiques (NHESS version, manuscript 36 pages and Supplement 36 pages, development
repository at dcfd8cd, public repository at e36b13e; reports in [reviews/](../../research/reviews/),
`2026-10-30_paper2_blind11_*.md`) gave the same decisions as the tenth: the editor returned the manuscript before review
(announcing a recommendation to reject if the contribution is not shown next time), and all three referees asked for major
revision. Points common to the four:

1. **Contribution**: the design cannot detect an effect confined to the surroundings of settlements (ceiling 0.027). Put
   this forward as part of the design rationale rather than as a discovery made after the result, and restate the
   contribution as the dilution caused by assignment through settlements (the ceiling, the records and the positive
   control) and as the bound on this scheme. Compare the names' increment with records assigned in the same way, as paired
   differences (editor M1, referee 1 M1, the conclusions of referees 2 and 3)
2. **Length and structure**: halve the main text; move the amendment history, the census and the 48 versions near
   settlements to the Supplement; add a glossary of the comparisons (all four)
3. **Title and abstract**: "area-wide" contradicts the within-area estimand (the Supplement's title says within-area);
   "water-related terms" is a composite that includes crossings, valleys, qualifiers and wells; the abstract omits the
   pre-specified results that point toward names, the dependence of the estimate on a few places, the design's inability to
   detect effects near settlements, and the fact that no person has checked the readings (all four)
4. **A pre-specified comparison not reported**: the row "Same, flexible models refitted (secondary)" of Fig. 4a (R 0.17)
   is in neither the text nor any table. Report every pre-specified comparison on every split, with its placebo and
   verdict (editor M9, referee 1 M7, referee 2 M6)
5. **Uncertainty of the post hoc upper limit**: the effective number of clusters for the names increment is about 2, so
   normal critical values are not credible; the robustness checks of Table S11 all hold the models fixed and do not apply
   to the refitted 0.031 of level (d); the upper limits range from 0.016 to 0.067 across constructions, and calling one of
   them "the result" was a choice made after all had been seen. Use super-clusters linking location groups and activations
   for folds and resampling (editor M6 and M8, referee 2 M3 to M5)
6. **Hypotheses of the tests**: the protocol's rule is a fixed-margin TOST in Brier units (|ΔBS| < 3.78 × 10⁻⁴); the test
   on the ratio is a test of a different hypothesis, and calling it a "corrected test" conflates the two (referee 2 M2)
7. **Readings**: languages assigned by country, regional and indigenous languages, romanisations that do not match the
   gazetteer, whole-word matching, wells and waterfalls, the composite class, the directions of the element groups, and
   language-model readings that can draw on knowledge of the places (referee 3, editor M4)
8. **Numbers and references**: many (every minor item of the four)
9. **What cannot be checked here**: the toponymic audit, the translations, the time of the freeze, the deposit of the data,
   the policy on generative AI (authors' decisions)

## A. Rewriting (no computation; changes to the paper and Supplement)

1. **Title**: say that the class was set by the lexicon, and name the assignment, the within-area estimand and the mapped
   river flooding; at most 25 words
2. **Abstract**: state the standing of each sentence (the test fixed in the protocol / reanalysis after the result /
   exploratory). Include: the estimand; the protocol's verdict (and that it is a fixed-margin test in Brier units); for the
   post hoc upper limit, the range over the constructions that refit models and the largest value under a reference for
   few clusters (B2, B3); that the design could not detect effects confined to the surroundings of settlements (the
   ceiling); the dilution of records assigned through settlements and the paired differences between names and records
   (B4); the pre-specified results that point toward names (the H4 total increment, and the flexible refit on floods after
   the language model's cut-off together with its placebo); that no person has checked the readings; and that the areas were
   selected for extensive mapped flooding. At most 200 words
3. **Structure**: about half the present main text (at most about 16 pages in the final layout). 1 Introduction (the
   ceiling as part of the design rationale); 2 Data and design (with a glossary table); 3 Results (3.1 the record of the
   protocol's test and every pre-specified comparison; 3.2 the range of post hoc upper limits; 3.3 information assigned
   through settlements; 3.4 what the readings capture (the relief of settlement sites) and their limits; 3.5 near
   settlements (the planned GFM estimate and the EMS estimate on the same cells, the observability of the rings, the rest
   in one sentence)); 4 Discussion; 5 Conclusions. The amendment history, the census, the 48 versions and the details of
   sensitivity analyses and checks go to the appendices and Supplement
4. **Description of the tests**: distinguish the protocol's rule as a test of H_abs (|ΔBS| < 0.05 × the realised
   yardstick, 3.78 × 10⁻⁴ in Brier units) and the test on the ratio as a test of H_ratio (|ΔBS / D| < 0.05); do not call
   the latter a "corrected test". The 32% and 31% are error rates against H_ratio (uniform shift, location groups, models
   fixed); give the 0% under a proportional shift beside them. Put the Y1 check (margin from the lower 95% limit of D, the
   95–95 method) in the main text. Label the two kinds of "TOST" separately in the tables
5. **How the post hoc upper limit is written**: none is called "the result". The one-sided 95% upper limits of the
   constructions that refit models ((b) to (d), amendment 11, pigeonhole, B3) are listed with their assumptions and given
   as a range, with the values under a reference for few clusters (B2). G* is defined as an outcome-based index of the
   concentration of variance (equivalent to a Satterthwaite effective degrees of freedom). The checks of Table S11 are
   described as concerning the interval with models fixed (l. 277–280, 498–500 and 525–526 are corrected). Injected odds
   ratios are said to "correspond to" a limit, with the scale (nominal injection scale) and the reference (the main
   reference) stated; the ranges go to the Supplement
6. **Results that point toward names**: the table of every pre-specified comparison (B1) in the main text or Appendix C.
   The flexible refit on floods after the language model's cut-off is reported in the text with its placebo (the placebo
   improved almost as much, so the verdict was not-meaning). The flexible H1 comparison (0.146) is described. The GFM
   increment in Brier units (5.7 × 10⁻⁴) and R − R(placebo) (EMS 0.005, GFM 0.028) are given as secondary quantities. The
   joint Wald test is presented as the answer to whether names add anything
7. **Scope of the readings**: "water-related" is defined as the composite of five lexicon classes (water, wetland,
   riverine, crossing, lowland) that includes wells, springs, waterfalls, qualifiers and valleys, stated near the title and
   abstract. The element groups are presented as hypotheses. The rows for "names read correctly" (Table 4, Fig. S3) are
   labelled as hypothetical calculations. State that languages are assigned by country, that regional and indigenous
   languages are read with the national block, that romanisations may not match the gazetteer, that whole-word matching
   limits recall, and that the African AOIs were almost unread (citing the comparison without Africa). The date and motive of
   naming are unknown (l. 538 is corrected). The language-model question asks for the referent and can draw on knowledge of
   the places. l. 506–512 are rewritten to say that neither result concerns the meaning of name elements. The relief of
   settlement sites gets a short section
8. **GFM**: drop "designated replication" and describe it as a check of the label processing (also in the heading of Sect.
   4). State that for many AOIs it comes from the same Sentinel-1 scene as the EMS map and is not independent. Correct the
   explanation of the character model's within-area gain ("language and region" cannot explain it, since the reference
   contains the AOI mean). Apply the non-identification of the odds ratios near settlements in both directions
9. **Near settlements**: the main text keeps only the planned GFM estimate and the EMS estimate on the same cells (with the
   intervals of B10) and the observability of the rings (B5). The figure of the 48 versions goes to the Supplement and is
   said to have no inferential standing. Write the estimator as an equation
10. **Numbers and references**: every minor item of the eleventh round (four reports). In particular: the median distance
    instead of "kilometres" at l. 3; "baseline gain" instead of "terrain model's gain"; the description of amendment 1; the
    treatment of cells under non-riverine polygons; the reference for 0.098 in Table C1; "factor of 3"; the two verdicts of
    the centred comparison; Eq. (2) against M1s → M3s; 48% against 49%; the denominator of R 0.060; the description of
    Table S28; the order of the post hoc intervals (0.013 first); "two in total"; -eau under BE; 167 against 139; 149
    against 153; 28 pairs with HAND against 324 with elevation; the character model does not "match"; references for "a
    common way"; the literature search; the protocol's predictions (H2 small) in the text; the positional join to the
    Supplement; the colour range of Fig. 3e; "dashed" in Fig. 4; the axis, bands and crosses of Fig. S5; the labels of Fig.
    S6; H2t in Table 2; minus signs in Table S15; the yardsticks of Table S24 and Table 3; long captions; references
    (OpenAI, the version of Claude, the authors of Ploton, the Supplement's references, URLs); the Supplement's title; the
    four ice-jam AOIs against 154; the description of Trigg et al.; continent counts; 0.5% against 1%; the three values of
    the upper limit with location groups refitted (0.013, 0.008, 0.016); the section number in Table S33; row 7 of Table C2;
    the placebo condition in Table S28; the two meanings of "pigeonhole"; 9.6% against 10.0%; Finnish yli-; Kittilae; the
    strata of Table S9 named by block; the n of κ; the release of readings of indigenous names (CARE); the page of Gelling
    and Cole
11. **References**: cite only verified works for the practice of assigning settlement attributes to the nearest grid
    cells (Thiessen or Voronoi assignment); if none is found, do not call it "common"
12. **Repository wording**: change "pre-registration" in the README and .zenodo.json to "protocol (not registered with a
    third party)"; the frozen file name (preregistration_v2.md) is kept, with a note, because renaming would change its
    sha256

## B. Local computations (script v76; seed 20261104)

1. **Table of every pre-specified comparison**: from the metrics of the confirmatory run (world_v2), for each of the four
   splits (primary, country, time, language-model cut-off): nested (M1c → M3n; lexicon only M3nx; character model only
   M3nl), total (M1s → M3s), flexible (M1 → M3), H1, flexible H1, H4 (M1c → M3r), H4 total (M1s → M3rs) and H4 flexible
   (M1 → M3r), with ΔBS (95%), R (90% and 95%), placebo R, R − Rp (95%) and verdict. No new computation (tabulation only)
2. **Reference for few clusters**: for every construction with saved resamples (the protocol's interval, (a) to (d), the
   two-way bootstrap of amendment 11, the pigeonhole of amendment 14 and B3 below), the one-sided 95% upper limit of R (the
   upper end of the 90% Fieller interval) with three critical values: z (1.645), t(G − 1) (G the number of top-level
   resampling units) and t(ν), where ν is the Satterthwaite effective degrees of freedom from the unit contributions l_g to
   L+ = Δ − 0.05D, ν = (Σγ_g)² / Σγ_g², γ_g = ((l_g − L̂+ n_g)/N)² (location groups on the primary split, super-clusters for
   B3). Also give ν for Δ and D separately. The Fieller variance is the resample covariance of each construction (as in the
   present tables); only the critical value changes
3. **Super-cluster reanalysis**: on the split that keeps activations whole (connected components of location groups linked
   by a shared activation; folds by 106 super-clusters; corrected code), save the stage-2 design locally (`world_eval.py
   --tag world_v2_b15sc_stage2_local --split primary --arms M0,M1 --stacks M1c,M3n --dump-stage2 --group-by-activation`;
   the code is changed so that the dump includes the AOI). A one-stage cluster bootstrap: 500 draws of super-clusters (seed
   20261104 + b), with stage 2 refitted with weights in each draw (the refit of level (b) of amendment 13). Comparisons:
   M1c → M3n (the estimand) and M1s → M3d (names centred within AOIs, no composition; as defined in item B5 of amendment 14).
   Point estimate, L± limits (percentile and normal approximation), 90% Fieller interval, one-sided upper limits (z and
   t(ν)) and TOST p. If the point estimate from the local dump differs from that of the cluster run (world_v2_b7actfix),
   give both
4. **Paired differences between names and records**: for each record of item B6 of amendment 14 (surface water from 2022,
   earlier EMS events, JRC 100-year flooding, the HAND of the settlement site added to M1), on that record's cells and
   resamples (4,000 location-group draws, models fixed, the seed of item B6 of amendment 14), R(names) − R(record assigned
   through settlements) on the same yardstick. Names: the main estimand of the corrected code (M1c → M3n of world_v2_lexfix
   restricted to the same cells). 90% percentile intervals
5. **Observability of the rings around settlements**: on the capped EMS table (as in v55 and v37), compare the cells within
   500 m of settlements with water-related names and of other settlements for lc_built and lc_tree (WorldCover 10 m), the
   share of pixels GFM cannot assess (GFM AOIs) and the share of partly permanent water (AOI fixed effects, distance bands,
   95% intervals clustered by location group). Reporting rule: if every difference is below 0.05 in absolute value with its
   95% interval within ±0.05, the text describes the EMS upper limits as bounds on the association if detection errors are
   the same in the two groups; otherwise it gives the differences and does not call the limits bounds
6. **GFM diagnostics**: (a) the yardstick D, with its 95% interval (4,000 location-group draws), of M0 and M1 trained on
   EMS (main run of the corrected code) scored on the GFM labels of the cells GFM can assess; (b) the median within-AOI AUC
   and the CSI with the mapped area matched to the flooded area, of M1 of the run trained and scored on GFM; (c) the GFM
   within-area names result (R and real − placebo) separately for AOIs with and without a GFM scene on the date of the EMS
   map
7. **Tables of the readings** (Supplement): (i) by block (languages included, countries, number of forms by class and
   position, forms that ever matched an analysed name, cells and settlements matched, share water-related, five most
   frequent forms); (ii) the gazetteer by country (feature codes, settlements per AOI, median nearest ground, share of
   coordinates on whole arc-minutes, share of names whose name and asciiname differ); (iii) N-read classes (including
   unknown) by block and continent; (iv) precision and misses by class for each of the two coders; (v) checks of
   romanisation and rules: numbers of Cambodian names with Phnum, Stoeng or Phumi against Phnom, Stung or Phum and their
   matches; Iranian names containing pol, darreh or rudkhaneh as words; Vietnamese names containing cau, con, song, ben,
   bai or dam as words and their matches; Finnish ala- and yli-; -ham and -hamm by country; Spanish sur ("south") among
   the -sur- names; (vi) analysed settlements whose Wikidata P138 (named after) points to a person or another place, and how
   many of them were read as water-related (the query is saved)
8. **Extension of the census**: apply the rules of the v69 census to the outputs of v70 to v76, adding the interval
   estimates of amendments 11 to 15; BH and BY by family, the number excluding the no-effect value, and the direction (the
   rules of item B10 of amendment 14)
9. **How the Norwegian correction reached a Pakistani group**: in the dump of B3 (corrected code), emulate the frozen code's
   reading (lexicon features of Norwegian names and the composition of their AOIs set to 0). Refit stage 2 with the frozen
   and corrected readings, and decompose the change in coefficients, the change in each location group's Brier
   contribution, and the change in the group that changed most, by feature. Also say whether the emulated frozen R matches
   that of the cluster run (world_v2_actgrp)
10. **Inference for the estimates near settlements**: for the planned GFM estimate and the EMS estimate on the same cells,
    (a) a wild cluster bootstrap (location groups, Webb six-point weights, 9,999 draws, intervals without imposing the null);
    (b) the median odds ratio over 200 permutations of the labels within AOIs (names fixed), to see whether the centre away
    from one under reassigned names is a bias of the estimator itself
11. **Checks**: the width of the interval of the fa stratum of Table S9; the magnitude terciles of amendment 7 (not
    reported); 49% against 48%; 167 against 139; 149 against 153; continent counts (53 against 55, 330 against 328); 9.6%
    against 10.0%; the three values of the upper limit with location groups refitted; the yardsticks of Table S24 and
    Table 3; the four ice-jam AOIs

## C. Cluster computations

None (the stage-2 design of B3 is saved locally).

## D. Authors' decisions and a separate study (not done in this amendment)

1. The toponymic audit (the v69 guide, 300 names), including the blocks that carry the results (Persian and Urdu, Hindi and
   Bengali, Khmer, Burmese and the main European blocks)
2. Checking the English translations of the protocol and every amendment
3. A third-party time stamp of the freeze (the authors decide, since the development repository is not to be shared)
4. Deposit of the code and derived data (DOI), available to referees during review
5. One statement of the use of generative AI, and the authors' contributions
6. Deposit of the conversation and prompts that produced the lexicon
7. How machine readings of indigenous names are released (CARE principles)
8. Article type (research article or brief communication) and journal

## Reporting rules

- Everything is post hoc (exploratory); the confirmatory verdict is unchanged. The protocol's verdict is written as the
  record of the test of H_abs
- No single post hoc upper limit is called "the result"; the range over the constructions that refit models is given with
  the values under a reference for few clusters
- Results are reported whatever their direction
- Rules, units and numbers of replicates are not changed after computing; implementation errors are corrected by addenda

## Addendum 1 (record of the implementation of part B; 2026-10-29T21:16Z, commit 30190d0, script v76, before B1 to B11 were computed)

1. **B2**: (a) redraws the weights of item B1(a) of amendment 13 with the same seed and order (the draw for the
   confirmatory run comes first and is used for the row "protocol run, models fixed"); (b) and (c) use the saved replicates
   (v70_b1_refit_reps_*.npy); (d) adds the Nadeau–Bengio variance and covariance of the fold assignments (item B1d of v70)
   to the resample covariance of (c); the pigeonhole uses v73_b4_pigeonhole_reps.npy; the two-way bootstrap of amendment
   11 saved no replicates, so only its percentile interval is given. Units for ν: location groups for (a) and (b),
   countries (the majority country of each activation) for (c) and (d), the smaller of location groups and activations for
   the pigeonhole. The contributions are the test-side cells of the stage-2 dump of the corrected code on the primary split
   (world_v2_b8stage2_local)
2. **B3**: the super-clusters are the dump's units (location groups linked by activation). M3d centres the name features
   within AOIs separately on the training and test sides (as the d_* columns of world_eval). lbfgs starts from each arm's
   unweighted fit. Resampling weights use seed 20261104 + b. Interval details use the same function as v73 (_levels_detail:
   percentile, normal approximation, TOST p)
3. **B4**: weights are drawn with the seed of item B6 of amendment 14 in the order of its records (a record whose rows
   differ between its two runs is skipped without a draw, as in B6). The names increment is the M1c and M3n of the record's
   own run on the same rows. For the HAND of the settlement site (added to M1, world_v2_b14oracle), the rows are checked to
   equal those of the corrected run, and 4,000 new draws are made with seed 20261104
4. **B5**: the table is v55.table() (capped, frozen reading); cells within 500 m. The share GFM cannot assess is 1 − gfm_v
   (AOIs with GFM only); partly permanent water is f_perm_water of the cells. The difference is the slope of a regression
   centred within AOIs, with CR1 variance by location group
5. **B6**: in (a) the EMS-trained M0 and M1 are the confirmatory run's predictions in the v41 table; 4,000 location-group
   draws, seed 20261104. In (b) the CSI marks as flooded the top-ranked cells of each AOI, as many as its flooded cells; the
   AUC is v70's _auc. In (c) the days are classified as in item B12 of amendment 12 (whether a GFM scene falls on any day of
   the window without its first and last day)
6. **B7**: the settlements are those nearest to an evaluable cell of the uncapped table (with their cells counted); frozen
   lexicon. Words checked: KH (phnum, stoeng, phumi, phnom, stung, steung, phum, kampong, kompong, kaoh, koh, prek,
   preaek), IR (pol, darreh, rudkhaneh, darra, dara, pul, ab, chah, cheshmeh, rud), VN (cau, con, song, ben, bai, dam), PK
   (chah, ab, pul, nala, nullah), MM (myit, kyun, chaung, ywa), Finnish prefixes ala and yli, -ham and -hamm, sur. The
   Wikidata query takes items with P1566 (GeoNames) and P138 whose target's P31 is a human (Q5), a water feature (Q4022,
   Q47521, Q355304, Q23397), a place (Q486972, Q515, Q3957, Q532, Q5084, Q1549591, Q15284) or a mountain (Q8502, Q54050)
7. **B8**: every JSON output v36 to v76 (the census CSVs excluded). Amendments: v68 → 11, v69 → 12, v70 to v72 → 13, v73
   to v75 → 14, v76 → 15. BH and BY within each family; directions by the function of item B10 of amendment 14. Run after
   B3 and B9
8. **B9**: Norwegian names are emulated by setting the lexicon features (lex_*, 14 columns) of rows whose nearest
   settlement is in Norway (cc NO) to 0 and recomputing the composition as AOI means (rows where only the second or third
   nearest settlement is Norwegian are not reproduced). The coefficient change is the fold mean of the standardised M3n
   coefficients
9. **B10**: the table is v41 (cells GFM can assess, frozen reading). The wild cluster bootstrap is β* = β̂ + S⁻¹ Σ_g v_g s_g
   (s_g the location-group sums of the scores with the AOI intercepts projected out, v_g Webb six-point weights), 9,999
   draws, seed 20261104; CR1 intervals with normal and t(G − 1) critical values beside it. Permutations use seed 20261104 +
   500,000 + b
10. **B11**: the strata of Table S9 come from H5_strata in the confirmatory run's metrics; magnitude terciles are searched
    for by word in the verification JSON files; the share of AOIs with a higher AUC is computed from v73_b15_per_aoi.csv
    with two definitions; the remaining checks (167 against 139, 149 against 153 and so on) are traced in the code and
    recorded when the text is corrected

## Addendum 2 (implementation corrections to B7 and B11; 2026-10-29T21:27Z, commit 89f8eeb; before their results were used)

1. **B7 (i)**: for blocks without matches, an empty list of booleans was read as a column selection and the script
   stopped. It was changed to a row selection (`.loc` with a boolean array); no other value changes
2. **B7 (vi)**: the Wikidata query failed twice because the response was cut off (the service's time limit); no result had
   been obtained. It was restricted to the GeoNames identifiers of the analysed settlements (122,128) and sent by POST in
   batches of 3,000, with the classes of item 6 of addendum 1
3. **B11**: the language strata of Table S9 are in the strata of v36, not in the confirmatory run's metrics; the AUC column
   is delta_auc; the "magnitude terciles" are the severity terciles of v55 (in a Supplement table under the name "flood
   size"), and the names are made consistent
