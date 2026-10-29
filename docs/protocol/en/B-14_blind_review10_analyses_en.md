# B-14 Restructuring and post hoc analyses after the tenth round of first-time critiques

Status: **plan (2026-10-29T18:13Z, commit e530f97; fixed before any quantity of parts B and C below was computed)**.
The confirmatory verdict (equivalent-null) is unchanged. The protocol's verdict is reported as the record, and the post hoc
re-analysis is headed by a one-sided upper bound (non-superiority). Everything here is post hoc (exploratory).
Governing protocol: [`preregistration_v2.md`](../preregistration_v2.md) (FROZEN 2026-10-25T20:58:37Z, commit e625924).
The Japanese file governs where the two differ; this translation was drafted with a language model and has not been checked
by the authors.

Values seen when this was written: the confirmatory analysis and the results of B-2 to B-13, except B9 (a) to (c) and (g) of
B-13 ((a) on the local machine, (b), (c) and (g) on the cluster; the replicates of (c) were saved but B9 was not yet
aggregated) and B11 of B-12 (a penalty on the name features only, being computed locally). None of the quantities of parts B
and C below had been computed.

## Background

The tenth round of first-time critiques (NHESS version, 33-page paper and 33-page Supplement, development repository at
de60b4f, public repository at a6dbf14; reports in [reviews/](../../research/reviews/), `2026-10-30_paper2_blind10_*.md`):
the editor returned the manuscript before review, and the three referees recommended major revision. Points common to the
four:

1. **Headline**: every upper limit is 0.031 to 0.036, and the post hoc re-analysis should be headed by a one-sided upper bound
   (non-superiority). The corrected test was adopted after the interval of amendment 11 (−0.055 to 0.036) had not fallen
   inside the margin, so the two should be given in the order they were defined, and "at every level" should go (editor M4a,
   referee 2 3.1)
2. **Rows (c) and (d) of Table S10**: percentile and normal-theory limits are mixed, and the limits and Fieller interval of
   (c) cannot come from the same variance (heavy tails). State how each row is computed and give TOST p values; the operating
   characteristics at the refitted levels are unknown (editor M5, referee 1 M7, referee 2 3.2)
3. **Cluster inference**: the effective number of influential clusters is about two. Jackknife, studentized bootstrap, the
   effective number of clusters, and a crossed (pigeonhole) resampling of places and events (editor M6, referee 2 3.3)
4. **Reference model**: it depends on a leaky composition term; give the decomposition from M1s and the composition-free
   comparison in the main text (editor M6b, referee 2 3.5)
5. **Positive results in the abstract**: H4's total increment is a between-area result; the placebo of the GFM within-area
   result was added afterwards; the planned GFM odds ratio of 0.55 points the other way. Show the direction or remove them
   from the abstract (all four)
6. **"About a fifth"** fits only two records; give the range with intervals for every record; the JRC maps come from a
   reanalysis that includes the event years (all four)
7. **Reading**: the description of the corrected lexicon does not match what was added and removed (Old English wæsse and
   hamm, Brühl, terp and river names were not added); the etymology sample contains almost none of the countries and forms
   that carry the results; the block mapping, the composite class, the relief at settlement sites (referee 3, editor M7)
8. **Numbers and references**: M1c's Brier score (0.04309 and 0.04305), 33 % and 31 %, two definitions of Y3, eight and nine
   sensitivity runs, 0.84 (0.75 to 1.33), -au 182 and -ab 73, section numbers, missing references (all four)
9. **What cannot be verified**: the policy on generative AI, the freeze record, the translation, the derived data (authors'
   decisions)

## A. Rewriting (no computation; changes to the paper and Supplement)

1. **Title** (at most 20 words, no abbreviation) naming the reading (a lexicon compiled by a language model), the assignment
   (nearest grid cells) and what the design can detect
2. **Abstract** (at most 200 words): what was tested (reading and assignment) → the limit of the design (effects confined near
   settlements cannot reach the margin; the margin corresponds to an injected odds ratio of about 2.3 on all nearest ground) →
   the protocol's verdict (record) → the post hoc one-sided upper bound (R, Brier units, injected odds ratio) → the dilution of
   records known before each event (range of the share kept) → that this verdict should not be cited as evidence against
   reading place names for local hazard awareness → open. State that the protocol was not registered with a third party, and
   "before any model was fitted to the cell-level labels". Remove the positive exploratory results (H4's total increment, GFM
   within areas) from the abstract and report them in the text with their direction
3. **Sect. 3** in two parts: Part I, the confirmatory analysis as specified, the record of its verdict and the defect of its
   rule; Part II, the post hoc re-analysis of the primary comparison. Part II gives the interval of amendment 11 and the test of
   amendment 13 in the order they were defined, drops "at every level", and states the statistical reasons for the choice
   (product weights over activations × countries are nearly nested and count variance twice; percentile intervals of a ratio
   are unstable). It is headed by the one-sided 95 % upper limit (R, Brier units, injected odds ratio), and gives TOST p values,
   whether each limit is a percentile or normal-theory limit, the levels whose operating characteristics are known ((a) and,
   from B3, (c) with models fixed) and those whose are not (the refitted levels). "TOST p < 1/B" in Text S1 is relabelled as
   the protocol rule's value
4. **Reference model** (Sect. 4.1): the decomposition BS(M1s) − BS(M3n) = [BS(M1s) − BS(M1c)] + [BS(M1c) − BS(M3n)] (B5) and
   the composition-free comparison (M1s → M3d) in the main text. Replace the Mundlak analogy by "AOI-level covariates". State
   that the L2 penalty with C = 1 is negligible (close to maximum likelihood) and, as a limitation, that one coefficient per
   feature is assumed across languages
5. **Positive exploratory results**: H4's total increment is a between-area result (the nested increment was equivalent-null),
   and the language model can know the places. The GFM within-area result is given with the facts that its placebo comparison
   was added after the result, that the placebo made predictions worse, and that the planned GFM odds ratio of 0.55 points the
   other way; its direction is reported (B7)
6. **Dilution** (Sect. 5.1): drop "about a fifth"; give the share kept for every record with an interval (B6). Surface water for
   events from 2022 kept nothing. Because the global JRC maps come from a reanalysis that includes the event years, "known
   before each event" is restricted to the earlier EMS events and to surface water for events from 2022; the JRC maps are
   described as a record that does not contain the map of the event. State that this lesson was found after the result
7. **Reading** (Sect. 6.1): call the "corrected lexicon" the "edited lexicon" and list the added entries (fr sur, cs and pl
   nad, de hamm, damm and deich, es rambla, hu fok, km trapeang) and the removed ones in a Supplement table (B8); state that Old
   English wæsse and hamm, Brühl, terp and river names were not added. Give the mapping of language blocks to countries (B8),
   the definition of the composite water-related class, and that its classes have different expected directions. State that
   the etymology sample contains almost none of the countries and forms that carry the results (B8). Quantities that depend
   on the misreading rates (row b of Table 4 and the dilution table of the Supplement) are illustrative in the main text and
   removed from the argument. The "Expected" column of the element groups is post hoc and written without specialists. The
   island group (without au, aue, auen) is reported as hypothesis-generating, with its caveats. The relief difference at
   settlement sites (Supplement Fig. S4a) becomes an observation of Sect. 6, with its overlap with HAND. Qualifier pairs on
   elevation and position along the flow (B9). A note in the Introduction that the Old English examples cannot be read by the
   lexicon. Delete the sentence "the simplest way to turn the advice to read names as warnings into a feature of a gridded
   model"; state in the Introduction what an adequate test would need
8. **Numbers and references**: one run for M1c's Brier score (corrected code); the protocol rule's shares from B1op of B-13;
   one definition of Y3 (the H4 rows of Table S10 are labelled as using the evaluation-side constant); how the sensitivity
   runs are counted (the eight runs leave out the defective first run of analysis 6; analysis 2 has two); -au and -ab counts
   matched to the table (B8); one value for Africa's share; "Sect. 2.5" and "Sect. 6.3"; the numbers missing from Table S34
   (the GFM odds ratio beyond 1 km, the Brier scores of the activation split, the between-area coefficient); "28 % of them" →
   of the 192 AOIs; the sign of GFM's change in Brier score; T1 is a 5 × 5-cell window (about 1.25 km); Table 3 heading (other
   EMS events exclude the mapped event). References: the exact identifier and reference of the second coder, the version of
   Claude Opus 5.5, Hawker et al. (2022), Zanaga et al. (2022), a URL for the GFM PDD, HydroRIVERS, the Supplement's references
   (cite Baugh, Murphy and Pekel; Holm, Wilson, Benjamini–Hochberg, Benjamini–Yekutieli), Berger and Hsu (1996), Hauschke et
   al. (1999), Carter et al. (2017), Davison and Hinkley (1997), Wooldridge (2010). Correct the description of Faccini et al.
   (2017) and Isoda et al. (2019) (not comparisons with observed floods). Give both the content commit and the freeze commit
   and tag
9. **Structure and display**: the record of amendments goes to Appendix A and the Supplement. No two different intervals under
   one name (the "refitted" intervals of amendments 8, 11 and 13 are distinguished). Fig. 4 adds the corrected intervals. Fig.
   1b boxes, a marker-size key in Fig. 2a, 500 m rings and CSI and AUC in Fig. 3, row labels of Supplement Fig. S6, the type
   size of Table S3, citations of Supplement Figs. S1, S3 and S4, row labels of Table 2, thousands separators
10. **Description of the labels**: the number of AOIs in which overview maps were preferred to detail maps, the treatment of
    generation 0 and 1 AOIs as observed throughout, the different denominators of EMS and GFM (whole cell and assessable
    pixels), the share of cells kept by the cap (about 8 %). The positive control for H3 (C2) and the country composition for
    H4 (C3)

## B. Local computations (script v73; seed 20261103)

1. **Consistency of the corrected test**: for each level (a) to (d), the L± limits both as percentiles (as planned) and in the
   normal approximation (variance of the replicates); the TOST p value (the larger of the two one-sided normal p values); the
   one-sided 95 % upper limit of R (upper end of the 90 % Fieller interval) in R, Brier units and injected odds ratio (curves of
   B2 of B-13). The interval of amendment 11 (two-way, refitted) in the same table. Monte Carlo error of the percentile limits
   from 1,000 redraws of the saved replicates
2. **Robustness to heavy tails** (models fixed, primary split of the corrected run): (i) jackknife variances leaving out one
   location group, and one activation, at a time, and from them one-sided 95 % limits of Δ, D and L± and an interval for R
   (normal approximation); (ii) studentized bootstrap limits of L± (4,000 location-group replicates, each studentized with the
   robust variance from the location-group sums); (iii) the effective number of clusters G* of Carter, Schnepel and
   Steigerwald (2017) for the location-group contributions to ΔBS, on the primary, activation-merged and country splits; (iv) R
   without the two most influential location groups (top two of B6 of B-13), R with location groups weighted equally, with AOIs
   weighted equally and with activations weighted equally (Fieller intervals from 4,000 location-group replicates)
3. **Operating characteristics** (models fixed, activations resampled within countries in two stages, countries drawn as in
   addendum 6 of B-13): the double resampling of B1 of B-12 (1,000 outer × 500 inner), uniform and proportional shifts, the true R
   placed at −0.05, 0 and +0.05; the shares in which the L rule declares equivalence and in which the upper-limit rule (one-sided
   upper limit of L+ below zero) declares non-superiority. Binomial Monte Carlo errors for every share
4. **Pigeonhole resampling** (Owen 2007): independent multinomial weights on location groups and on activations, cell weight =
   product; second stage refitted (design and codes of B1(b) of B-13), 500 replicates, seed SEED + 200,000 + b. L± limits
   (percentile and normal) and the Fieller interval
5. **Decomposition and composition-free comparisons**: for each split of the corrected run with predictions (primary,
   activation-merged, country, time, language-model cut-off), BS(M1s), BS(M1c), BS(M3n) and the decomposition; the corrected test
   at level (a) (4,000 location-group replicates) for the composition-free comparison M1s → M3d and for the total M1s → M3s
6. **Share kept**: for each record with both versions in one run (surface water for events from 2022, earlier EMS events, JRC
   100-year, EMS events before or after (run of B-7)), ΔBS (through settlements) ÷ ΔBS (per cell) on the same replicates (4,000
   location groups, 90 % percentiles; replicates with a denominator at or below zero are counted and reported). For a record
   whose two versions are in different runs with the same cells (surface water, all events), cells are matched first
7. **Direction**: held-out coefficients of the name features with 95 % intervals (cluster-robust by location group; estimator
   of B8 of B-13) on the EMS labels (corrected run) and on the GFM labels (planned GFM table, with the held-out M1c prediction of
   the run trained and scored on GFM labels as the baseline)
8. **Toponymic tables**: (i) the edits of the edited lexicon (language, position, class, entry); (ii) the mapping of blocks to
   countries; (iii) the number of analysed settlements containing omitted elements (-was, -ham, brühl, terp or wierde, -sur- or
   nad, and settlements recorded in Wikidata as named after a river); (iv) matches of char, -au and -ab by country; (v) the
   deposited codes of both coders by country and stratum, and precision and miss rates by block and by the top 20 forms (with n
   and Wilson 95 %); (vi) a gazetteer table by country: shares of feature codes, settlements per AOI, the median area of ground
   nearest a settlement (cells × 0.0625 km²), the share of coordinates on whole arc-minutes
9. **Qualifier pairs**: for the pairs of B15 of B-13, the difference in site elevation (mean GLO-30 elevation of the cell
   containing the settlement point) and the ratio (lower / upper) of the drainage area of the nearest river (the river size of
   B-8). Sign tests and 90 % location-group intervals
10. **Census**: the direction (for or against names) of the estimates that exclude no effect, by family; the family and the BH
    q of the highland, GFM, non-Sentinel-1 estimate (2.45) in the final census
11. **Check of 0.84**: the distribution of the pairs replicates of C2's uncapped EMS estimate (median, skew on the log scale,
    share above one)
12. **Units of the margin**: the margin and the one-sided upper limit as shares of the Brier score of the climatology (training
    constant), i.e. as increments in Brier skill
13. **Checks of the estimates near settlements** (estimator of Sect. 5, capped EMS table): (i) a fractional logit with the
    flooded cover share from the 25 m maps as the outcome (same design, cluster-robust by location group); (ii) AOIs with a
    flood type only; (iii) with log(1 + population) of the nearest settlement and the 10 m built-up share added; (iv) a
    leave-one-location-group-out jackknife interval for the planned GFM estimate
14. **Alignment of the positional join**: for the 500 analysed AOIs, whether cells, features and landcover have the same number
    of rows, and whether the distance to the nearest settlement P_d1 in features is reproduced from the cell coordinates and the
    GeoNames coordinates (a row-count check is added to the code, with a test; if an AOI fails, the primary comparison is also
    given without it)
15. **Per-AOI distributions**: per-AOI changes in Brier score and ΔAUC by continent and flood type; R with each of the five most
    influential activations (activation jackknife of B2) left out in turn (models fixed)

## C. Cluster runs (primary split, corrected code, the default five seeds)

1. **EMS labels without unobservable ground**: cells whose 10 m built-up plus tree-cover share (lc_built + lc_tree in features)
   is 0.5 or more are left out of evaluation (`--label-variant observable`). The primary comparison with its verdict, and the
   near-settlement odds ratio on that table (estimator of C2 of B-13, capped table restricted by the same rule)
2. **Positive control for H3**: site HAND assigned through settlements like the names, at T0 (new arms T0co → T0o) and at M1
   (M1co → M3o); a run with `--oracle`
3. **Country composition for H4's total increment**: on the language-model cut-off split, country indicators in the reference
   (M1k → new arm M3rk), beside M1s → M3rs

## D. Authors' decisions and a separate study (not done in this amendment)

1. The policy on generative AI (NHESS does not allow it for text or interpretation): which parts a language model drafted, and
   what the authors wrote, checked and changed
2. Checking the translation of the protocol and of every amendment
3. Depositing the freeze record with a third party (development history and cluster records); the development repository is not
   to be shared, so the authors decide
4. Depositing the derived data (to referees only during review)
5. The toponymic audit (300 names); its allocation should include the forms and countries that carry the results and the names
   of EMSR052 and EMSR130
6. Depositing the conversation and prompt that compiled the lexicon
7. For the settlement-scale study: finer grids or fractional labels, detail maps, rings outside built-up footprints, assignment
   on the same bank and catchment, assignment to territories, multi-event recurrence as the outcome, the European JRC maps, HAND
   thresholds on a bare-earth DEM, defences, image timing, regional blocks, river names, expert databases, and a reading of
   Japanese aza names by toponymists
8. Article type (brief communication) or choice of journal

## Reporting rules

- Everything is post hoc (exploratory) and the confirmatory verdict is unchanged. The protocol's verdict is reported as the
  record, and the post hoc re-analysis is headed by the one-sided upper bound
- Results are reported whatever their direction
- Rules, units and numbers of replicates are not changed after computing; implementation errors are corrected in addenda

## Addendum 1 (record of the implementation of C; 2026-10-29T18:19Z, commit e738e7c, before C1 to C3 were computed)

1. **C1**: the variant `observable` of `topsem.world.dataset.riverine_labels`. The label (y) is unchanged; cells whose lc_built +
   lc_tree in features (missing values taken as 0) is 0.5 or more are left out of evaluation (the cap draws from the remaining
   cells; inclusion of AOIs is unchanged). Run world_v2_b14obs (`--label-variant observable --arms M0 --stacks M1c,M3n,M3np`,
   primary split, 10 folds)
2. **C2**: arms T0co (T0 + AOI mean of the site HAND) and T0o (T0co + site HAND), run only with `--oracle`. Run
   world_v2_b14oracle (`--oracle --arms M0 --stacks M1co,M3o,T0co,T0o`). Both R values are divided by the primary baseline gain
   (BS(M0) − BS(M1))
3. **C3**: arm M3rk (base M1 + country indicators + N-read features), run only when N-read is given, with the country indicator
   columns built as for M1k of B-10. Run world_v2_b14h4k (`--split llmcut --nread <N-read parquet> --arms M0 --stacks
   M1s,M3rs,M1k,M3rk`)
4. **Code check**: `load_aoi_table` stops when cells, features and landcover have different numbers of rows (B14); three tests
   were added

## Addendum 2 (record of the implementation of B; 2026-10-29T18:37Z, commit 38d41c7, script v73, before B1 to B15 were computed)

1. **Common**: "the corrected run" is world_v2_lexfix (primary split, 10 folds). Yardstick D = BS(M0) − BS(M1),
   Δ = BS(M1c) − BS(M3n), L± = Δ ∓ 0.05·D. Location groups are the unit of the predictions; an activation is the first part
   of the AOI name (the EMSR number); an activation's country is the country holding most of its cells (amendment B-13,
   addendum 6). Random numbers restart from seed 20261103 in each item (B4 alone uses SEED + 200,000 + b as planned).
   Normal-approximation limits are one-sided 95% (z = 1.645)
2. **B1**: the 4,000 weights of B-13 B1(a) are redrawn with the same seed (20261102) in the same order (after the
   confirmatory run) and used after checking that they reproduce the saved interval. (b)(c) use the 500 saved replicates.
   (d) is the variance of (c) plus the Nadeau–Bengio variance, so it has a normal approximation only (no percentile limit).
   TOST p = max(Φ(L̂+/se+), 1 − Φ(L̂−/se−)), with se the replicate standard deviation. The one-sided 95% upper limit of R is
   the upper end of the 90% Fieller interval; in Brier units it is multiplied by D̂; injected odds ratios use the curves and
   the four no-effect references of B-13 B2 (the v70 computation is called). The Monte Carlo error of a percentile limit
   is the standard deviation of the limit over 1,000 resamples (with replacement) of the saved replicates. For amendment
   11's two-way interval only the percentile interval of the ratio was saved (not the replicates); that interval is shown
   and the L± columns are left empty
3. **B2**: (i) jackknife variance (G−1)/G·Σ(θ₋g − θ̄)² from the leave-one-out estimates (Σ − contribution)/(N − n); L± limits
   L̂ ± 1.645·se; Fieller interval for R from the jackknife variances and covariance; over location groups and over
   activations. (ii) 4,000 location-group replicates; the variance of L* in each is the robust variance that counts each
   copy of a drawn group as a separate cluster, G/(G−1)·Σ w_g (l_g − L* n_g)² / (Σ w_g n_g)², and t* = (L* − L̂)/se*. Upper
   limit of L+ = L̂+ − q0.05(t*)·sê+, lower limit of L− = L̂− − q0.95(t*)·sê−. (iii) G* = G/(1 + Γ), Γ = mean((γ_g − γ̄)²)/γ̄²,
   γ_g = ((e_g − Δ̂ n_g)/N)² (the formula of Carter, Schnepel and Steigerwald 2017 with the observed location-group
   contributions to the estimator of ΔBS). Primary split: corrected run. The activation-grouped split (world_v2_actgrp) and
   the country split (country split of world_v2) have no corrected run, so the confirmatory code's runs are used (with the
   confirmatory primary split alongside for comparison). (iv) the two groups with the largest |ψ| for the corrected run in
   B-13 B6. Equal-weight R: the ratio estimator with cell weights 1/(cells of the group, AOI or activation); Fieller
   intervals from 4,000 location-group replicates
4. **B3**: activation sums of the corrected run (fixed models). Outer: countries with activations drawn with replacement
   (as many as there are), then within each drawn copy of a country its activations drawn with replacement (as many as it
   has). Inner: the copies of countries in the outer sample drawn with replacement, then within each copy the activations
   it holds drawn with replacement (the same two stages). Shifts as in B-13 (proportional: e − ê + r0·d; uniform:
   e − ê + r0·d̂; ê and d̂ from the full sample). The L rule and the upper-only rule need no condition on the sign of the
   denominator, so no replicate is set aside (the share of replicates whose denominator has a 95% lower limit at or below
   zero is reported separately)
5. **B4**: the design and signs of B-13 B1(b)(c) (stage-2 dumps of B-8 C2, same starting values). In each replicate the
   location-group weights (multinomial) and then the activation weights (multinomial) are drawn from the same stream; the
   cell weight is their product. Replicates are saved to data/verification/v73_b4_pigeonhole_reps.npy. The point estimate
   is the unweighted refit (as in B-13)
6. **B5**: in each split of the corrected run with predictions (primary, language-model cut-off): BS(M0), BS(M1), BS(M1s),
   BS(M1c), BS(M3n), BS(M3d), BS(M3s) and the increments divided by the yardstick; corrected test (a) (4,000 location-group
   replicates) for M1s → M3d and M1s → M3s. For reference, the confirmatory code's splits (activation-grouped, country,
   temporal) in the same table with the arms that have predictions (M3d from world_v2_dm and world_v2_dm_actgrp)
7. **B6**: ΔBS by the arms compared in each record (through settlements M1cX → M3oX, directly in the cell M1cY → M3oY).
   Surface water from 2022 = G (and Gn) versus W of world_v2_b13gsw22; earlier EMS events = E versus Ec of
   world_v2_b13recprior; JRC 100-year = J versus Jc of world_v2_b13jrcset; EMS events before and after = E versus Ec of
   world_v2_b7recur; surface water of all events = G (and Gn) of world_v2_b6gsw versus W of world_v2_b8gswcell_local (joined
   after checking that aoi and y agree row by row). Interval of the share: 5% and 95% points of the replicates with a
   positive denominator, with the share of replicates whose denominator is zero or negative alongside
8. **B7**: EMS: the name coefficients and 95% intervals of the B-13 B8 estimator (held-out cells of the stage-2 dumps of
   B-8 C2, standardised [base logit, composition, name features], cluster-robust by location group), and alongside the form
   with the logit of the corrected run's held-out M1c prediction as the base ([logit M1c, name features]). GFM: the planned
   GFM table (v41 table + GFM labels), y = GFM label, base = logit of the held-out M1c prediction of world_v2_gfm (cells
   joined by aoi, longitude and latitude), name features = the dumped values of those cells (joined to the dumps by row
   position in the v41 table). Coefficients are log odds per standard deviation of the standardised name feature
9. **B8**: (i) the additions and removals of v70_lexicon_corrected_record.json counted with the class of the original
   lexicon; every item in docs/verification/v73_lexicon_edits.csv. (ii) country_langs of the lexicon. (iii) ASCII names of
   the analysed settlements (nearest settlements in the uncapped table) split into lower-case words: -was, -ham (including
   -hamm), terp and wierde at the end of a word; brühl as a word containing bruhl or bruehl; sur and nad as words. Wikidata
   is queried once through the public SPARQL service (items with GeoNames ID P1566 whose P138 points to a river Q4022 or a
   subclass) and matched (left empty, with a note, if the query fails). (iv) for the frozen and the corrected lexicon, the
   number of settlements whose matched entry has form char, au or ab (any position), by country. (v) coded rows of the two
   published coders (v70_two_coders.csv); truth = any of the five water-side classes in the codes. Precision = share true in
   the wet stratum, miss rate = share true in the not stratum (each coder), by block and by the 20 most frequent matched
   forms of the wet stratum. (vi) by country of the analysed settlements: shares of feature codes, median settlements per
   AOI, median number of evaluable cells nearest to the settlement × 0.0625 km², share on whole arc-minutes (rule of v37)
10. **B9**: the pairs of B-13 B15 (v70_b15_pairs.csv) joined by gid to the settlement river table of v57 (elev_site = GLO-30
    mean of the cell holding the settlement point; upland_nearest = drainage area of the nearest reach); settlements in
    several AOIs take the median. Difference = lower − upper; ratio = lower / upper (tested on the log). Sign test and
    location-group bootstrap (2,000) 90% intervals
11. **B10**: rows of v69_estimates.csv (the final census) that are unique, valid, about a name effect and exclude the null.
    Direction: in the flood family, positive R and ΔBS favour names; odds ratios and differences favour names when above the
    null for water-side classes and below it for the high class. In the site family, wet-side quantities (low_hand,
    low_half, flooding, water occurrence occ and seas) favour names when larger for water-side classes, dry-side quantities
    (hand_mean, hand_pct, log_river, elevation) when smaller, reversed for the high class. Nuisance rows get no direction.
    Rows that cannot be classified are counted and their paths listed. The 2.45 row is found by its path
12. **B11**: the 2,000 saved log odds ratios (data/verification/v68_pairs_v71_c2_uncapped.npy)
13. **B12**: climatology = Brier of the training-side constant of B-12 B17 (primary split)
14. **B13**: the capped table of v55 and the v55 estimator (AOI fixed effects, distance bands, logit of the frozen run's M1,
    cluster-robust by location group), frozen lexicon reading (as for the capped estimate of Sect. 5). (i) response =
    riverine cover (other flood types added for untyped AOIs, as in v55 item 6); (ii) untyped AOIs
    (unspecified_as_riverine) left out; (iii) covariates log(1 + population of the nearest settlement) (GeoNames) and
    lc_built (10 m), with the estimate without them on the same rows (rows with missing values left out); (iv) the planned
    GFM estimate (v41) refitted without each location group, 95% interval of the log odds ratio from the jackknife variance
15. **B14**: row counts of cells, features and landcover in the 500 analysed AOIs (selection of v55). P_d1 recomputed as the
    distance between the cell's x and y in cells and the GeoNames coordinates of the P_gid1 of the same row of features
    projected into the AOI's LAEA; |difference| ≤ 1 m counts as agreement (rows with P_gid1 = −1 only counted)
16. **B15**: by AOI of the corrected run: ΔBS = BS(M1c) − BS(M3n), ΔAUC = AUC(M3n) − AUC(M1c) (AOIs with both classes). The
    continent of an AOI is the most frequent continent of its cells; flood type is typed or untyped
    (unspecified_as_riverine). Median, quartiles, share positive. R without each activation: the five largest |ψ| of the
    activation jackknife of B2(i)

## Addendum 3 (record of the aggregation of C; 2026-10-29T18:53Z, commit dccd8f1, script v74, before C1 to C3 were aggregated)

1. Decision table with the function of the B-13 aggregation (v71): _verdicts of v56 (4,000 location-group replicates,
   fixed models, protocol rule), with corrected test (a) alongside (4,000 location-group replicates, seed 20261103, L± and
   the Fieller interval)
2. **C1**: comparison M1c → M3n (placebo M3np). The odds ratio near settlements uses the estimator and the corrected-lexicon
   reading of B-13 C2 on the capped table of v37, restricted to cells whose lc_built + lc_tree in the AOI's features
   (missing values taken as 0) is below 0.5 (joined by (aoi, ix, iy)); the frozen-lexicon reading is given alongside for
   comparison with the capped estimate of Sect. 5
3. **C2**: T0co → T0o and M1co → M3o, both divided by the primary yardstick (addendum 1)
4. **C3**: on the language-model cut-off split (one fold), M1s → M3rs and M1k → M3rk; the yardstick is BS(M0) − BS(M1) on
   that split

## Addendum 4 (descriptions added after the result of B4 had been seen; 2026-10-29T19:44Z, commit fcb4112; estimates and verdicts unchanged)

The 90% Fieller interval of B4 (pigeonhole, refitted) was −0.096 to 0.067; its percentile limits declared equivalence and
its normal approximation did not (TOST p 0.13). After seeing this, we add the following descriptions to interpret it (not
planned; the paper says they were added after the result):

1. The overlap of location groups and activations: the share of activations that lie within one location group, the
   share of location groups that hold one activation, the number of distinct pairs, and the corresponding shares of cells
   (corrected run, primary split)
2. The number of refitted B4 resamples whose ΔBS has a robust z (median and MAD) above 4, and the standard deviation
   without them
3. The paper's headline gives the upper limit of (d), 0.031, alongside that of B4, 0.067 (B4 is not dropped); the largest
   upper limit of the other levels and versions is given separately, without B4

## Addendum 5 (versions of the near-settlement odds ratio without EMSR052 and EMSR130; 2026-10-29T19:56Z, commit 85ef317, script v75, before they were computed)

The tenth round of first-time critiques (referee 1, item (f); referee 3, item (d)) asked for the water-side odds ratio
within 500 m of settlements without activation EMSR052 (Cambodia) and without EMSR130 (Myanmar). The caption of Fig. S6 of
the Supplement said that "versions that leave out one activation at a time" were not shown, but the only version in the
census of estimates (v69_estimates.csv) that leaves out an activation is one that drops EMSR097 from another class of
names (the riverine class); no such version exists for the water-side odds ratio. We correct the caption and compute the
following (not planned; it answers critiques made after the results were seen, and the paper calls it exploratory).

1. Versions: (a) the two versions that amendment 11 made the main ones, EMS labels on all evaluable cells without the cap
   (B2 ems_all of amendment 11) and GFM labels on all cells GFM can assess (B2 gfm), with the corrected code's readings;
   (b) the GFM estimate planned in amendment 3 (the v41 table, capped, frozen readings; refitted as in item B13 (iv) of
   this amendment)
2. For each, three versions: without the AOIs of EMSR052, without those of EMSR130, and without both. If a table holds no
   AOI of an activation, the version is reported as not applicable
3. The estimator is that of each version (AOI fixed effects, distance bands, logit of M1). Intervals are 95% Wald
   intervals clustered by location group (as in v55.odds; no pairs bootstrap). For comparison the version without
   exclusions is given with the same Wald interval. Also the number of AOIs removed, the cells kept and the flooded cells
   within 500 m of water-side names
4. Reporting: added to the Supplement table of B13, and the caption of Fig. S6 is corrected. Sect. 6 of the paper gives the
   range of the point estimates without the two activations in one sentence. The interpretation and the headline do not
   change with the result
5. Script v75_b14_leaveout.py; output docs/verification/v75_b14_leaveout.json

## Addendum 6 (record after the result of addendum 5 had been seen; 2026-10-29T20:01Z, commit 54dfa87: the premise was wrong, and the form of reporting; estimates unchanged)

1. The premise of addendum 5 was wrong. Versions of the uncapped EMS estimate (all evaluable cells) without each activation
   in turn (157 activations, points 0.86 to 0.97) and without both EMSR052 and EMSR130 (1.10, 95% pairs-bootstrap interval
   0.83 to 1.41) had been computed in item B9 (b) of amendment 13 (v70_b13_local2.json) and were given in Sect. 6 of the
   paper (NearChecksText) and in the Supplement table of amendment 13 on the estimates near settlements. They postdate the
   v69 census, so the census does not hold them, and the caption of Fig. S6 ("not shown") meant that the figure does not
   draw them. We overlooked this when we wrote addendum 5
2. v75 reproduced them (without EMSR052 0.932, without EMSR130 0.972, without both 1.098, equal to the values of amendment
   13). The only new information is that neither GFM table (all cells GFM can assess; the capped estimate planned in
   amendment 3) holds an AOI of either activation
3. Reporting replaces item 4 of addendum 5 (to avoid duplication, whatever the result): the GFM column of the two rows of the
   amendment-13 table says that the GFM tables hold no AOI of either activation; the rows with Wald intervals are not added
   to the amendment-14 table (they stay in the output of v75); the caption of Fig. S6 points to the amendment-13 table; and
   the existing sentence of Sect. 6 adds that the GFM estimates do not include these two activations
