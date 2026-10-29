# B-11 Rewriting and exploratory analyses in response to the seventh round of first-time critiques (NHESS version, 2026-10-29)

English translation, drafted with a language model; the Japanese file `docs/protocol/amendments/B-11_blind_review7_analyses.md`
is the record and governs where the two differ.

Status: **plan (2026-10-29T04:54Z, commit 3709288, before any quantity in B and C below was computed)**. The confirmatory verdict is not
changed. Everything below is **post hoc** and responds to an editor and three referees who read only the PDFs of the
manuscript and its supplement for the first time (language-model agents of a kind that does not load the project's
CLAUDE.md, told that they could consult the public repository at 380e9ca; reports in
`docs/research/reviews/2026-10-29_paper2_blind7_*.md`, commits 68e8ae2, 1536ba0 and c77de62). The paper describes these
analyses as "fixed in this document after the confirmatory result had been seen and before they were computed".
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924). Code: with the correction of amendment B-4
(corrected code). Items C1 to C4 of amendment B-9 and item C4 of amendment B-10 were still running on the cluster; this
amendment does not change how they are computed (their reporting follows A and B3 below).

Values seen when this document was written: all results of the confirmatory analysis and of B-2 to B-10 (except B-9 C1–C4
and B-10 C4), and values that the reports point out or estimate. In particular: item B10 of B-10 (within-AOI evaluation on
GFM labels, positive-negligible) and the placebo of its addendum 2; the results of B-10 C1–C3 and C5; the near-settlement
odds ratio on the uncapped sample (0.91, pairs 0.76–1.15, frozen-code readings, B-8); the odds ratio of 0.57 on maps drawn
from optical images; the joint test of addendum 3 of B-7 (13 coefficients, p = 0.40, which the paper had not printed).
None of the quantities in B and C below has been computed.

## What the four reports have in common

1. What the test answered concerns a scheme, "giving the reading of a name to the nearest 250 m cells"; the title's "area-wide
   information" entered after the result and can be read as information between areas, which the estimand holds fixed. The
   title and abstract should name the scheme and the within-area estimand; the ceiling (that an effect confined to the ground
   near settlements cannot pass the margin) is a property of the design and belongs in the methods, with the statement that
   it was found after the result; the interpretation plan (reading an equivalent-null as a "settled null") was abandoned after
   the result
2. The reference model M1c does not transfer outside the primary split (it beats M1s only where AOIs of one activation fall
   in different folds). Estimands with a transferable reference (names centred within AOIs over M1s; scoring within AOIs)
   should be shown alongside, and both runs of the activation split (frozen code, not-meaning, R +0.018; corrected code,
   +0.004) given in the main text with the reason for the change
3. No expert has checked the readings, and their errors concentrate where the near-settlement evidence comes from (Pakistan,
   Bangladesh). The nature of the gazetteer (GeoNames), ASCII folding, feature codes, rounded coordinates, errors (kallar,
   kacha, ried, ...) and gaps (wæsse, ...) in the lexicon, the assignment of languages and the bias of the etymology sample
   should be stated accurately
4. The near-settlement evidence is presented selectively (the range sentence, the interval types, subsets added later, bounds
   and direct estimates on different scales, and a cap that discards nine tenths of the flooded cells)
5. The pre-specification cannot be verified (unchecked translation, private development history, derived data not deposited),
   and some statements are inaccurate ("a rule fixed in advance" was fixed in amendment 4; water-side is defined in the frozen
   code, not in the protocol; the inclusion criteria are missing from the main text)
6. Length and structure: the few quantities that carry the paper are hard to find; the paper should be focused

## A. Corrections and restructuring (no computation; in the paper and supplement)

1. **Title, abstract, introduction and conclusions**: a title in the form of the scheme and the within-area estimand (for
   example "Landscape terms in settlement names, assigned to the nearest 250 m cells, do not improve within-area predictions of
   mapped river flood extent: a pre-specified test in 500 Copernicus EMS areas"). The abstract states (i) what was tested
   (the scheme, the within-area estimand with the regional composition fixed, the inclusion criteria), (ii) the result (R, the
   achieved equivalence bound, Brier units), (iii) the property of the design (an effect confined near settlements cannot pass
   the margin; found after the result; the interpretation plan was abandoned), (iv) the controls (records that improve the
   predictions when given to each cell fall inside the margin when assigned through settlements, but are not settled on the
   AOIs where the record exists), (v) near-settlement odds ratios marked exploratory, on both EMS and GFM labels, (vi) the
   readings have not been checked by experts, (vii) the settlement-scale question is open. "Area-wide" is not used. The ceiling
   (a label-free quantity) moves to the design section of the methods. The heading of Sect. 3.2 takes the form of the estimand
2. **Overstatements**: "could not have registered any such effect, however strong" → "could not have called it meaningful",
   with the bounds within 500 m (fixed 6.25, folds 7.77, refitted 18.4). "The limit applies to any information attached to
   settlement points" → "information that concerns only the ground near settlements" (the same-event flooded-share control gave
   R 0.51). "The same verdict for information that does predict flooding" → limited to the full sample with the record set to
   zero where absent; on the 227 AOIs with the record it is inconclusive (R 0.081). "Names did not improve on position alone
   (H1)" → in the nested form, with the flexible comparison (R 0.146, 0.004–0.355). In Sect. 7.3, "as toponymy implies" →
   "a check that the readings pick up proximity to water, not evidence of memory". Remove from the Conclusions the unplanned
   quantity (HAND at settlement sites), or mark it exploratory
3. **Selective presentation**: the range sentence for the near-settlement odds ratio either covers every variant (GFM 0.55,
   optical 0.57, river size 1.61) or is limited to "EMS variants of sample and reading". Intervals in Sect. 5 and Fig. 6 use
   the pairs bootstrap, other types go to the Supplement. The GFM subset "cells that GFM can assess in full" is marked as added
   after the GFM result. The abstract's odds ratios are marked exploratory, with the GFM value beside them. Subsets defined
   later (islands without au etc., the floodplain reading) move to the Supplement or are flagged wherever they appear
4. **Factual corrections**:
   - The expectations for element groups were written after the class odds ratios near settlements (addendum 3 of amendment 2:
     wetland, riverine, crossing) had been seen (the wetland group equals the wetland class). Bengali char is not in the island
     group. In the spring and well group, bourne and borne are words for streams
   - "By a rule fixed in advance, the confirmatory run stays primary" → by a rule fixed in amendment 4 (after the confirmatory
     result, before the corrected rerun)
   - "The protocol calls a name water-side" → the code frozen with the protocol (`topsem.world.lexicon.WET_SIDE`, commit
     e625924)
   - "Two-letter affixes that the lexicon's own rules excluded" → the policy stated in the lexicon's header (short forms prone
     to false matches are left out, e.g. English -ey, Vietnamese ho); there is no rule on a minimum affix length
   - "Without access to any name, label or event" → without being given them (its training data may include gazetteers and
     Wikipedia)
   - The rule that calls H4 "uninformative" entered with amendment 9 and is added to the departures in Appendix A. The 85 % is
     from a simulation that injected an effect of the lexicon; there is no power analysis for N-read
   - The explanation of why the placebo favours real names (mismatch with the regional composition) goes into the main text
     (pointing to Table S28 and Fig. S1b). The total increment (H4 total: M1s → M3rs against M3rps) has no composition terms,
     so this mechanism does not apply to it
   - The H5 range includes the continent strata. The coder-agreement figures distinguish the binary water-related code (350)
     from all 462 excerpts. The joint test (p = 0.13) names its model and coefficients. The range over fold assignments includes
     the reassignment of Table S18 (+0.005) or says why not. The rows of Table 3 use one run. The numerical mismatches
     (6.15/6.14; 883/882 and the interval; M1c Brier 0.0430/0.0431; llmcut 0.0497/0.0499; cell counts of Table S53; the
     Groups column of Table S47) are corrected or explained. Y3 in Tables S7 and S8 is "from a constant to M1s". The zero
     coefficient of the indicator "the lexicon reads a language of the country" in Table S14 is explained (whether it is
     constant in the sample). Table S3 adds addendum 2 of amendment 10 and the five critique rounds with access to the
     records. The numbers of superseded tables in the Supplement text are corrected. −0.74 is from 12 pairs. The list of
     amendments matches PROVENANCE in the public repository. "The second" → "Amendment 4". The sentence on delivered and
     monitoring products; the condition under which untyped flooding is treated as riverine (centroid more than 10 km from the
     coast); R −0.062 with continent interactions described as coefficients that do not transfer; the character model is not a
     language model; "from none" → "from zero"; wetland names did not lie low as expected (HAND +0.1 m); Fig. 3 (d, i) and
     Fig. 7b show M1c; the wetland value in Fig. 6 comes from two models; the overlapping dashed lines of Fig. 5 (margin and form
     AC); soften "quantitative tests are few and local"; cite Hawker et al. where FABDEM is used; model identifiers and dates;
     the number of sensitivity analyses (eight runs of seven analyses against the rows of Table C2); how lexicon entries are
     counted (form × position); the wording on "forms that supply most water-side cells"; the point about minor names; the
     different recall weightings of Tables S67 and S68; the Wikipedia headings by language; the 163 digraph matches; what
     "water-side" means in Fig. 3 and Table S60; the riverine class missing from the abstract; the balance of "not local
     knowledge" between abstract and body; "mostly in Europe" in the abstract
   - Restore the joint test of addendum 3 of amendment 7 (13 lexicon coefficients on held-out predictions, p = 0.40) to the
     main text (it was dropped when the paper was shortened in round 6). Add rows for amendments 9, 10 and 11 to Table S4
     (planning status) and update its section numbers
5. **Methods**: the inclusion criteria (at least 400 evaluable cells, at least 20 flooded cells, and others) and the fact that
   452 of the 499 excluded AOIs had fewer than 20 flooded cells (selection on the outcome, which removes small floods); the
   estimand is conditional on AOIs with substantial mapped flooding. The cap (allocated by AOI, training also capped, the
   composition computed on the capped table). The denominator of the flooded fraction is the whole cell (stated for both label
   sets; GFM uses assessable pixels). The nature of the gazetteer (national gazetteers, the NGA GNS, user edits; romanisation;
   the ASCII field; feature class P including PPLX, PPLH, PPLQ and PPLW; 11 % of coordinates rounded to whole minutes). The
   encoding of the JRC maps (flooded fractions at RP10, 100 and 500 and the mean RP100 depth; about 90 m; undefended design
   floods). The list of terrain statistics and the sub-grid. River distances are Euclidean; HydroRIVERS derives from 15-arcsec
   HydroSHEDS. The three nearest settlements are a choice of the protocol (no rationale written). Early stopping uses a random
   tenth of the training cells (not spatially blocked)
6. **One account of the use of language models, in the methods**: the lexicon (Claude Opus 5.5, 2026-10-25, one conversation;
   whether the record of its instructions can be released is for the authors), the second coder (gpt-6-astra, Codex CLI 0.154.0,
   2026-10-27), N-read (Gemma-4-31B-it, prompt in Table S64), drafting of code and how it was checked (tests, independent reruns;
   of the two implementation errors found after the confirmatory run, one was found by a language-model critique), and that the
   authors have not checked the English translation
7. **Structure**: focus. Main text: (1) introduction; (2) data and methods (the design and its ceiling, labels and inclusion,
   names and readings, covariates, models, estimand and decision rule, exploratory analyses and multiplicity); (3) pre-specified
   results (verdict, achieved bound and Brier units, secondary hypotheses in the protocol's terms, how the reference transfers
   and the estimands with a transferable reference, both runs of the activation split); (4) what assignment through settlements
   can carry (ceiling, controls with their limits, bounds on two scales, checks of the assumptions); (5) near settlements (the
   water-side odds ratio on the uncapped sample with EMS and GFM labels, AOIs not mapped from Sentinel-1, detectable effects,
   label errors; groups, floodplain reading and sensors to the Supplement); (6) the readings (what was read, known errors,
   language, limits of the check; short); (7) discussion (scheme and names, limits of readings and labels, what a
   settlement-scale test needs, historical geography and living land terms, limitations); (8) conclusions (no unplanned
   quantities). The amendment chronology (Table A1) moves to the Supplement; Appendix A keeps a short list of departures. The
   Supplement keeps the tables cited by the main text or the Supplement text and moves the rest to "Additional tables"
8. **Terms**: "terrain model/terrain gain" → "baseline model/baseline gain" (terrain, HAND, rivers, land cover, settlement
   layout); "yardstick" → "measure of the baseline gain"; verdict labels (not-meaning and others) glossed at first use
9. **References** (checked against Crossref before they are added): Jones 2016 (wæsse and traditional ecological knowledge),
   Isoda et al. 2019 (exaggerated tsunami place names), Mason et al. 2018 and Chini et al. 2019 (SAR flood detection in urban
   areas), Stephens et al. 2014 (problems with binary pattern measures), Dottori et al. 2022 (JRC European hazard maps),
   Mundlak 1978, Schuirmann 1987. A reference list for the Supplement (Owen 2007, Murphy 1973, Holm 1979, Benjamini and
   Hochberg 1995, and others)
10. **Public repository**: README and PROVENANCE up to B-11; the paper also says that the translation is unchecked

## B. Computations from existing artefacts (local, script v68; seed 20261031)

1. **Change of R on the activation split from frozen to corrected code**: from the held-out predictions of world_v2_actgrp
   (frozen) and world_v2_b7actfix (corrected), split ΔBS(M1c → M3n) by location group and the difference between the runs into
   location groups that contain Norwegian AOIs and the rest. Check whether the M1 predictions of the two runs are identical
   (largest difference) and list the change of the second-stage coefficients where the fold records hold them (descriptive)
2. **Water-side odds ratio near settlements on every evaluable cell, without the cap** (corrected-code readings; the estimator
   of Sect. 5: AOI fixed effects, distance bands, held-out logit of M1, water-side × distance band; pairs cluster bootstrap with
   2,000 replicates; smallest odds ratio detectable with 80 % power): (a) EMS labels, 500 AOIs (conditioned on the M1
   predictions of the uncapped run world_v2_s2a_nocap), (b) GFM labels, every cell that GFM can assess, (c) EMS labels on the
   cells of (b), (d) (b) and (c) on the AOIs whose EMS map is not from Sentinel-1. If GFM labels do not exist for uncapped
   cells, the extent is recorded and (b)–(d) are computed on the capped table. (a) becomes the main estimate of Sect. 5, beside
   the capped value (1.01) (decided after seeing that both include one)
3. **Bounds and margin on the scale of the Sect. 5 estimator**: from item B11 of B-10 (forms A and B, nominal odds ratios 1.25,
   1.5, 2, 3, 20 replicates each), interpolate monotonically on the log scale the curve from nominal to estimator odds ratio
   (relative to the no-injection reference), and translate the margin of form A (2.33) and its bounds (fixed 1.23, folds 1.31,
   refitted 1.57; and the bounds of B-9 C1–C4 and B-10 C4 by the same rule when available). Form B (within 500 m) likewise.
   Values beyond nominal 3 are written as "> the estimator odds ratio at nominal 3"
4. **Test of the GFM domain**: score the held-out predictions of the corrected run (world_v2_lexfix, trained on EMS) on the EMS
   labels of the cells that GFM can assess (the cells where the B-3 GFM labels are defined), and give the baseline gain (M0 →
   M1) and R (M1c → M3n) with a cluster bootstrap over location groups (4,000 replicates), beside the gain of the run refitted
   on GFM (−6.6×10⁻⁴). The number of GFM images per AOI (from the GFM summary, descriptive)
5. **Denominator of the flooded fraction**: (a) the number and share of evaluable cells with 0 < permanent water < 0.5 or
   observed < 1 (overall, within 500 m of water-side settlements, within 500 m of other settlements). (b) A label with the
   denominator "observed area that is not permanent water" (f_observed − f_perm_water): y = 1 ⇔ riverine flooding / denominator
   ≥ 0.5, y = 0 ⇔ flooding of any type / denominator < 0.5, other cells and cells with a denominator below 0.25 left out (counts
   recorded). (c) With this label, the water-side odds ratio near settlements (capped and uncapped, pairs) and R of the primary
   comparison (rescoring the existing predictions of the corrected run; cluster bootstrap over location groups, 4,000
   replicates; the baseline gain on the same labels)
6. **Spread from refitting the base models**: from R under the 13 fold assignments of the primary split (the confirmatory
   assignment and the 12 of B-5; frozen code), a 90 % interval with the corrected variance of Nadeau and Bengio (2003),
   (1/J + n_test/n_train) S² (J = 13, n_test/n_train = 1/9), with t(J − 1). Descriptive; the verdict is unchanged
7. **A combined interval**: extend C2 of B-8 (the second-stage designs of world_v2_b8stage2_local, corrected code) so that each
   replicate (i) refits the second stage, (ii) weights cells two-way (activation × country; the cell weight is the product of
   the two multinomial weights, the units of the protocol's two-way bootstrap), and (iii) computes the denominator (M0 → M1)
   with the same weights. 500 replicates, seed 20261031 + replicate. Report the 90 % interval of R and whether it lies within
   ±0.05. For comparison, the same computation (refit + ratio) with multinomial weights over location groups
8. **Census of exploratory estimates**: scan the JSON files v36 to v67 in docs/verification and collect mechanically, by key
   rules, every estimate of a name effect (R and ΔBS of a names increment; odds ratios of name classes, groups and readings;
   differences in HAND or distance by name class) that has a 95 % interval (record the rules and the excluded keys; baseline
   gains, Brier scores and the information of the controls — GSW, recurrence, oracle, hazard maps — are not name effects and are
   excluded). Each row gives the file, the JSON path, the kind, the estimate, the interval, the amendment, the planning status
   and whether the interval excludes the no-effect value (0 for R and differences, 1 for odds ratios). Report the counts (overall,
   by amendment, by kind), the share that excludes the null value, a histogram of approximate two-sided p values from the
   intervals, and BH q values in the family of name odds ratios near settlements. Rows with the same estimate and interval count
   once. The machine-readable table (CSV) goes to the public repository
9. **Description of the readings**: the nearest settlements of the analysed cells and their water-side matches by GeoNames
   feature code (PPL, PPLA…, PPLX, PPLH, PPLQ, PPLW, …); the share of settlements whose coordinates are rounded to whole minutes,
   by country
10. **Composition in the pilots**: from the run artefacts, the increments (ΔBS and R) of M1s → M3s and of the placebo at δ = 0
    in pilot3 (2 replicates) (descriptive)
11. **A new sample for the toponymic check**: from the frozen lexicon and GeoNames, a sample of names stratified by matched form
    × country, weighted by the number of water-side cells and of flooded cells within 500 m (at most 3 names per stratum, up to
    150 names, from the heaviest strata). The coding guide is rewritten so that coders consult standard place-name dictionaries
    rather than transcribe Wikipedia, and judge homographs and whether the name refers to the settlement site or to a feature
    beside it. The study's readings (the answers) are withheld until the check is complete

## C. Cluster computation (utens, corrected code)

1. **Permutations with recomputed composition on the country split**: the computation of item 1 of B-7 (`--perms-recomp`;
   reference M1 + c(q) and names M1 + c(q) + q) with `--split country`, permutation seeds 1 to 100. Position (percentile) of the
   real R in the permutation distribution and the one-sided p. The H4 total increment has no composition terms and is not
   computed here (A4 gives the reason)

## D. Not done in this amendment (for the authors, or another study)

1. A check by toponymists and language specialists (requested by the authors with the sample and guide of B11) and a rerun with
   a reading restricted to checked forms
2. Checking the translation, depositing the development history with a third party, depositing derived data and code (reserved
   DOI) with a tagged release, releasing the Japanese test, releasing the record of the conversation that compiled the lexicon
   (the authors)
3. The settlement-scale test (flooding per settlement from GFM since 2015, built-up extents, element sets fixed by experts,
   third-party registration): a separate study
4. Rebuilding labels and baseline: detail maps, the lag between flood peak and image, flow-path HAND (FABDEM with thresholds),
   GFPLAIN250m, GHSL built-up layers, FLOPROS, the JRC European maps, early stopping blocked by location, random effects or
   partial pooling, per-name language identification and alternate names, single-image GFM labels. Each would be a new analysis
   close to the primary comparison; they are inputs to the design of the settlement-scale test. The paper states them as limits
5. Injecting through a calibrated probability (rebuilding the plasmode): B3 shows the difference of scale, and the paper states
   the assumption (M1's prediction taken as the true probability) and its effect
6. Changing the primary split (making the activation split primary): the protocol is not changed. Both runs of the activation
   split are given in the main text beside the primary

## Reporting rules

- B and C are all post hoc and were decided after the results of the confirmatory analysis and of B-2 to B-10 (except B-9
  C1–C4 and B-10 C4) had been seen. Quantities judged by the decision table keep the table's rules (verdict names glossed)
- Making B2(a) the main estimate of Sect. 5 was decided after both the capped and the uncapped values had been seen (both
  include one)
- Bounds are given both as nominal odds ratios and as estimator odds ratios
- The census of B8 replaces the paper's "several hundred intervals"
