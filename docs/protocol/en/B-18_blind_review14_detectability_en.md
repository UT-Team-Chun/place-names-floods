# B-18 Restructuring after the fourteenth round of first-time critiques: the protocol's test first and briefly, the design's detectability as the subject, and corrections of a wrong statement on upper limits and of inconsistent numbers

> English translation drafted by a language model after the analysis and checked by the authors on 2026-11-06. The Japanese
> file (`docs/protocol/amendments/B-18_blind_review14_detectability.md`) is the record and governs where the two differ.

Status: **plan (2026-11-02T01:10Z, commit bba26f7; fixed before the quantities of B were computed)**. The confirmatory verdict (equivalent-null) is unchanged.
Everything here is post hoc (exploratory). Protocol: [`preregistration_v2.md`](../preregistration_v2.md) (FROZEN
2026-10-25T20:58:37Z, commit e625924).

Values seen when this was written: the confirmatory analysis and every result of B-2 to B-17, and the values the four
reports of round 14 read from the tables. The quantities of B below have not been computed (of the Fieller intervals of
B1, those for earlier floods and the JRC map on all ground and within 500 m were computed in B9 of B-17).

## Background

Round 14 (NHESS version, 31-page paper and 43-page Supplement, development repository d693dda, public repository
9f7c5d1; reports in `docs/research/reviews/2026-11-02_paper2_blind14_*.md`): the editor would not send the paper to
review in its present form (a new submission is welcome); the three referees recommended major revision (the statistician
on the boundary with rejection). Common points:

1. **Policy on generative AI** (editor M1): NHESS does not allow generative AI for text or interpretations; this alone
   prevents review (authors; D1).
2. **What the paper is** (editor M2 (i); referee 2, minimum requirement 1; referee 3, alternative): the confirmatory test
   cannot discriminate. The supported contribution is the calibration with records of known value, which shows that
   allocation through settlements dilutes even such records. Make it the subject and report the confirmatory test first
   and briefly (editor M6 (e): the interpretation plan said the Results would open with the primary verdict).
3. **Wrong statement on upper limits** (editor M2; referee 2, M6 (c)): the coverage check of B1 of B-17 concerned the
   protocol's construction (location groups, models fixed). The sentence "the only upper limit whose coverage was checked
   is the 0.141" transferred it to the designated construction (super-clusters, second stage refitted). The checked
   construction's limits are 0.003 (z) and 0.007 (t(ν − 1), frozen reading; none with the corrected reading).
4. **Exploratory results in the abstract and conclusions** (editor M7; referee 2, M12 (c); referee 3, M8 (d)): the
   results after the language model's cut-off, on the GFM labels scored within areas, on the split that keeps activations
   whole under the frozen reading, and of centred names are not reliable enough for the abstract or the conclusions.
5. **Inconsistent numbers** (editor M8, 13 items; minor points of referees 2 and 3).
6. **Protocol text against frozen code** (editor M6 (c); referee 2, M8 (a)): the protocol's text does not say that the
   bootstrap holds fitted models fixed, that the placebo is drawn once, or whether the placebo rebuilds the regional
   composition. The frozen code holds the models fixed, and its placebo keeps the real names' composition, the scheme
   withdrawn in amendment 5 from the permutation distribution because it favours the real names.
7. **GFM** (referee 1, M3; referee 2, M11 (f)): not independent; for each result it is unclear whether the models were
   trained on EMS or on GFM labels; the GFM-trained baseline is worse than position alone (D < 0); scoring shifted within
   areas uses the evaluation labels.
8. **Estimates near settlements and map production** (referee 1, M4): they change strongly with sensor, product
   generation, hydrography layer and delay, and the text does not say so.
9. **Readings** (referee 3, M1 to M3 and M8): some blocks are composite (fa is Persian/Urdu/Pashto/Sindhi, etc.), and the
   text misdescribes them; assignment by country; matching on ASCII strings without diacritics; different romanisation
   conventions. The dry-point interpretation does not fit the paper's own data (lower site HAND for water-related names,
   a deficit that depends on the sensor). The qualifier-pair sentence is in no table and conflicts with the Supplement's
   HAND result.
10. **What cannot be checked**: the time of the freeze, the translation, the toponymic audit, label QA and data deposit
    (authors; D).

## A. Rewriting (no computation)

1. **Structure**: Results 3.1 sample → 3.2 the protocol's test (a brief record: estimate and interval, verdict, the odds
   ratio the margin corresponds to, the rule's operating characteristics) → 3.3 what the design can detect (the subject:
   ceiling, dilution of known records, injected effects) → 3.4 evidence within areas after the result (randomization
   inference as the designated post hoc test; limits of each construction as descriptions) → 3.5 what the readings
   capture → 3.6 labels near settlements (including the dependence on map production). The Methods follow the same order.
2. **Title**: the design's detectability first, the test of names as the case, without overstatement ("insensitive").
   Proposal: "Nearest-settlement allocation dilutes settlement-scale information in gridded flood models: a protocol-based
   test with machine-read place names".
3. **Abstract** (at most 200 words): question → data → the confirmatory test (R and 90 % interval; the verdict as a
   record; the margin corresponds to an odds ratio of about 2.3 on all Voronoi ground) → detectability (ceiling within
   500 m, share of value kept by known records, the JRC map allocated through settlements also equivalent, injected odds
   ratios of 1.25 and 1.5 also equivalent) → after the result: real names do not stand out among reassignments within
   areas → what the test does not answer (minor names, early forms, indigenous names, local knowledge) → the lesson. No
   exploratory result in either direction; no bare p values or coined verdict labels.
4. **Upper limits**: correct the statements on coverage (background 3). State that no construction has checked
   properties (the simulations vary less than the data and do not reproduce one group carrying most of the variance), and
   give each construction's upper end in Table 5 as a description, with D̂/SE(D̂) and the t critical value ("none" comes
   from the uncertainty of the denominator). Do not convert upper ends into injected odds ratios (removed from the text).
   Remove TOST p values with normal critical values (Text S1).
5. **Protocol text and frozen code** (background 6): add to Table 2 a row on what the frozen code decided and the
   protocol's text did not (models fixed in the bootstrap; one placebo draw that keeps the real names' composition), and
   state in the Methods that this placebo favours the real names. Replace "real minus placebo" in the text by the rank among
   reassignments (C2 of B-17), and give placebo-dependent verdicts only as frequencies over redraws.
6. **No exploratory results in the abstract and conclusions** (background 4). Keep them in the Results with their
   qualifications (scoring shifted within GFM areas recalibrates each AOI to its evaluation labels; it is a diagnostic,
   not an evaluation of a forecast; etc.). The conclusions give the confirmatory record, the design's detectability, what
   the test does not answer and the lesson.
7. **GFM**: remove "independent" everywhere (abstract, heading of the data statement). For each result, state the labels
   the models were trained on (EMS-trained: baseline gain 1.4 × 10⁻² in Table S27, etc.; GFM-trained: D < 0). State that
   GFM's seasonal reference water is not assessed, and give the results of V1 (seasonal water assessable; v50).
8. **Estimates near settlements**: report in Sect. 3.6 the differences by sensor (optical 0.57, radar 1.19, ratio 0.50),
   product generation, hydrography layer and delay. State that the estimator of Eq. (2) was checked for calibration only
   above one. Remove the label permutation (cells permuted within AOIs) from the text; the name reassignment is the main
   inference. Remove "flood labels were also scarcer near water-related names" on GFM from the conclusions.
9. **Readings**: correct the description of the blocks (fa, hi, ng and others are composite lists; each name is matched
   against every language in a list). State in the Methods the assignment by country, matching on ASCII without
   diacritics and the different romanisation conventions (Cheshmeh against chashma, etc.). Keep dry-point siting only as
   one hypothesis beside the evidence against it (water-related settlements stand 2.1 m lower in site HAND; the deficit
   depends on the sensor), or drop it. For the qualifier pairs, cite both elevation (324 pairs, 85 %, archived table) and
   site HAND (Supplement, no difference), and do not claim position relative to water.
10. **Numbers and references** (editor M8; all minor points of referees 2 and 3): hedged derivations with count and rates
    of the same coder; the GFM n-grams from Table S44; word boundaries in the EMS n-grams; the κ table; one definition of
    the share of read names (9.6 %, 10.0 %, 10 %); numbers in no table (composition 0.037 and 0.430 with intervals, R with
    the tuned penalty, distance-band odds ratios, share of the largest uncapped AOI, qualifier pairs) into Table S47 or an
    explicit reference to the archived table; the definitions behind 150 and 154 in Table S47; one interval type for the
    share kept (Table 3 and Fig. 4b; B1); the sum of the ab words; cross-references (1.83 to 2.89 in Table S15; 1.0 %,
    13 % and 5.6 % in Table S45 and Table 6; n-grams in Table S44; R of read names in Table S43).
11. **What the test does not answer**: minor names, early forms, indigenous names and knowledge held in communities, in
    the abstract or conclusions and in the Discussion.
12. **Introduction and references**: discuss Mulyadi et al. (2025) in the Introduction. Describe the earlier test in Japan
    (three events, nine basins: source of names, how they were read, result) in Text S1. Limit the novelty claim to the
    search by a language-model agent (468 J-STAGE records, the first 20 read). Flag translated titles consistently; (Ed.)
    and (Eds.); version identifiers of models; source of the gloss of wæsse.
13. **Layout**: the overflow of Supplement Table S4 (amendments) and the cut-off row of amendment 17; Table S30 at the foot
    of the page; the caption of Fig. S7 over the page number; panels (b) and (g) of Fig. 3; the asterisk of Fig. 5a;
    the arrow of Fig. 4b replaced by the interval; a shorter caption for Table 3; the ‡ row of Table 4 as a footnote;
    every model label and a plain description of the verdicts in Table 1.
14. **Public repository**: correct the reference "paper, Appendix A" in PROVENANCE.

## B. Local computations (script v83; seed 20261107; run as utens jobs, not on the Mac)

1. **Share kept with Fieller intervals throughout**: for every record of Table 3 (GSW for events from 2022, earlier EMS
   floods, the JRC map, GSW for all events, EMS floods before or after), the 90 % Fieller interval of the share kept on
   all ground, built as in B9 of B-17 (predictions of the record runs, models fixed, 4,000 draws of location groups), so
   that Table 3, Fig. 4b and Table S45 give the same interval.
2. **Effective number of clusters for the odds ratios near settlements**: for the planned GFM estimate and the EMS
   estimates with and without the cap, the effective number of location groups contributing to the estimate of δ (the
   square of the sum of contributions over the sum of their squares), and the shares of the three largest groups.
3. **Table of reading variants** (referee 3, M5): for the existing variants (strict, floodplain, cleaned, edited, without
   wells and -ab), what changed, who decided and when (after which results), number of entries and of matched settlements,
   primary R and the odds ratio within 500 m, collected from existing outputs (no refit).

## C. Cluster computations

None (no refits). The inversion of the randomization test under injected effects (referee 2, M7 and M10), a jackknife
over super-clusters with every model refitted (M3) and settlement-level outcomes (referee 1, M1) are outside this paper
or too large, and the text says so.

## D. Decisions for the authors (not made in this amendment)

1. **Policy on generative AI** (precondition of review; decisive for the editor): the authors rewrite the text
   themselves or choose a journal that allows declared drafting.
2. **Toponymic audit** (referee 3's condition): the audit of the v69 guide (300 names), with rates by block and form
   carried into the corrected margins and bounds.
3. Checking the English translations of the protocol and every amendment.
4. Verifiable evidence of the freeze (e.g. a third-party time-stamped deposit of the development history; the
   development repository is not to be shared).
5. Deposit of the derived data (before review, at least a reserved DOI) and a DOI for the code.
6. Deposit of the conversation and prompt that compiled the lexicon.
7. Label QA (a visual check of source polygons against gridded labels in about 30 random AOIs, with an error rate).
8. Treatment of machine readings of indigenous names (CARE, inclusion statement).
9. Format and journal (including the policy on generative AI).
10. A test at the settlement scale (observed flooding around settlement footprints over many events), as a separate study
    registered in advance by the authors.

## Reporting rules

- Everything is post hoc (exploratory) and the confirmatory verdict is unchanged; it is reported as a record with its
  operating characteristics.
- The abstract and conclusions carry no exploratory result in either direction. They give the confirmatory verdict, the
  design's detectability, that real names do not stand out under randomization inference, what the test does not
  answer and the lesson.
- State that no upper limit has checked properties; give each construction's upper end as a description in a table, and
  do not convert upper ends into injected odds ratios.
- Rules, units and numbers of draws are not changed after computing. Implementation errors are corrected in addenda.

## Addendum 1 (implementation of B; 2026-11-02T01:17Z, commit 8bf2508; before the quantities of B1 and B2 were computed)

1. **Where B3 is built**: B3 (the table of reading variants) involves no computation and only collects existing
   outputs, so it is built in the paper's make_figures, not in v83 (strict reading: v52; floodplain reading and cleaned
   lexicon: the sources of the existing macros; edited lexicon: C2 of v71; without wells and -ab: C6 of v81). R and the
   odds ratio within 500 m are given only where they exist; missing values are marked as not computed.
2. **How B1 is built**: the record runs and pairs of arms are the five of B6 of B-14 (RECORDS of v73), without the
   versions excluding permanent water. All records use the same location-group weights (4,000 draws, seed 20261107). The
   interval is the 90 % Fieller interval from the covariance of the resampled (Δ through settlements, Δ per cell).
3. **How B2 is built**: the Eq. (2) estimator (AOI intercepts, distance bands, water-related name × band, logit of M1)
   is fitted with the AOI intercepts profiled out, and the influence of δ per location group is ψ_g = (S⁻¹ u_g)_j. The
   effective number is the inverse of the sum of squared variance shares. Tables: the planned GFM and EMS same-cell
   estimates (v41 table, frozen reading, GFM-assessable cells), the capped EMS estimate (v41 table, every analysed cell)
   and the uncapped EMS estimate (as B2 of v68, corrected reading). Computed in a utens job (script v83, output
   docs/verification/v83_b18.json).
