# B-19 Reformatting for Natural Hazards: in response to the fifteenth round of first-time critiques (editor, referee 1), the title and abstract return to the test and its limits, and the labels near settlements and the contribution of each activation are described

> English translation of `docs/protocol/amendments/B-19_natural_hazards_blind15.md`, drafted by a language model (Claude).
> The Japanese file governs. The authors have not yet checked this translation.

Status: **plan (2026-11-06T03:35Z, fixed before the quantities of B below were computed; the commit that added this file is
the record of the freeze)**. The confirmatory verdict (equivalent-null) does not change. Everything is post hoc
(exploratory).
Protocol: [`preregistration_v2.md`](../preregistration_v2.md) (FROZEN 2026-10-25T20:58:37Z, commit e625924).

Values seen when this was written: the confirmatory analysis and every result of B-2 to B-18, and the values that the two
reports of the fifteenth round read from the tables. The quantities of B below (the contribution of each activation, the
ceiling at 2 km and by settlement density, the flooding lost by binarisation, R and odds ratios on observable ground, the
distance-band odds ratios by how the maps were made) had not been computed. The ceilings at 250 m, 500 m and 1 km were
computed in B-6.

## Background

The fifteenth round of first-time critiques (NHESS version, 34-page paper and 44-page Supplement, development repository
108536c, public repository db90229; reports in [reviews/](../../research/reviews/), `2026-11-02_paper2_blind15_*.md`): the
editor returned the manuscript before review, and referee 1 (flood mapping) recommended major revision. The statistics
referee (2) and the toponymy referee (3) were stopped on 2026-11-06, when the user asked that no further critiques be run
with sub-agents; there are no reports from them.

The user's decisions (2026-11-06): the target journal is Natural Hazards (Springer). There is no research funding. The
co-author (P.-J. Chun) has approved. The DOIs of the data and code will come from Zenodo (the user does the browser steps).
The user reviews the manuscript.

Points common to the two reports:

1. **Title and main claim** (editor M1; referee 1, M5 and M9): the title's "nearest-settlement allocation dilutes
   settlement-scale information" rests on post hoc analyses that became the main subject after the confirmatory result,
   for one construction (250 m, nearest settlement, 500 m summaries, one score, one baseline), with wide intervals. The title
   and abstract return to the test that was done and its limits (the editor's route (a); referee 1, (a)). Dilution is placed
   among the exploratory diagnostics as a quantified property of this design.
2. **Qualifiers in the abstract** (editor M2 and M3; referee 1, M4 and M8): what the margin represents (about half of the
   per-cell improvement from the JRC 100-year hazard map; an odds ratio of 2.33 with the base models fixed, about 3.9 for truly
   water-related names, more with refitted models); that the decision rule declares equivalence about one time in three when
   the true R is at the margin; that inference rests on about two effective clusters; that no upper limit is stated; the
   intervals of the shares kept; that the JRC map is a reanalysis that includes the event years; that the name-reassignment
   result holds on the EMS labels.
3. **Activations** (editor M4; referee 1, M7): give the split that keeps activations whole the same weight as the primary
   split; the contribution of each activation to Δ and D (a forest plot) and results without the main contributors; the flood
   type of the main activations (whether it was inferred).
4. **Labels near settlements** (referee 1, M1 and M2): repeat the main comparison and the distance-band odds ratios on
   observable ground; the flooding lost by binarisation by distance band; the deficit by how the maps were made; make the
   behaviour of the labels near settlements one of the main results.
5. **The ceiling** (referee 1, M5 (b), minor 14): "depends only on the layout" is wrong. The ceiling depends on the reference
   model's errors and on the labels, and it varies with radius and settlement density.
6. **One standard for exploratory results** (editor M6; referee 1, M8): the flexible models after the language model's
   cut-off are described by the distribution of 200 placebo redraws, not by the single placebo (frozen reading: real minus
   mean 0.102, one-sided p 0.095; with Norway read, inconclusive in all 20 redraws). The within-area GFM result is described
   neutrally as an unexplained discrepancy between the two label sets, and the abstract's "did not stand out" is restricted
   to the EMS labels.
7. **Numbers and cross-references** (editor's minor points 1 to 26; referee 1's minor points 1 to 28 that can be fixed in
   the text).
8. **Readability** (editor M10; referee 1, M9): fewer model codes in the text, with descriptive names; the Supplement as
   supporting material.
9. **Journal format** (Natural Hazards): the Springer Nature LaTeX template (sn-jnl); an abstract of 150 to 250 words; 4 to
   6 keywords; line numbers; at most three heading levels; author-year references (every DOI as https://doi.org/); Statements
   and Declarations after the text (funding, competing interests, author contributions, data and code availability); the use
   of language models recorded in the Methods (Springer's policy: language models are not authors, their use is documented in
   the Methods, and humans are accountable for the final text); supplementary files called "Online Resource".

## A. Rewriting (no computation)

1. **Title**: the test that was done and its limit. Draft: "Machine-read settlement names did not detectably improve 250 m
   predictions of mapped river flooding, in a test insensitive to effects near settlements". "Protocol-based" leaves the title
   (a third party cannot verify the freeze); it stays in the text.
2. **Abstract** (150 to 250 words): question → data → the confirmatory test (R with its 90 % interval, the verdict as the
   record of a rule, what the margin means in terms of the JRC map and the range of odds ratios, the rule's operating
   characteristics, about two effective clusters, no upper limit) → what the design can detect (stated as post hoc
   diagnostics: the ceiling within 500 m, the shares kept with intervals, the JRC map allocated through settlements also
   equivalent) → the deficit of flood labels near settlements (post hoc) → the name-reassignment result on the EMS labels →
   what the test does not bear on → the lesson.
3. **Structure of the Results**: 3.1 sample → 3.2 the confirmatory test (the primary split beside the split that keeps
   activations whole, contributions of activations, the rule's operating characteristics, what the margin means) → 3.3 what
   the design can detect (exploratory diagnostics: ceiling and radius, records of known value, injected effects) → 3.4 labels
   near settlements (observable ground, binarisation, how the maps were made) → 3.5 evidence within areas after the result
   (name reassignment, the GFM discrepancy, few clusters, no upper limit) → 3.6 what the name readings capture.
4. **Ceiling**: correct "depends only on the layout"; state that it depends on the reference model's errors and on the
   labels, and how it varies with radius and density (B2).
5. **Margin**: "about half of the per-cell improvement from the JRC map (R = 0.095, 90 % 0.021 to 0.162)" and the range of
   odds ratios (2.33; 3.93 for truly water-related names; with refitted models an injected odds ratio of 2 gave a mean R of
   0.025), in the abstract, results and conclusions. State that the verdict says nothing about effects of plausible size.
6. **The JRC map**: wherever it is invoked, give its per-cell gain and interval, that it is a reanalysis that includes the
   event years, and that it covers basins of at least 500 km², at 90 m, without defences.
7. **After the language model's cut-off**: describe by the distribution of 200 redraws (values of Tables S41 and S40).
8. **GFM**: not called a "check". A second label set with a different definition and domain, stating first what it cannot
   observe. State that comparisons trained on GFM and the GFM odds ratio near settlements cannot be used to decide anything.
   The within-area result after shifting is an unexplained discrepancy, with two explanations side by side: strings (region,
   settlement type) and GFM's domain.
9. **Numbers and references**: the editor's minor points (the reference to Table A1 versus Table 3; 99 % is the TOST rate and
   the decision table's split is 91, 7, 1 and 1; numbers in the text that are not in the tables (98 % of Δ, 11 %, M1p's 51 %)
   moved into tables; the interval type of 0.91; naming the interval type of the shares kept; generation numbers as in
   Table S38; the inland rule and the inclusion criteria as implemented; 348 versus 344; "did not transfer" to "could not be
   shown to transfer"; 2.1 % versus 2 %; the reading in Table 6; a flow of the etymology sample counts; references to Fig. 1
   and Eq. (1); the word "secondary" in Fig. 4a; the lower end of JRC in Fig. 5b; the T0 − T4 row of Table 3; the numbering
   clash of the Supplement's tables; checking the Coates (2006) citation; citing Owen (2007); data citations; number format)
   and referee 1's minor points that can be fixed in the text (counts of map products, the treatment of traces, thresholds
   as fractions of the whole cell, hydrography layers, how GFM scenes are combined and the GFM versions, the positional error
   of rounded coordinates, bilinear resampling of HAND, Brier skill and within-AOI CSI next to R, how the ratio was estimated,
   multiplicity of strata, the "other" class of κ).
10. **Readability**: the only codes in the text are the position model (M0), baseline model (M1), reference model (M1c) and
    names model (M3n); the others get descriptive names. Table 1 keeps every code (the Supplement uses them).
11. **Language models in the Methods**: following Springer's policy, what they were used for and what the authors checked
    (the extent of the authors' checks is for the authors to state; D).
12. **Statements and Declarations**: no funding ("No funding was received for conducting this study."), no competing
    interests, author contributions (CRediT), data and code availability (Zenodo DOIs; placeholders until reserved).

## B. Computation (script v84; seed 20261108; run as a utens job, not on the Mac)

1. **Contribution of each activation**: with the held-out predictions of the primary split (frozen reading world_v2,
   corrected reading world_v2_lexfix) and of the split that keeps activations whole (frozen world_v2_actgrp, corrected
   world_v2_b7actfix), for each activation a the contribution to R, r_a = Σ_{i∈a}[(p_M1c − y)² − (p_M3n − y)²] /
   Σ_i[(p_M0 − y)² − (p_M1 − y)²] (which sums to R over all a), and its share of D. For the 12 activations with the largest
   |r_a|, 90 % intervals from redrawing only that activation's location groups (all else fixed; 4,000 draws). R without each
   activation (predictions fixed) and without the top three. Shares of |r_a| held by the top 1, 3 and 10. Each activation's
   country, event date, number of AOIs and number of AOIs whose flood type was inferred.
2. **Ceiling by radius and density**: R_max(S) = Σ_{i∈S}(p_M1c − y)² / N / D (models fixed), S = cells whose nearest
   settlement lies within 250 m, 500 m, 1 km or 2 km and has a water-related name, or any settlement; split into the parts
   from flooded and dry cells. By tercile of AOI settlement density (the median distance from the AOI's analysed cells to
   their nearest settlement): the share of cells within 500 m of water-related names and the contribution to the ceiling.
   Corrected and frozen readings.
3. **Flooding lost by binarisation**: of the summed riverine flood cover on the 25 m sub-grid (fractions of the whole
   cell), the share in cells labelled dry, by distance band (within 500 m, 0.5 to 1 km, beyond 1 km), on the analysed cells
   (capped) and on all evaluable cells (the uncapped run). For each nearest settlement, the share whose cells within 500 m
   have a flood label, any riverine flooding, a cell with cover of at least 0.10, and summed cover of at least half a cell.
4. **Observable ground**: cells with an ESA WorldCover (2021, 10 m) built-up share below 0.10 and a tree share below 0.50.
   On that ground, R of the confirmatory comparison (models fixed; with that ground's D and the overall D; 4,000 draws of
   location groups; 90 % Fieller; frozen and corrected readings), the distance-band odds ratios (AOI fixed effects, logit of
   M1, clustered by location group) and the odds ratio of Eq. (2) for water-related names within 500 m (frozen reading),
   beside the same on all ground.
5. **How the maps were made**: the distance-band odds ratios within strata of sensor (optical, radar, unknown), product
   generation, presence of a hydrography layer, tercile of the delay to the image (cut points made as in B6 of B-16) and
   whether the AOI's labels contain flood traces.

## C. Not done (outside the scope, or too large a computation)

Comparing allocation schemes on the records (distance-decay weights, settlement footprints, evaluation at settlement level),
another resolution, finer support on EMS detail maps, same-day single-scene GFM labels, OpenStreetMap embankments, the HAND
accumulation threshold, site contrasts on bare-earth heights, decomposing the strong baseline, and characterising the 452
excluded AOIs. The text says what was not done and what is left to the settlement-scale test (a separate, pre-registered
study).

## D. For the authors (not done in this amendment)

1. The authors check the English translations of the protocol and of every amendment (Natural Hazards reviewers cannot read
   the Japanese documents).
2. The description of language-model use in the Methods: what the authors themselves wrote and checked (Springer's policy
   requires it to be documented, and the authors are accountable for the final text).
3. A toponymic check (the v69 guide and 300-name sample).
4. Zenodo DOIs (the user creates the token and stores it in the Keychain; Claude makes the deposit after confirmation).
5. ORCID, cover letter, suggested reviewers.
6. Third-party evidence of the time of the freeze (the development repository is not to be shared).
7. Publication of machine readings of indigenous names (CARE).

## Reporting rules

- Everything is post hoc (exploratory) and the confirmatory verdict does not change. The confirmatory verdict is reported as
  a record, with the rule's operating characteristics, the effective number of clusters and what the margin means.
- The abstract gives the confirmatory verdict with its qualifiers, what the design can detect (stated as post hoc, with
  intervals), the deficit of flood labels near settlements (post hoc), name reassignment on the EMS labels, and what the test
  does not bear on. Exploratory results that favour or disfavour names are treated by one standard whichever way they point,
  and are not in the abstract.
- No upper limit is stated. The quantities of B are descriptions, not tests (intervals are given, no verdicts).
- Rules, units, counts and thresholds (built-up 0.10, trees 0.50, radii, terciles) are not changed after the computation.
  Implementation errors are corrected in addenda.

## Addendum 1 (implementation of B; 2026-11-06T03:54Z, before the intervals of B1 for the split that keeps activations whole were seen)

1. **Input on the NAS**: the first job stopped at B5 (the EMS inventory data/global/ems/inventory.json was not on the NAS, so the
   image times could not be read). The local inventory (sha256 f5c55eff…8699, retrieved on 2026-10-25, before the freeze) was
   uploaded and the job resubmitted. No rule was changed.
2. **The unit of the B1 intervals**: in the runs that keep activations whole (world_v2_actgrp, world_v2_b7actfix), the unit
   stored with the predictions is the activation group, not the location group, so the first output gave no interval for the
   leading activations of those two runs (each counted as "one group"). As planned, the intervals redraw location groups: every
   run now maps each AOI to its location group of the primary split (world_v2), checking that this equals the stored unit in
   the primary runs. Point estimates, R without each activation and the shares do not change. Only B1 is rerun.
