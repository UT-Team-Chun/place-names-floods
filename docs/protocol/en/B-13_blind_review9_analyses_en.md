# B-13 Restructuring and exploratory analyses in response to the ninth round of first-time critiques (NHESS version, 2026-10-29)

English translation, drafted with a language model; the Japanese file `docs/protocol/amendments/B-13_blind_review9_analyses.md`
is the record and governs where the two differ.

Status: **plan (2026-10-29T12:19Z, commit c096169, fixed before any quantity in B and C below was computed)**. The
confirmatory verdict (equivalent-null) is not changed. The protocol's verdict is reported unchanged as the record, and a
post hoc test that corrects its error-rate problem (B1) is reported beside it. Everything below is post hoc (exploratory).
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924).

Values seen when this document was written: the confirmatory analysis and the results of B-2 to B-12, except items C1 to C4
of B-9 (being computed on the cluster), item B11 of B-12 (penalty on the name coefficients only, being computed locally) and
item C1 of B-12 (injection within 1 km, being computed on the cluster). None of the quantities in B and C below has been
computed.

## Background

The ninth round of first-time critiques (NHESS version, 33 pages of main text and 32 of supplement, development repository
at 40b0a3d, public repository at 5622df7; reports in `docs/research/reviews/2026-10-29_paper2_blind9_*.md`) came back with
the editor returning the manuscript before review and three referees recommending major revision. Points common to the four
reports:

1. **The decision rule does not control its error rate**: the rule that divides the fixed-model interval by the point
   yardstick declares equivalence in 32–33 % of samples when the true R lies at the margin (B-12, item B1). Refitting makes
   the standard deviation of R about 2.3 times larger. A valid equivalence test that carries the covariance of numerator
   and denominator, refits the second stage and resamples at the higher units (a Fieller-type TOST against the ratio
   margin) should be the main basis of interpretation (referee 2, M1 and M4; editor, M3).
2. **Numbers attached to the wrong quantities**: the odds ratio 1.83 in the Discussion belongs to the refitted
   location-group interval (upper limit 0.013). For the two-way interval (upper limit 0.036) it is about 2.1 on the nominal
   scale (about 2.6 on the estimator scale). The 18 % p0 is the value of coder 1, the model that compiled the lexicon; for
   coder 2 it is about 23 %. The dilution example uses coder 1's rates (referee 2, M6; referee 3, M5; editor, M6).
3. **Results in the positive direction under-reported**: the total increment of H4 (positive-small) and the within-area R of
   0.011 on GFM labels (above the placebo) are missing from the abstract and conclusions (editor, M5; referee 2, M9).
4. **Leakage in the controls**: surface water (1984–2021) and the recurrence record include the event itself and later
   events. Show the dilution with controls known before the event (referee 1, M1 and M7; editor, M1).
5. **Near-settlement estimates**: much of the text goes to estimates that do not identify an effect; the uncapped estimate
   was selected; intervals rest on few clusters; a corroboration rule; label weakness is put into the conclusions (referee 2,
   M7; referee 1, M2, M3 and M6; referee 3, M7).
6. **Lexicon**: homographs (-ach and -bach, -ore and Store, mont-, Ried, kacha, taung), blocks attached to countries (all of
   Italy read with the German block, Belgian -eau read as Au), matching rules that differ between languages, relative
   qualifiers, and a check that is only agreement between coders (referee 3, M2 to M6; editor, M2).
7. **Length and coined terms**: many uses of water-side, not-meaning and M symbols; a long title (34 words). Four parts for
   the main text (the confirmatory result with all its uncertainty, a few pre-specified checks, the design ceiling and
   dilution, validity of the readings).
8. **What cannot be checked**: the freeze record, the unchecked English translations, the data deposit, the policy on
   generative AI (for the authors).

## A. Rewriting (no computation; corrections in the paper and supplement)

1. **Title** (about 20 words; "protocol-based", because third parties cannot check the advance fixing): "Water-related
   terms in settlement names, assigned to the nearest grid cells, did not improve within-area predictions of mapped river
   flood extent: a protocol-based test in 500 Copernicus EMS areas".
2. **Abstract**: the limits of the design before the numbers (an effect confined to the surroundings of settlements cannot
   reach the margin; assignment through settlements also dilutes information known before the event). Define R, activation
   and the baseline gain in one sentence. Describe the outcome labels (a single mapping; areas chosen for extensive
   flooding). Report the protocol's verdict as the record and give the result of the valid test of B1. State that the
   two-way interval does not lie inside the margin and that the reading as a bound was chosen after the result. Report the
   total increment of H4 and the within-area R on GFM labels. Either state that the odds-ratio translation refers to "an
   effect spread uniformly over all ground nearest a name" or drop it. Write "the Brier score increased slightly".
3. **Sect. 3**: the protocol's verdict (record) beside the valid test of B1; the split that groups activations, with both
   codes, next to the main split. The total increment of H4 (positive-small) and the flexible H1 (positive) in the text and
   in Table C1. State that "H1 failed" refers to the prediction about the nested form.
4. **Sect. 4**: gather the main uncertainty into B1, with the two-way and refitted intervals briefly as its components.
   Alternative denominators (Y3 and ΔBS, the increment in Brier skill), the influence of South Asia (B6), leakage in the
   reference model (B5). For the relaxed second stages, give the verdict of each version (the continent interaction is
   negative; correct the inconsistent sentence on JRC).
5. **Sect. 5**: the leak-free controls (C1) as the main ones; keep the controls that include the event as "what the design
   can pick up". Write that bounds "correspond to" injected effects and tie each bound to its interval (B2). Give the ranges
   for each no-effect yardstick with equal weight. State that Monte Carlo ranges are not confidence intervals. Give the
   dilution example with the rates of both coders (B3). Rewrite "past EMS floods were picked up in the AOIs with a record"
   with numbers and a verdict.
6. **Sect. 6**: (a) Readings: describe Table 5 as "agreement with language-model codes". Correct the coder attribution.
   The table of systematic lexicon errors (B10) goes to the Supplement, with the main points in the text. Describe N-read as
   a reading that includes knowledge of places. (b) Near settlements: one paragraph. The planned GFM estimate with its
   planned interval; the uncapped EMS estimate as an exploratory result ("no positive association is visible on the mapped
   labels"), with the statement that it does not identify an effect. The account by class goes to the Supplement. Drop the
   corroboration rule. Write label weakness as an untested hypothesis. Treat the GFM estimates the same way in both
   directions.
7. **Discussion and conclusions**: practical → simple. What NHESS readers gain (the general lesson that assignment through
   settlements also dilutes information known before the event; what a settlement-scale test needs: GFM does not evaluate
   built-up areas, so coherence and double-bounce methods, and the Global Flood Database).
8. **Terms**: water-side → water-related (state in the definition that it includes crossings, valleys and "lower"
   qualifiers). Replace not-meaning by "not distinguishable from the placebo". Turn Table 1 into a glossary (M0 to M3s, T0 to
   T4, Y0 to Y4). Define activation, location group and ground at first use.
9. **Supplement and appendix corrections**: amendment 12 and the eighth and ninth critiques in the timeline. Correct the
   section numbers in Table S5 (status of the plan) and add rows for amendments 12 and 13. The inconsistent title and the
   missing note b of Table S6. Countries in the activation rows of B2 in Table S17. The garbled characters in the row of B4.
   How the numbers of "AOIs other than Sentinel-1" (68, 111, 93, 108) were counted, in the caption. The caption and headings
   of Table 3. "T4 = M1 + HAND" in Fig. 1b. The legend for symbols and colours in Fig. 4. Plain row names in Fig. S6. Numbers
   not in any table (median distance 1.5 km, the highland odds ratio, calibration intercept and slope, the near-settlement R
   of 0.060, the values of H4, the correlation 0.04, the 9.7 % for Australia, the 1.45 for floodplains, the calibration
   gain, the quantities of Appendix C) into tables. References (ESA WorldCover, Natural Earth, the GSW version, literature on
   spatial cross-validation and flood susceptibility mapping, Owen and Nadeau–Bengio in the text, uncited references in the
   Supplement). sha256 in full. A summary of the earlier test in Japan in the Supplement.
10. **Corrections of numbers and wording**: 33 % → 32–33 %; "more than 98 %" → "at least 98 %"; the wording of one-sided p;
    the inconsistent sentence on 0.017; land cover for events from 2022 is the 2020 map; the gazetteer was retrieved after
    all events; the breakdown of the 154 typed AOIs; the scale of "two-thirds" in Appendix C; the 641 estimates come from
    amendments 2 to 10; 0.5 % for Africa against 1 % in Table 5; the identical flexible upper limits of T0 and T1 in Table 2;
    the identical 90 % and 95 % values for JRC in Table S13 (all checks and corrections of existing values).
11. **Limits of reproduction**: simulations differ by up to 1.4 Monte Carlo SE between arm64 and amd64; the rebuilt GFM
    labels agree for 189 of 192; in the text or the Supplement.

## B. Computation from existing products (local, script v70, random seed 20261102)

1. **A valid equivalence test** (Fieller-type TOST): compute L± = ΔBS ∓ 0.05·Δ_base (ΔBS = BS(M1c) − BS(M3n), Δ_base =
   BS(M0) − BS(M1)) on the same resamples; equivalence if the one-sided 95 % upper limit of L+ is below 0 and the
   one-sided 95 % lower limit of L− is above 0 (the 90 % Fieller interval of R lies inside ±0.05). Resampling: (a) location
   groups, models fixed (4,000 replicates); (b) location groups, the second stage refitted in each replicate (500, the same
   refit as item B7 of B-11); (c) activations as the unit (two-stage resampling of activations within countries, second
   stage refitted, 500); (d) (c) plus the variance from refitting the base models, from the 13 fold assignments (v47) with
   the Nadeau–Bengio correction, as a normal approximation. For the main comparison (the confirmatory run and the corrected
   run), the nested and total increments of H4, and the country split. Operating characteristics: with the double
   resampling of item B1 of B-12 (models fixed), the share of samples in which the L rule declares equivalence when the true
   R is set to ±0.05 and to 0 (uniform and proportional shifts).
2. **Translation of the bounds**: the injected odds ratio (nominal and on the estimator scale) corresponding to the upper
   limit of each interval (fixed 0.003, refitted, two-way 0.036, each upper limit of B1), from the curves of v64; for each of
   the four no-effect yardsticks; with Monte Carlo ranges.
3. **p0 of coder 2**: the share of water-related derivations in the stratum of names not read as water-related, by coder 2's
   codes (with an interval, counted directly). The dilution table and the misreading correction (the probabilistic bias
   analysis of ten_misclass) with the rates of both coders.
4. **Strata by size near settlements**: the size terciles of v55 (the table of bounds) as numbers in the text (no new
   computation; a table only).
5. **Leakage in the reference model**: on the main split (second stage of the corrected code), R between M1c without the AOI
   mean of the character-model score and M3n.
6. **Influence**: for the main R (the confirmatory run) and for (b) of B1, the values with each of the 10 most influential
   location groups removed in turn, and the effective number of clusters based on influence (relative to Σ influence²).
7. **Alternative denominators and absolute quantities**: for the main comparison and the nested and total increments of H4,
   R with Y3 (the gain over a constant prediction), ΔBS (in Brier units), and the increment in Brier skill relative to
   climatology.
8. **Joint test of the name coefficients**: fit the second-stage M3n on all cells and test jointly that all 14 name
   coefficients are 0 (Wald, cluster-robust by location group). Find and report the "held-out coefficient test" of addendum 3
   of B-7 (heldout_coef_test and wald_tests in v55).
9. **Checks near settlements**: (a) the uncapped EMS odds ratio (pairs) with settlements whose coordinates are not on whole
   arc-minutes only; (b) the uncapped EMS odds ratio (point) with each activation removed in turn, and the pairs interval
   without the top two (EMSR052, EMSR130); (c) the odds ratio without products of generations 0 and 1; (d) p for the planned
   GFM estimate and the uncapped EMS estimate against the distribution under permutation of names within AOIs (200);
   (e) the 2×2 agreement of EMS and GFM by distance band and built-up share; (f) a test of the interaction water-related ×
   500 m × sensor (optical and radar) (the table of bounds); (g) the uncapped EMS odds ratio and the planned GFM odds ratio
   without settlements whose GeoNames feature code is PPLX (section of a populated place), PPLH, PPLQ or PPLW (historical,
   abandoned, destroyed).
10. **Lexicon checks**: counts of homographs and language assignments (-ach names that are not -bach names, -ore and -ora
    with Store, monastery names matched by mont-, Ried, kacha, taung, Belgian -eau read as Au, names read by the German block
    in Italy, how whole-word and suffix matches of ab are counted in Iran, ala-/yli-). A table of blocks and countries, and a
    table by language of match positions (whole word, start, end) (in the Supplement). The numbers of matches gained and lost,
    by language, when the local GeoNames name column (Unicode) is read with the rule ä→a, ö→o, ü→u, å→a, ø→o, æ→ae, ß→ss
    and combining marks removed, compared with asciiname.
11. **N-read**: shares of classes and of unknown by block and by country.
12. **The two coders**: agreement (κ), and a file listing GeoNames ID, country, stratum and both coders' codes (public,
    separate from the audit sample).
13. **AOI-balanced metrics**: medians over AOIs of ΔAUC and ΔPR-AUC (M1c → M3n), and the increment in Brier skill with equal
    weight per AOI.
14. **Diagnosis of the collapse when refitting on GFM** (as far as local products allow): flooded shares of GFM and EMS by
    AOI, and the baseline gain by AOI on GFM.
15. **Pairs of relative qualifiers** (referee 3, M7; Supplement): pairs of settlements in the same country whose names are
    the same, with three or more letters, once the qualifiers of C2 (d) are removed (for languages where the qualifier
    attaches at the start, once that start is removed), where one carries a qualifier from the lower side and the other one
    from the upper side, within 10 km of each other, both with coordinates not on whole arc-minutes. For each pair:
    (a) the difference in HAND at the sites (lower − upper); (b) the difference in the EMS flooded share of the evaluable
    cells within 500 m (when both have such cells); (c) the same within the domain of the planned GFM analysis. The number of
    pairs, the median difference, the share of pairs in which the lower member is lower (floods more) with the p of a sign
    test, and a 90 % interval resampling pairs by location group.

## C. Computation on the cluster

1. **Leak-free controls** (entered in the second stage, not in the base models): (a) surface-water occurrence (GSW
   1984–2021) only in the AOIs of events from 2022 onwards (directly per cell, and summarised within 500 m of settlements and
   assigned through settlements); (b) the recurrence record of other EMS events, from earlier events only (per cell and
   through settlements); (c) the JRC global river flood hazard map (inundation at a 100-year return period) summarised within
   500 m of settlements and assigned through settlements. Main split, corrected code. With the verdict from the verdict
   table.
2. **A lexicon with its documented errors corrected** (a derived file `docs/verification/v70_lexicon_corrected.yaml`,
   written mechanically from the frozen lexicon by `backend/verification/r_okauchi/v70_lexicon_corrected.py`; the frozen
   file is not changed). Only the entries named by referee 3 and the three classes the referee asked to separate are
   corrected; whether other entries are right is left to the toponymist audit (D5). The matching rules (split at
   non-letters, at least two letters after an affix), the six classes and the 14 name features are unchanged (removed forms
   are not read). Forms are removed in all positions (word, start, end, inside), except -ach and mont- in (b), which are
   removed only in the position stated.
   - (a) The cleaning of v57 (remove start and end entries of two letters or fewer, and the English endings -burn, -bourne,
     -borne, -pool, -well, -wells).
   - (b) Homographs (referee 3, M3 (iii) and (ii)): the de ending ach (the ending bach stays, so -ach is read only as -bach);
     de ried; the fr start mont (the words mont and monts stay); no ore and ora; fa kacha; my taung; pt and es monte; en bank
     and banks.
   - (c) Riverine entries named as having no flood meaning in toponymy (M3 (ii)): sv ang, angen, angar, vall; no eng, enga,
     engene, voll, vollen; da eng, enge (the same word as no eng); no oy, oya; fr ile, iles, isle, ilot; sv and no sand,
     sanden. Island words stay in the other languages (they are also used for islands in rivers, and the analysed areas are
     river floods).
   - (d) Relative qualifiers (words that mark one member of an upper/lower pair; M3 (ii), M7). Lower side: en low, lower,
     nether; de unter, nieder; nl laag, neder, beneden; fr bas, basse, sous; it basso, bassa, bassi, sotto, inferiore; es
     bajo, baja, bajos, bajas, abajo; pt baixa, baixo; el kato; sh donji, donja, donje; sl spodnji, spodnja, spodnje; bg
     dolna, dolni, dolno, dolen; sq poshte, poshtme; ro jos; hu also; cs dolni, dolny, dolne, dolna; pl dolny, dolna, dolne;
     uk nyzhnii, nyzhnia, nyzhnie, nizhniy, nizhnyaya, nizhneye, nizhne, nizhni; lt zemutinis, zemutine, zemieji; et alam;
     fi ala; sv, no and da nedre; tr asagi; ar tahta, tahtani, tahtaniya, sufla, sufli; fa nichla, nichli, nichle, payan,
     payin, paeen, pain, payeen; hi nichla; id bawah, hilir, ilir; tl ibaba, baba; sw chini; mg ambany; ka kvemo. Upper side:
     en high, upper; de ober; nl hoog, boven; fr haut, haute, hauts; it alto, alta, alti, sopra, superiore; es alto, alta,
     altos, altas, arriba; pt alto, cima; el ano, epano; sh gornji, gornja, gornje; sl gorenji, gorenja, gorenje, zgornji,
     zgornja, zgornje, gornji; bg gorna, gorni, gorno, goren; sq siper, siperme; ro sus; hu felso; cs horni, horny, horne,
     horna; pl gorny, gorna, gorne; uk verkhnii, verkhnia, verkhnie, verkhniy, verkhnyaya, verkhneye, verkhne; lt
     aukstutinis, aukstoji; et ulem; fi yla; sv, no and da ovre; tr yukari; ar fawqani, foqani, fouqani, fawqa, oulia, ulya,
     ulia; fa bala; hi upar; id atas, hulu, ulu; tl itaas; sw juu; mg ambony; ka zemo. Adjectives of height (de hoch and
     tief, vysok- in cs, pl and uk, and so on) and nouns for valleys and lowlands stay.
   - (e) Springs and wells (M3 (ii); the referee asked for them to be taken out of the water-related composite): en well,
     wells, spring, springs, soak; de born, brunn, brunnen, quell, quelle; fr fontaine, fontaines, source, sources, puits; it
     fonte, fontana, fontane, sorgente, pozzo, pozzi; es fuente, fuentes, pozo, pozos, manantial, ojo; pt fonte, fontes, poco,
     nascente, cacimba; el pigi, pigadi, pigadia, vrysi, vrisi, vryses, vrises, krini, kefalovryso, vrys, vris; sh vrelo,
     vrela, izvor, izvori, bunar, studenac, vrutak; sl studenec, izvir, studenci; bg izvor, izvori, kladenets, kladenci,
     bunar, cheshma; sq burim, burimi, pus, puse; ro izvor, izvorul, izvoare, fantana, fantanele; hu kut, forras; cs studnice,
     studna, studnicka, pramen, prameny; pl zrodlo, zrodla, studnia, studzianka; uk krynytsia, krinitsa, kolodez, kolodiazi,
     dzherelo, istok; lt saltinis, saltiniai; lv avots, avoti; et allikas; fi lahde; sv kalla, kallan; no and da kilde; tr
     pinar, kuyu, bulak, kaynak; ar ain, ayn, aioun, ayoun, bir, biar, abar, hassi; fa chashma, chashmeh, chah, chahi; vi
     gieng; id sumur; tl bukal; sw chemchemi, kisima; bantu chitsime; af fontein, put; mg loharano; qu puquio, pukyu; nah
     ameyal; ka tskaro, tskarosi; kk bulak, bulaq; so ceel; am minch, minchi.
   - (f) Rapids and waterfalls (M3 (ii)): en falls; fi koski; sv fors, forsen; no foss, fossen.
   - (g) Add the named omissions (M3 (iv)): the fr word sur (riverine); the cs and pl word nad (riverine); the de words and
     endings hamm, damm, deich (riverine); the es word rambla (water); the hu word and ending fok (water); the km word
     trapeang (water).
   - (h) Countries and blocks (M4 (i)): remove the de block from Belgium and Italy (keep it for Germany, Austria,
     Switzerland, Liechtenstein and Luxembourg only).

   With this lexicon: the main comparison (world_v2_b13lexfix, arm M0, stacks M1c, M3n, M3np) and the near-settlement odds
   ratios (capped and uncapped, EMS and the planned GFM). The differences in the numbers of matches from the frozen lexicon
   and from the lexicon of v57, by language. The rules are as written in this plan and are not changed after the results
   are seen.
3. **Early stopping by location group**: choose the number of iterations of M0 and M1 by an inner validation that splits the
   training side by location groups (instead of early stopping on a random 10 %), and refit the main split. R, the
   denominator, and the comparison with the constant prediction.
4. **Injection within 1 km** (B-12, C1): enter the results of the run in progress.

## D. Decisions for the authors and a separate study (not done in this amendment)

1. The policy on generative AI (B-12, D1): where it is stated, which parts a language model drafted and what the authors
   checked.
2. Checking the English translations (the protocol and all amendments).
3. A verifiable freeze record (development history and cluster records for third parties). The instruction is not to share
   the development repository, so the authors decide.
4. Deposit of the derived data (Zenodo, in a form referees can use).
5. The toponymist audit (300 names and the guide). Once it is done, redo the main comparison and the near-settlement odds
   ratios with the audited readings (rules fixed before the audit).
6. Rebuilding the labels and the baseline (detailed maps, maximum extent from monitoring maps, image times, generations,
   scene-by-scene GFM rebuilding with advisory flags, FABDEM and HAND along the flow path, the JRC European maps,
   GFPLAIN250m, a proxy for defences) belongs to the settlement-scale study. In this paper, claims that depend on them (for
   example that labels are weak near settlements) are written as untested hypotheses.
7. The article type (Brief communication) or the choice of journal.
8. Rebuilding the set of settlements and the readings (redo the nearest assignment without PPLX, PPLH, PPLQ, PPLW and
   records that are personal names; read Unicode names and alternate names; judge the language and layer of each name;
   registers of minor names) belongs to the settlement-scale study.

## Reporting rules

- Everything is post hoc (exploratory), and the confirmatory verdict is not changed. The protocol's verdict is reported as
  the record, and the valid test of B1 is the main basis of interpretation (added to the deviations in Appendix A).
- Results are reported whatever their direction. If some form of B1 cannot declare equivalence, that is written as such.
- Rules, units and numbers of replicates are not changed after computation. Implementation errors are corrected in an
  addendum.

## Addendum 1 (record of the implementation; 2026-10-29T12:37Z, commit a7b0f33, before any quantity in B and C was computed)

1. **C1(a)**: `--min-event-year 2022` (a table of the AOIs of events from 2022 onwards only, with the split, base models
   and second stage rebuilt). Surface water through settlements uses the existing `settlement_gsw.parquet` (mean
   occurrence and seasonality within 500 m of the settlement point; B-6 item 3(b)) in M1cG → M3oG, and the mean without
   permanent water in M1cGn → M3oGn. Per cell, the existing `cell_gsw.parquet` (B-8 item C4) in M1cW → M3oW. The main
   comparison of names (M1c → M3n, placebo M3np) on the same table. Run world_v2_b13gsw22.
2. **C1(b)**: `--recurrence --recurrence-prior`. Other events are those whose extent overlaps, whose activation differs,
   and which are more than 30 days earlier than the event (item 4 of B-7 used both earlier and later events). Through
   settlements M1cE → M3oE, per cell M1cEc → M3oEc. Run world_v2_b13recprior.
3. **C1(c)**: `--jrc-set`. Three features assigned from the nearest settlement: the mean flooded share at a 100-year
   return period over the cells within 500 m of the settlement point (cells with at least half of the JRC tile pixels),
   whether any of those cells floods (the "presence" of the plan), and an indicator that a value exists (M1cJ → M3oJ). As
   a comparator, the same run also gives the per-cell flooded share at 100 years (M1cJc → M3oJc); this comparator is not
   in the plan and is added before the results.
4. **C2**: run world_v2_b13lexfix (`--lexicon docs/verification/v70_lexicon_corrected.yaml`, sha256 c1ed4f18300a…). The
   derived lexicon has 2,148 entries against 2,538 (an entry is a pair of form and position; removed: 15 by the cleaning
   of v57, 15 homographs, 37 riverine entries, 181 relative qualifiers, 145 springs and wells, 10 rapids and waterfalls;
   added: 13). Record in `docs/verification/v70_lexicon_corrected_record.json`.
5. **C3**: `--grouped-es`. In every fit of the base models (M0 and M1; all inner cross-fits and the fit on the whole
   training side), 10 % of the training location groups (rounded, at least one) are chosen with the fit's seed as an inner
   validation; the model is fitted on the rest for up to 200 iterations without early stopping, and the whole training
   side is refitted with the number of iterations that minimises the validation log loss (the previous default stopped
   when 10 iterations brought no improvement on a random 10 %). Run world_v2_b13gres.
6. All C runs use the main split (10 folds, one fold per job, 4 CPU and 14 Gi), the corrected code and the five default
   seeds.
7. **B1**: (a) from the held-out predictions of the confirmatory run (world_v2) and the corrected run (world_v2_lexfix),
   resampling location groups 4,000 times. (b) and (c) refit with the stage-2 design of the corrected code
   (world_v2_b8stage2_local), as item B7 of B-11 did (the stage-2 design of the confirmatory run was not saved). The
   two-stage resampling of (c) draws countries with equal probability as many times as there are countries, and for each
   drawn country redraws its activations as many times as it has activations. (d) computes L± from the held-out
   predictions of the 13 fold assignments of v47 (world_v2_lexfix and world_v2_fs1b, fs2b, fs3 to fs12), corrects their
   variance by Nadeau–Bengio and adds it to the variance of (c) (normal approximation, one-sided 95 %). L± are on the scale
   of the Brier score per cell (divided by the sum of weights). Fieller intervals are also given by normal approximation
   from the variances and covariance of the replicates of (a) to (c).
8. **B15**: HAND at the site as in v37 (mean hand_min of the cells within 177 m of the settlement point), and whether
   coordinates lie on whole arc-minutes by the rule of v37. Settlements are those nearest to at least one cell of the
   analysis table (the table of v41, capped). Pairs take lower and upper settlements of the same country and base name in
   order of distance, one at a time (each settlement once). Flooded shares use the cells of the same table within 500 m;
   GFM uses the cells of v41 with GFM labels. Intervals resample the location group of the lower settlement (2,000
   replicates).

## Addendum 2 (correction of the implementation of C1(a); 2026-10-29T12:40Z, commit 39f2082, before any quantity in C was computed)

The `--min-event-year 2022` of addendum 1, item 1, selects AOIs before the cells are drawn within the caps of the location
groups, so cells that are not in the main table enter, and the per-cell surface-water table (`cell_gsw.parquet`, only the
1,507,486 cells of the main table) has no value for them. It is therefore replaced by `--post-min-event-year 2022`: the main
table (cells drawn over all AOIs) is built first, the rows of AOIs whose event is before 2022 are then removed, and the
location groups and folds are rebuilt afterwards (events from 2022: 210 AOIs, 581,512 cells, 69 activations). Everything
else is as in addendum 1.

## Addendum 3 (record of the implementation of B1 to B5, B7 and B13; 2026-10-29T12:47Z, commit 03c120e, before any quantity in B was computed)

1. **B1(c)**: an activation that spans two or more countries is assigned to the country with most of its cells (for the
   nested resampling).
2. **Intervals of B1**: the one-sided limits are the 5 % and 95 % points of the replicates of L±. Fieller intervals use the
   normal approximation (90 %) with the variances and covariance of the replicates of ΔBS and Δ_base. Percentile
   intervals of the ratio are given beside them. The point estimates of (b) and (c) are the unweighted refit (the
   corrected run).
3. **B1(d)**: the variance of L± is the variance built from Var(ΔBS), Var(Δ_base) and their covariance in the replicates
   of (c), plus the Nadeau–Bengio variance of L± over the 13 assignments (n_test/n_train = 1/9). The Fieller interval uses
   the same summed variances.
4. **B2**: upper limits on the scale of R are multiplied by the corrected run's yardstick to reach the Brier scale; the
   bound is the odds ratio at which the mean gain reaches that value minus each no-effect reference ΔBS_null (the "bound"
   reading of v59). Shapes A and B. The upper limits of B1 are the Fieller upper limits of (a) for the corrected run, (b),
   (c) and (d). The estimator scale follows the rule of v68 B3 (the curves of B-10 B11).
5. **B3**: coded rows of coder 2; water-related means any of the five water-related classes (as in v57b). Wilson 95 %
   intervals. The probabilistic bias analysis uses the stored pairs replicates of v57 with rates drawn from
   Beta(k+1, n−k+1), for each coder (20,000 draws).
6. **B5**: the character model's AOI mean is removed from both the reference (M1c) and the names model (M3n) (second stage
   refitted on the design of B-8 C2). The standard R from the same refit is given for comparison.
7. **B7**: Y3 = Brier of the constant prediction minus Brier of M1. The constant is the training-side flooded share of
   each fold on the main split, and the test share on the single-fold llmcut split (B-12 B17). Increment in Brier skill =
   ΔBS divided by the Brier of that constant.
8. **B13**: for each AOI with both flooded and dry cells, the difference in AUC and PR-AUC between M3n and M1c. The
   increment in Brier skill uses the AOI's own flooded share p, with p(1 − p) as climatology, per AOI; mean and median
   with equal weight per AOI.

## Addendum 4 (record of the implementation of B6, B8 to B12, B14 and B15; 2026-10-29T12:57Z, commit c394b67, before these quantities were computed)

Script `backend/verification/r_okauchi/v70b_b13_part2.py` (output `docs/verification/v70_b13_local2.json`).

1. **B6**: influence is the jackknife over location groups (models fixed), ψ = R − R(−g). Effective number of clusters =
   (Σψ²)² / Σψ⁴. R without each of the 10 most influential groups, for the confirmatory and the corrected run. For B1(b),
   the 10 most influential groups of the corrected run are removed one at a time with the second stage refitted (design
   of B-8 C2).
2. **B8**: logistic regression over the held-out cells of all folds (all cells) on the standardised [base logit,
   composition, name features] (one intercept, v39.fit), cluster-robust by location group, constant columns dropped;
   joint Wald test that all name coefficients are zero.
3. **B9**: the uncapped EMS estimate (table, corrected readings and 2,000 pairs replicates as in B-11 B2) and the planned
   GFM estimate (table of v41, cluster-robust). (a) the rounding rule of v37 applied to the GeoNames coordinates of the
   nearest settlement. (b) points (frequency-weighted IRLS) without each activation that has water-related cells within
   500 m, and the pairs interval without EMSR052 and EMSR130. (c) AOIs whose record has schema_generation 0 or 1 removed.
   (d) within each AOI, the water-related flags of the (AOI, nearest settlement) pairs permuted among the pairs (200
   times); points of the uncapped EMS estimate (IRLS) and of the planned GFM estimate; p two-sided on |log OR| as
   (1 + count) / (1 + 200). (e) cells of the planned GFM table, by distance band and by tercile of the 10 m built-up share
   (WorldCover lc_built). (f) the capped table (v41, EMS labels), AOIs with a known sensor (Sentinel-1 and other SAR =
   radar, optical), term water-related × within 500 m × optical. (g) cells whose nearest settlement has the feature code
   PPLX, PPLH, PPLQ or PPLW removed (both estimates).
4. **B10**: analysed settlements = settlements nearest to at least one cell of the uncapped table. Frozen lexicon with the
   corrected loader. Positions counted per matched entry. Unicode names read with ä, ö, ü, å, ø folded to one letter,
   æ → ae, ß → ss and combining marks removed.
5. **B11**: shares of the N-read classes for the settlements of B10, by the blocks assigned to the country (joined with +)
   and by country.
6. **B12**: a public file `docs/verification/v70_two_coders.csv` with both coders' codes beside the GeoNames ID, name,
   country and stratum. Names in any audit sample (v51, v68, v69) are left out (so that no audit answer is revealed).
   Agreement from v57b.
7. **B14**: flooded shares of EMS and GFM by AOI on the planned GFM table; the terrain gain (Brier of M0 minus Brier of M1)
   by AOI in the run trained and scored on GFM labels (world_v2_gfm).
8. **B15**: as in item 8 of addendum 1. Settlement values come from the capped table (the table of v37, the same cells as
   v41), using the cells whose nearest settlement is that settlement.

## Addendum 5 (record of the aggregation of C and of the near-settlement odds ratios with the corrected lexicon; 2026-10-29T13:14Z, commit 44fbbfc, before these quantities were computed)

Script `backend/verification/r_okauchi/v71_b13_cluster.py` (output `docs/verification/v71_b13_cluster.json`).

1. **Aggregation of C1 to C3**: from the held-out predictions of each run, R, intervals and the verdict of the decision
   table with the location-group bootstrap of the confirmatory analysis (4,000 replicates, `_verdicts` of v56). C1 gives
   the main comparison of names (M1c → M3n, placebo M3np) and each control (per cell and through settlements). C3 adds the
   Brier score of the constant prediction (the training-side flooded share of each fold) and the distribution of the
   chosen numbers of iterations.
2. **Near-settlement odds ratios of C2**: with the water-related flag of the nearest settlement read again with the
   corrected lexicon, the capped EMS estimate (table of v37, estimator of Sect. 5, cluster-robust), the planned GFM
   estimate (table of v41) and the uncapped EMS estimate (table and 2,000 pairs replicates of B-11 B2).
3. **C4 (B-12 C1)**: from the gain curve of shape BK, the injected odds ratio that reaches the margin (the "reach" reading
   of v59) and the odds ratios reached by the upper limits of the fixed and refitted intervals (for each of the four
   no-effect references), read as in C3 of v64.

## Addendum 6 (correction of the implementation of B1(c) and record of the aggregation of C; 2026-10-29T14:20Z, commit 07c9788, before B1(b)(c)(d) and B2 were computed)

1. **How countries are drawn in B1(c)**: the table of country codes (the codes of v68 B7) has 70 countries, but once each
   activation is assigned to the country with most of its cells (addendum 3, item 1), 55 countries hold an activation;
   the other 15 only have cells in AOIs of activations assigned to other countries. The first implementation drew from
   the 70 countries and stopped with a division by zero in a replicate that drew a country without activations (after
   about 95 minutes of computing; no results were saved). "Countries drawn with equal probability, as many as there are
   countries" (addendum 1, item 7) is corrected to drawing 55 times with equal probability from the 55 countries that hold
   an activation (`topsem.world.b13.nested_weights`). Activations are redrawn within each drawn country as in addendum 3.
2. **Saving as the computation proceeds**: the replicates of (b) and (c) are saved by scheme
   (`data/verification/v70_b1_refit_reps_{groups,nested}.npy`) and read back if present. The numbers of replicates, the
   seeds and the refitting rule are unchanged (seed SEED + b for the 500 replicates of (b), SEED + 100,000 + b for (c)).
3. **Record of the aggregation of C1 to C4**: the aggregation of C3 in v71 stopped because the prediction files carry no
   fold number; the fold number is built from the lengths of the files in the order of load_split (file names), and the
   constant prediction is as in addendum 5. In C1(b), fold 8 was truncated when the predictions were fetched from the NAS;
   its two files were fetched again and every fold was checked to load before the aggregation.

## Addendum 7 (part of the B9 pairs resampling on the cluster; 2026-10-29T16:40Z, commit f25fffe, before B9(b)(c)(g) were computed)

1. **Where**: the pairs resampling (v68.pairs_or, 2,000 replicates, the seed of v68) of B9 (b) the table without EMSR052 and
   EMSR130, (c) the table without schema generations 0 and 1, and (g) the table without PPLX, PPLH, PPLQ and PPLW takes
   more than two hours per resampling on the local machine with two workers, so it is computed on the cluster. Tables,
   subsets, seed and number of replicates are those of v70b (`backend/verification/r_okauchi/v72_b9_pairs_cluster.py`
   builds the same subsets with the functions of v70b and calls v68.pairs_or with the same keys). The uncapped table is
   `data/verification/v68_uncapped.parquet` on the NAS (the same file as on the local machine; its sha256 is recorded).
2. **Aggregation**: the replicates of each resampling (`data/verification/v68_pairs_v70_b9{b_no052_130,c_gen2plus,g_fcodes}.npy`)
   are fetched from the NAS, and B9 of v70b on the local machine reads them (it does not recompute a key that is saved).
   (a), the leave-one-activation-out points of (b) and (d) are computed on the local machine.
3. **Machines**: the cluster is amd64 and the local machine arm64. IRLS results can differ slightly between machines; the
   difference from (a), computed locally, is expected not to show in the digits of the interval limits but is not
   checked (in the primary comparison the two agreed to the fifth decimal of R).
4. **Implementation record**: the first submission (image 1793292010673138000) stopped while loading, before computing
   anything, because v36 finds the repository root by CLAUDE.md and the image had no CLAUDE.md. The Dockerfile now copies
   it, and the jobs were resubmitted with the same arguments.
