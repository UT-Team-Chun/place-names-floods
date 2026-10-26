# B-5 Post hoc analyses in response to the first blind review (2026-10-26 evening) — English translation

> English translation (2026-10-27) of `docs/protocol/amendments/B-5_blind_review_analyses.md`. The Japanese file is
> the record; where the two differ, the Japanese file governs. Translator's notes are in [square brackets].

Status: **planned (2026-10-26 evening, before the runs and computations below)**. The confirmatory verdicts are not
changed. Everything below is **post hoc**, added in response to three referees and an editor who read only the PDF
for the first time (reports in `docs/research/reviews/`, `2026-10-26_paper2_blind1_*.md`). The referees only read
the manuscript and the public repository; none of the quantities below has been computed by anyone. The paper
describes them as "post hoc analyses fixed in an internal amendment before computation", distinguished from
registration with a third party (as the referees asked).
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924).
The analyses use the code with the B-4 correction (the Norwegian lexicon), and comparisons with the confirmatory
analysis are set beside the B-4 rerun (world_v2_lexfix).

The author's decision (2026-10-26): paper_2 narrows its claims to what the test answered and goes to a specialist
journal. The settlement-scale question is addressed in a separate study registered in advance with a third party.
The analyses below are for paper_2.

## Analyses (the numbers are used in reporting)

1. **Permutation distribution** (referees 2 and 3): on the primary split, a placebo that permutes names within AOIs,
   200 times (seeds 1–200, distinct from the seed of the confirmatory analysis). The base models are fitted once and
   the second stage (the M1c → M3n form) is refitted for each permutation. Report where the real increment lies in
   the permutation distribution (one-sided rank) and the 2.5–97.5% range of the permutation distribution
2. **Second-stage coefficients** (referees 1, 2 and 3): save the coefficients of the second-stage regressions of each
   fold of the primary split (M1c, M3n, M3np, M3d) and report, for each name feature, the mean and range over folds
   (standardised scale)
3. **Repeated fold assignments** (referee 2): on the primary split, fold-assignment seeds 3–12 (10 assignments; 13
   together with the existing 0, 1 and 2). Arms M0, M1, M1c, M3n, M3np. The distribution of R and an interval
   combining the bootstrap SE with the SD over assignments
4. **Upper bound on the odds ratio by test inversion** (referees 1 and 2): with the existing synthetic labels (the
   grid of δ on the primary split), the upper bound is the largest δ for which the observed ΔBS does not fall below
   the 5th percentile of the distribution of ΔBS over the replicates at that δ (linear interpolation between δ).
   Both for the effect over the whole area and for the effect within 500 m of the settlement only. Stated to be
   conditional on the synthetic model
5. **Oracle control** (referee 1): instead of the name features, enter the terrain at the site of the nearest
   settlement (the mean, over the settlement's points, of the minimum HAND of the cells within 177 m of the
   settlement point, assigned to cells from the nearest settlement as for names) and compute an increment of the same
   form as M1c → M3n on the primary split. A guide to how much a "name" that is perfectly readable and correctly
   describes its place could add in this design. Reporting rule: if this increment does not reach the margin
   (R = 0.05), write that the area-wide design can hardly pick up information about the settlement site
6. **Within-area metrics and influence of places** (referees 1 and 2): from the existing held-out predictions, (a)
   the distribution of the per-AOI differences in AUC and PR-AUC between M1c and M3n, (b) the difference in
   conditional log-likelihood with the intercept refitted per AOI (bootstrap over location groups), (c) the range of
   R leaving out one location group at a time, and the places that contribute most to R
7. **Label variants** (referee 1): on the primary split, labels rebuilt in two ways — (a) flood traces are not
   counted as flooded, (b) the flooded share is measured over the area excluding permanent water — and run with the
   same arms as the confirmatory analysis
8. **Directional analysis by element** (referee 3): split the entries of the frozen lexicon into groups whose
   expected direction differs in toponymy — wetland (the wetland class), island and river island (island, holm,
   werder, au, waard, etc.; dry rises), spring and well (well, born, brunn, etc.), crossing (crossing), dike and
   reclamation (dijk, dam, polder, etc.). For each group, report the within-AOI odds ratios within 500 m of the
   settlement and over all cells (the same regression as v39, adjusted for terrain) and the positive control (mean
   and percentile of HAND). Expectations: wetlands flood more and lie low, islands sit high, springs are unrelated to
   flooding. Treated as description; conclusions are limited to groups whose Holm-corrected intervals exclude 1
9. **Sensitivity excluding settler countries** (referee 3): the main comparison without the AOIs of the countries of
   the Americas and Oceania (from the held-out predictions, bootstrap over location groups). A guide to the influence
   of transferred and commemorative names
10. **Comparison with etymologies recorded by people** (referees 1, 2 and 3, in place of an expert check): collect the
    etymologies of the names of the settlements of the analysed AOIs from settlements whose Wikidata "named after"
    (P138) points to a water body (river, lake, marsh, spring, etc.) and from the etymology sections of Wikipedia
    (including language editions). An LLM is used only as a tool to extract the etymology from the text, keeping the
    source URL and a quotation. Per language, the precision (share of readings that agree with the record) and recall
    (share of names whose recorded etymology is water-related that were read as such) of the lexicon and of the
    language model. State the bias towards names with a human record, and that this is not an expert check

## Reporting rules

- All the quantities above are reported (SI). The text gives a summary of 1, 3, 4, 5 and 10 and the result of 8
- The verdicts of the confirmatory analysis are not changed. If the results above disagree with the confirmatory
  analysis, this is written as it is
- The result of 5 is the basis of the discussion of "the scope of the test" in the text (in either direction)
- The order of computation, thresholds and group definitions are not changed after seeing results. If they are
  changed, an addendum is added to this document and the results before the change are also reported

## Addendum 1 (2026-10-26 evening, fixing the forms of the groups of item 8 before computation)

The groups of item 8 take only the following (class, form) pairs from the entries of the frozen lexicon. Wetland,
crossing and lowland are the lexicon classes as they are. Island, spring and dike are forms selected by meaning
within the lexicon classes (the same spelling across languages is treated the same). Forms with split meanings
(English bank, German ried, etc.) are not included. German au, aue and auen are included in island, matching the
registered text above and the referee's list (they also mean riverside meadow, which is noted when reading the
result of this group).

- **Island and river island** (from riverine): eiland, djezira, gezira, jazira, jazirat, jezira, chilumba, chisumbu,
  ostrov, holm, holmen, insel, werder, wert, werth, woehrd, woerth, wohrd, worth, au, aue, auen, nisaki, nisi, nisos,
  ait, eyot, holme, island, isle, isla, islas, saar, saare, luoto, saari, ile, iles, ilot, sziget, pulau, isola,
  isole, isolotto, kaoh, koh, aral, araly, sala, salos, salas, nosy, kyun, waard, waarden, weerd, oy, oya, ostrow,
  ostrowek, ostrowo, kepa, kepy, ilha, ilhas, ilheu, insua, insula, ostrovu, ada, ostrvo, otok, jasiirad, ishull,
  ishulli, holmar, kisiwa, ko, pulo, adasi, ostriv, ostrovok, con
- **Spring and well** (from water): fontein, put, ain, aioun, ayn, ayoun, abar, bir, biar, hassi, cheshma, izvor,
  izvori, kladenci, kladenets, bunar, pramen, prameny, studna, studnice, studnicka, kilde, born, brunn, brunnen,
  quell, quelle, kefalovryso, krini, pigadi, pigadia, pigi, vris, vrises, vrisi, vrys, vryses, vrysi, borne, bourne,
  spring, springs, well, wells, waterhole, soak, fuente, fuentes, manantial, ojo, pozo, pozos, allikas, chah, chahi,
  chashma, chashmeh, lahde, fontaine, fontaines, puits, source, sources, forras, kut, sumur, fontana, fontane, fonte,
  pozzi, pozzo, sorgente, bulak, bulaq, saltiniai, saltinis, avoti, avots, loharano, ameyal, studnia, studzianka,
  zrodla, zrodlo, cacimba, fontes, nascente, poco, pukyu, puquio, fantana, fantanele, izvoare, izvorul, studenac,
  vrela, vrelo, vrutak, izvir, studenci, studenec, ceel, burim, burimi, pus, puse, kalla, kallan, chemchemi, kisima,
  bukal, kaynak, kuyu, pinar, dzherelo, kolodez, kolodiazi, krinitsa, krynytsia, gieng
- **Dike and reclamation** (from riverine): wal, hraz, dige, levee, penger, digue, tanggul, argine, argini, bonifica,
  pylimas, dambis, dam, dijk, polder, sluis, dig, nasip, set
- **Wetland**: all forms of the class wet. **Crossing**: all forms of the class crossing. **Lowland**: all forms of
  the class low (reference)

Matching uses the same rules as the confirmatory analysis (entries of the country's languages, position rules) with
the reading corrected by B-4. Expectations by group (as registered): wetlands flood more and lie low, islands sit
high, springs are unrelated to flooding. No direction is expected for crossing, dike and lowland.

## Addendum 2 (2026-10-26 evening, fixing the procedure of the Wikipedia part of item 10 before computation)

The first half of item 10 (Wikidata's named-after P138, v43) is running, and its interim progress shows a record of
the namesake for only about 0.8% of the Wikidata items (only counts were seen; the comparison with the readings has
not been computed). The second half (the etymology sections of Wikipedia) is therefore done by the following
procedure (v44). The comparison between the readings (lexicon, language model) and the records is not computed until
all records have been classified.

1. **Population**: of the nearest settlements (gid) of the analysed table (the v37 table), those with a Wikidata item
   (found by v43 via P1566) and an article in the English Wikipedia or in the edition of the first language of the
   country in the lexicon
2. **Sampling**: 300 at random each (seed 20261026) from the settlements that the lexicon read as water-side
   (lex_wetside, with the B-4 reading) and from those it did not. From the articles of both editions, take the text of
   the sections whose heading denotes the etymology or origin of the name (a list of headings per language is fixed in
   the code) (at most 2,000 characters per article). Settlements without such a section are "no record" and are not
   classified
3. **Classification** (extraction by an LLM): the excerpts are placed in a file in random order with all information
   other than a serial number hidden (the readings of the lexicon and language model, the stratum, the country), and
   Claude Opus 5.5 (the model of this session) classifies them by the following rules. The meaning described as the
   origin of the name is recorded as a set of the lexicon's classes: water (river, stream, lake, pond, spring, well,
   waterfall, the name of a river itself), wet (wetland, marsh, peat bog, mud, reed bed, wet meadow), riverine (river
   island, riverside meadow, floodplain, river bank, sandbank, reclaimed land, dike), crossing (ford, bridge, ferry,
   landing place), low (valley, hollow, lowland, "lower"), high (hill, mountain, height, rock, "upper"), other (a
   person, a saint, a tribe, transfer from another place name, a plant, an animal, a building, an occupation, etc.).
   A record that does not state the origin of the name is "no derivation", and one whose theories are split across
   the water side (water, wet, riverine, crossing, low) and the other side is "unclear"; both are excluded from the
   metrics and their numbers reported. For each classification, a quotation of the original text (at most 300
   characters) and the URL are kept
4. **Metrics**: with "the recorded origin is water-side (any of the 5 classes above)" as the truth, the precision of
   the lexicon (share water-side in the stratum it read as water-side, Wilson interval), the share water-side in the
   stratum it did not read, and the recall weighted by the size of the strata in the population (interval from 2,000
   bootstrap resamples within strata). For the language-model reading (nread_wetside), the stratum-weighted precision
   and recall in the same sample. By language group (only groups with at least 30 classified records). For each
   lexicon class, the share whose recorded class agrees (description)
5. **Bias**: state that the sample is limited to settlements with an article and an etymology section (large, old
   towns), and that the classifying LLM is the same model that assisted the analysis and knows the content of the
   lexicon (the individual readings were hidden). This is not an expert check

## Addendum 3 (2026-10-26 evening, before submitting the repeated fold assignments of item 3)

Because B-5 is done with the code including the B-4 correction (the Norwegian lexicon), the 13 assignments of item 3
are also aligned to the corrected code. Seed 0 is the B-4 rerun (world_v2_lexfix), seeds 1 and 2 are rerun with the
corrected code (world_v2_fs1b, world_v2_fs2b), and seeds 3–12 are new (world_v2_fs3 … world_v2_fs12). The arms of the
second stage are the same set as the existing runs of seeds 1 and 2 (M0s, M0c, M1s, M1c, M3n, M3np, M3s, M3ps)
(including the registered M1c, M3n and M3np; the others only add the same regression and do not change the registered
quantities). The values of seeds 1 and 2 with the code before the correction (shown in the earlier version of the
paper) are also reported alongside. None of the results has been computed yet.

## Addendum 4 (early 2026-10-27, enlarging the sample of the Wikipedia part of item 10 before any classification)

With the procedure of Addendum 2, the population (settlements with a Wikidata item and an article) was 16,150 (1,021
in the stratum read as water-side, 15,129 in the other stratum), and only 102 (17%) of the 600 drawn had an etymology
section. This is too few for the precision by stratum and the values by language group, so **at this point, before
any excerpt was classified (the content of the file to be classified had not been opened)**, the sample is enlarged as
follows. Only the counts above were seen.

- The stratum read as water-side is taken in full (1,021); 1,500 are drawn at random from the other stratum (seed
  20261027; this does not necessarily include the first 300). Settlements whose excerpts were obtained in the first
  retrieval use the same excerpts (they are not retrieved again)
- Classification, measures and reporting follow Addendum 2 (the stratum weights are computed from the sample sizes and
  the sizes of the strata in the population)
- No values are computed on the first sample of 600 (only the enlarged sample is reported). The reason for the
  enlargement and the first counts are reported
