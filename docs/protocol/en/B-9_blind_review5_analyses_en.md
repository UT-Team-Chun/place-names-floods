# B-9 Corrections to the reporting and exploratory analyses in response to the fifth round of first-time critiques (NHESS version, 2026-10-28)

English translation, drafted with a language model; the Japanese file `docs/protocol/amendments/B-9_blind_review5_analyses.md`
is the record and governs where the two differ.

Status: **plan (about 2026-10-28T19:00Z, before any quantity in B and C below was computed)**. The confirmatory verdict
is not changed. Everything below is **post hoc** and responds to an editor and three referees who read only the PDFs
of the manuscript and its supplement for the first time (language-model agents; reports in
`docs/research/reviews/2026-10-29_paper2_blind5_*.md`, commit cb1995d). The paper describes these analyses as "fixed in
this document after the confirmatory result had been seen and before they were computed". Protocol:
`preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924). Code: with the correction of amendment B-4
(corrected code).

Values seen when this document was written: all results of the confirmatory analysis and of B-2 to B-8 (as printed in
the version of the paper at commit e969593), and values that the reports point out or estimate. In particular: re-reading
the bound with the refitted interval gives an odds ratio of about 1.57 and an odds ratio detectable with 80 % power of
about 2.0; combining the 13 fold assignments gives a bound of about 1.3; the bound lies below the smallest simulated
effect (δ = 0.25), and with a locally quadratic curve 1.23 would become about 1.25 and the zero-reference bound 1.12
about 1.18 (referee 2's estimates); across the 95 % interval of the yardstick the margin corresponds to odds ratios of
about 1.8 to 2.9 (referee 2); on cells that GFM can assess in full, both label sets give 0.67 (Table S39, computed under
amendment B-3, addendum 2). None of the quantities in B and C below has been computed.

## Main points of the critiques (editor: return before review; three referees: major revision)

1. The title, abstract and conclusions present the null of the pooled test as an answer to the whole question. State in
   Sects. 1–2 (not first in Sect. 4) that the pre-specified test cannot register an effect confined to the ground near
   settlements, and lead with the bounds (as ranges) rather than the verdict label
2. The bound and the detectable odds ratio understate uncertainty (fixed-model interval, one fold assignment, a straight
   line below the grid, injections that only add flooding, the same effect in every group). The margin is conventional
   and depends on the uncertain yardstick
3. No person has validated the readings of the names (toponymist coding, rates by language). The water-side class mixes
   elements whose expected direction differs
4. The labels and the terrain model are weakest near settlements (permanent water treated differently between AOIs, a
   surface model, HAND from a dense drainage network). Report the terrain model's skill in flood-mapping terms (CSI etc.)
   in the main text
5. The provenance cannot be checked by third parties, the post hoc analysis is too large, and the manuscript is too long
   and dense (halve it). Deposit the derived data now
6. Mechanical errors: every supplementary table number is off by one (the contents longtable takes Table S1); the count
   of Fig. 4b (the GFM decision is a dictionary); the extent of the Norwegian defect (2.2 % of the cells, not 0.34 %);
   mixed interval types for GFM; the expected directions of the element classes were stated after the estimates
   (different from those registered in amendment B-5, item 8); the result on cells that GFM can assess in full (0.67 on
   both label sets) is not in the text; the flexible refit on the language-model split (not-meaning) is not in the text

## A. Corrections to the reporting (no computation; in the paper and the supplement)

1. **Mechanical errors**: supplementary table numbers (reset the table counter after the contents); the count of Fig. 4b
   (read the label of a decision stored as a dictionary: 17 rows, 14 with a verdict, 2 other than equivalent-null); the
   extent of the Norwegian defect (cells whose nearest name the lexicon did not read: 2.2 %, in 11 AOIs; the indicator that
   the lexicon covers the language also acted as an indicator of Norway); GFM intervals within 500 m from the pairs
   bootstrap throughout; notes in the supplementary tables on the three EMS estimates on the GFM sample (0.78, 0.79, 0.80);
   the code version in the supplement's preface (Tables S7 and S8 come from the frozen code); the grid in the note
   of the supplementary table (> 20)
2. **Framing**: the title, abstract, introduction and conclusions state that (i) the test bounds an area-wide effect,
   (ii) of these readings of the names, (iii) on flood extent as EMS mapped it. Sects. 1 and 2 state that the
   pre-specified test cannot register an effect confined to the ground near settlements. The abstract uses no verdict
   label ("within the equivalence margin" instead), states that the protocol was not registered with a third party and
   that the estimand was chosen after a relation between the composition of names and flooded fractions per area had
   been seen, and gives few numbers (the main result, the range of bounds, the design lesson) with at most one or two
   exploratory quantities
3. **Expected directions by class**: the estimates by element class are reported against the expectations registered in
   amendment B-5, item 8 and its addendum (wetland: more flooding and lower; islands: higher; springs: unrelated;
   crossings, dikes and lowland: no direction). Directions stated after the estimates (the third paragraph of the
   introduction) are not used. The evidence is summarised by the joint tests (classes × distance bands, p = 0.13;
   water-side × distance bands, p = 0.89)
4. **GFM**: the result on cells that GFM can assess in full (amendment B-3, addendum 2, planned after the GFM result had
   been seen; 0.67 on both label sets) is reported in the text, and "on the ground that both label sets can assess, the
   automated labels record relatively less flooding near water-side names" is restricted to cells that GFM sees in part
5. **Every pre-specified comparison in the text**: the flexible refit on the language-model split (lexicon and character
   model, not-meaning), and the share of the terrain information recovered at T0 in Table 2. H4 is described as carrying
   little information (the verdict stays). Table C2 also gives 95 % intervals
6. **Wording**: "equivalent to zero" → "within the equivalence margin"; "changed the predictions" → "changed the held-out
   Brier score"; the refitted bootstrap refits only the second stage (500 replicates); the Spearman trend test over five
   points has almost no power; the coverage check does not include the variability from which groups the models are
   trained on; inference concerns AOIs of this kind (EMS products chosen by analysts); the flood-information controls
   change both the information and its assignment; "the information that names carry once position is known is
   regional" is tempered; the use of a list of names as a flag concerns the water-side class; the extension to oral
   histories and community maps is tempered; the heading of Sect. 7.3 is narrowed; the abstract and Sect. 7.4 state that
   gazetteer names are not local knowledge in the sense of the literature (practices, oral memory, minor and field names)
7. **Methods details**: the default early stopping of gradient boosting (M0 200 iterations, M1 78–99; v40); the 50 m
   sub-grid of large AOIs (number of AOIs); cells at least half covered by non-riverine flooding are excluded (number);
   where a product has no hydrography layer, channel cells enter as dry cells; the definition of the settlement-layout
   features; the penalty of the second stage (L2, C = 1, standardised features); the target of the character model; the
   identifier and date of the second coder's model (from the records); where the first coder's rules are (public
   repository)
8. **Presentation**: shorten the main text by about half (exploratory details to the supplement; the sentences on the
   multiple of location groups needed leave the main text); figures and tables near their first citation; model codes
   in plain words; the rules choosing the AOIs of Fig. 3 in the supplement; the different colour scales of Fig. 3; grid
   points in Fig. 6; the colour of the marker in Fig. 8a; citations for the methods (Wilson, Kish, Holm, pairs cluster
   bootstrap, jackknife with t quantiles, probabilistic bias analysis, Murphy decomposition); mark the preprint and the
   conference abstract; check and cite Fanta et al. (2019)

## B. Computations from existing artefacts (local; script v63)

1. **Uncertainty of the bound and the detectable odds ratio** (referee 2, 2(a)(b); referee 3, 5): on the 50-replicate
   curves of amendment B-8, item B1 (on the combined grid once the low-δ grid of C1 is available), read the bound again
   with the upper 90 % limit of the observed increment replaced by (i) the upper 90 % limit of the interval with the second
   stage refitted inside the bootstrap (B-8, C2, corrected code, 500 replicates), and (ii) the upper limit of the 90 %
   interval that combines the spread of R over the 13 fold assignments (corrected code, v47) with the bootstrap error (the
   value of Sect. 3.4). The odds ratio detectable with 80 % power uses the standard deviation of (i). The no-effect
   references are the four of B-8 (primary: the mean of the reassignments with rebuilt composition of B-7). Monte Carlo
   intervals as in B-8, B1. **Rule**: the bound in the text and in Table 3 is given as a range from the fixed-model value
   to the refitted value, not as a single value
2. **Dependence of the margin's translation on the yardstick** (referee 2, 3): with 5 % of the 2.5th, 50th and 97.5th
   percentiles of the bootstrap distribution of the yardstick (ΔBS from M0 to M1) as the margin, read the odds ratio on
   the curve of form A
3. **Calibration and information** (referee 2, minor point 8): the Murphy decomposition (reliability, resolution,
   uncertainty; the same binning as B-8, B8) of the Brier scores of M1c and M3n in the primary comparison (corrected code),
   and 95 % intervals of their differences from a cluster bootstrap over location groups (4,000 replicates, seed 20261029)
4. **Within-AOI evaluation** (referee 1, 4(d)): shift the logit of every prediction within each AOI by one constant so that
   the mean prediction equals the observed flooded fraction (the AOI intercept set from the answer), for M0, M1, M1c and
   M3n; R_within = ΔBS(M1c → M3n) / ΔBS(M0 → M1) with cluster-bootstrap 90 % and 95 % intervals. The increment about where
   flooding occurs, with the area-level amount removed
5. **Flood-mapping skill** (referee 1, 4(c); referee 2, 3): CSI, hit rate and false-alarm ratio when each AOI is mapped on
   as many cells as flooded, in order of the prediction (as in B-8, B9). For HAND alone, cells are ordered by minimum HAND,
   ties by mean HAND, remaining ties at random (seed 20261029), replacing the stable sort of B9. Add the AOI median of the
   difference in CSI (M3n − M1c) with a cluster bootstrap over location groups (4,000 replicates, 95 % interval). **Rule**:
   CSI, hit rate and false-alarm ratio go into the main text, and the difference in CSI is reported beside R as the effect
   of names in absolute units
6. **Leaving out one activation at a time** (referee 2, minor point 9): refit the main model of B-8, B6 (water-side ×
   500 m, EMS) leaving out each activation that has flooded cells in turn; report the minimum and maximum and the
   activations that give them
7. **Multiplicity of classes**: Benjamini–Hochberg q values of the six class coefficients (within 500 m, the class model of
   B-8, B6) beside the Holm p values
8. **Detail products** (referee 1, 2(d); referee 3, 3(a)): from the frozen product inventory, the number of analysed AOIs
   with a detail-scale product and with both detail and overview products (descriptive)

## C. Computations on the cluster

1. **Low-δ grid** (referee 2, 2(c)): the plasmode with base models fixed (the implementation of B-6, item 2, and its rule
   for the replicate salts) at δ = 0.05, 0.10, 0.15 and 0.20, forms A and CA, replicates 0–49. The bounds of B1 above and of
   B-8, B1 are read again on the grid combined with that of B-6
2. **Injections that remove flooding** (referee 2, 2(d)): set target cells with label 1 to 0 with probability 1 − p′/p
   (p′ = p e^δ / (1 − p + p e^δ), δ < 0; the odds of flooding given terrain are multiplied by e^δ). δ = −0.25, −0.5,
   −0.75, −1.0, −1.5, −2.0, −3.0, forms A and B, replicates 0–49. The δ at which the gain reaches the upper 90 % limit of the
   observed increment minus the reference is read as a lower bound (the lower side of the odds ratios that the data
   exclude). **Rule**: the bounds in the text are two-sided (for example, "excludes effects with odds ratios below x or
   above y")
3. **Effects that differ between settlements** (referee 2, 2(d); referee 3, 5(b)): form AH: inject on all ground nearest a
   water-side name with δ_s = δ w_s, where w_s is an exponential weight per settlement (mean 1, set by the settlement's gid
   and the replicate's salt). δ at the seven values of B-6, replicates 0–49. The bound is read on the mean δ: how much the
   gain changes for the same mean effect when the effect is concentrated in some settlements
4. **Simulated power with every model refitted** (referee 2, 2(d)(e)): inject odds ratios 1.25, 1.5 and 2.0 (δ = ln) on
   form A and refit every model (the implementation of B-8, C3), 40 replicates each (replicates 200–239). Apply the
   decision table without the placebo condition. Report the share of positive verdicts (lower 95 % limit above zero) with
   Wilson 95 % intervals, the share of equivalent verdicts and the mean R. **Rule**: the smallest odds ratio whose share of
   positive verdicts is at least 0.8 is reported as the simulated detectable odds ratio, beside the normal approximation
   (1.56); if none reaches 0.8, write "above 2.0"
5. **One definition of permanent water for all AOIs** (referee 1, 2(a)): exclude, in every AOI, cells whose mean JRC GSW
   occurrence over the 25 points of v58 (cell_gsw.parquet) is at least 50 % (the exclusion by the products' hydrography
   layers stays), and recompute the primary comparison (corrected code) and the near-settlement odds ratio (the main model
   of B-8, B6, pairs bootstrap with 2,000 replicates). Report the number of cells excluded, separately for AOIs with and
   without a hydrography layer. **Rule**: a verdict other than equivalent-null does not replace the confirmatory verdict,
   and is reported as such
6. **Seed of the subsampling under the cap** (referee 2, minor point 13): the primary comparison (corrected code) with the
   seed of the within-group subsampling (frozen: 20261025) set to 1, 2 and 3. Report the range of R
7. **Bare-earth elevation (FABDEM v1-2)** (referee 1, 4(a)(b); editor, 3.3): **condition for running**: FABDEM elevations
   are available for at least 95 % of the analysed cells (tiles fetched on the cluster and kept on the NAS). If met, the
   elevation features (1 km mean, mean, depression, minimum, standard deviation, slope) are rebuilt from FABDEM, and in
   place of ASF HAND the model receives the height above the nearest vertex of a HydroRIVERS reach with an upstream area of
   at least 100 km² and of at least 1000 km² (the cell's mean elevation minus the minimum elevation of the cell containing
   that vertex, as in B-8, B7). With this terrain model, recompute the primary comparison (corrected code; the position
   model M0 unchanged), the near-settlement odds ratio (pairs bootstrap with 2,000 replicates) and the flood-mapping skill
   (B5 above). **Rule**: report R and the verdict, the change in the yardstick and the change in CSI. The near-settlement
   odds ratio is an exploratory estimate; an interval excluding one is reported as uncorroborated and does not change the
   rule of B-8, B6 (the ground near settlements remains unresolved). If the condition is not met, report that, with the
   share of cells covered

## D. Not done in this amendment (authors' decisions, or a separate study)

1. Coding of the 100 names by toponymists, and expert review of the main lexicon blocks (to be commissioned by the
   authors)
2. Third-party archiving of the development history including the freeze; release of the pilot and design documents and
   of the protocol and results of the Japanese test (authors' decision)
3. Deposit, with DOIs, of the derived data (cell tables, readings of the names, held-out predictions, run summaries, the
   GeoNames snapshot) (by the authors)
4. Checking the English translation of the protocol (by the authors)
5. Defining settlements by built-up extent (GHSL, WSF), building labels from detail products, national multi-event records,
   Holocene alluvium, stratifying by protection standard (FLOPROS), and comparing pairs of names that differ only by an
   upper/lower qualifier: for the settlement-level study (a separate study registered with a third party). Adding them
   here would enlarge the families of post hoc estimates further

## Implementation and resources

- `world_eval.py`: `--plasmode-deltas` (list of δ), `inject_labels` with δ < 0 (labels 1 set to 0) and with an array of
  δ, form AH (weights per settlement), `--exclude-gsw-perm` (threshold of the mean GSW occurrence), `--table-seed` (seed
  of the subsampling under the cap), `--terrain fabdem` (FABDEM terrain features and heights above rivers). Tests confirm
  that injections with δ > 0 give the same values as the present code
- B is computed by `v63_b9_local.py`, C is aggregated by `v64_b9_aggregate.py`. Random seed 20261029
- C1–C3 are split into jobs by fold and range of replicates as in B-8, C1; C4 as in B-8, C3

## Addendum 1 (about 2026-10-28T22:03Z, before any result of C1–C4 or C7 was seen): implementation record

The results of C5 and C6 were seen on 2026-10-28. This addendum records implementation details; it changes neither the
computation nor the reporting rule of any quantity.

1. **Stopped runs of C5**: world_v2_b9gswperm stopped twice before computing anything on its first submission: (i) the
   GSW cell table (`gsw/cell_gsw.parquet`, the 25 points of v58) was missing on the cluster's storage (copied there from
   the local machine), and (ii) with `--exclude-gsw-perm` the code did not read the grid position of the cells (ix, iy;
   KeyError, fixed in ccb0a56). It was resubmitted with the same arguments.
2. **Retrieval of FABDEM for C7**: the provider throttles each connection to about 0.5 MB/s, so instead of whole 10° zips
   only the needed 1° tiles (821 tiles intersecting each AOI ± 0.2°) were extracted from the zips by HTTP range requests
   and stored in `fabdem/tiles/` on the cluster's storage (805 tiles; 16 sea tiles absent from the zips; CRC and size
   checked against the zips). The tiles are byte-identical to those in the zips, so the computation is as planned
   (21ac5b4). The run condition (FABDEM elevation available for at least 95 % of the analysed cells) was met at 99.96 %
   (1,506,916 of 1,507,486; `fabdem/coverage.json`).
3. **Cluster jobs recreated**: the 159 jobs of C1–C4 that had never started were recreated with smaller memory requests
   (6 GiB for injection runs, 8 GiB for power runs; limits unchanged) and the same commands, arguments and image. The old
   jobs stay suspended (mapping in `docs/verification/v64_b9_jobs_recreated.json`).
4. **Additions to v64**: the flood-mapping skill of C7 (the same function as B5, applied to the held-out predictions of
   the FABDEM run; the HAND row uses the same GLO-30 HAND as v63 as a reference) and the curves of Fig. 6 (mean gains and
   Monte Carlo standard errors on each form's grid for B1, and standard errors for C2 and C3). Run on the corrected run,
   the B5 function reproduced the values of v63 (critical success index of the terrain model 0.129).
5. **Reporting** (no computation changes): in the main text, amendment numbers are dropped for exploratory analyses that
   were planned before they were computed, and each section cites the table of planning status (Supplement Table S4)
   once; marks for analyses planned after results were seen, not planned, rules fixed in advance and recorded
   expectations are kept.
