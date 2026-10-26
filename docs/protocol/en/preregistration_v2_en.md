# Protocol v2 — Do place names remember floods? A global test against satellite-observed inundation

> **English translation of the frozen protocol.** The frozen record is the Japanese file
> `docs/protocol/preregistration_v2.md` (sha256 `4d9870ad888ebf88457cb0c8c1f2283264168fbab59ebad3db9599008bfd6dde`),
> whose content was fixed in commit `e62592444170907104d8375473d40b41cc053097` and marked frozen at
> 2026-10-25T20:58:37Z (commit `1c1faf6`, tag `prereg-v2-frozen`). This translation was made after the analysis
> (2026-10-27) to let readers check the paper against the protocol; where the two differ, the Japanese file
> governs. Code identifiers, file names and numbers are unchanged. Translator's notes are in [square brackets].

Status: **FROZEN** 2026-10-25T20:58:37Z / commit e62592444170907104d8375473d40b41cc053097
(sha256 of the frozen inputs, the incidence table, N-read and the power analysis are in §10. All later changes are
recorded in `docs/protocol/amendments/` with B numbers.)

Freeze condition: the protocol is frozen when the power simulations of §8 (synthetic labels only) are finished and
**no relation between the real global labels and the place-name or terrain features has ever been computed**. Steps
allowed before the freeze are data acquisition, rasterisation, feature extraction and label-only descriptive
statistics (flooded fraction and flood-type breakdown per AOI).

Related documents: v1 (Japan, paper_1) remains unchanged as `preregistration_v1.md` with amendments A-1..A-9. The
design rationale of v2 is `docs/research/2026-10-25_v2_plan.md` (including the exploratory diagnostics of v30).

---

## 1. Questions and hypotheses

| ID | Question | Main comparison | Type |
|---|---|---|---|
| H1 | Do settlement names encode flood-proneness (controlling for position only)? | M0c → M0n | secondary (descriptive) |
| **H2** | Even after controlling for terrain, hydrology, land cover, settlement geometry and the regional composition of names, does the place-specific information of names improve flood prediction in unseen places? | **M1c → M3n** (within-place increment, §4) | **primary** |
| H2s | Does that improvement come from the meaning of the names (difference from a falsifying placebo)? | M3np → M3n | falsification condition of the primary |
| H2t | Including the regional composition, do names improve prediction as data (including proxies for language and region)? | M1s → M3s | secondary |
| H3 | Is the value of place names larger the poorer the available terrain information? | ladder T0–T4 × with/without names | secondary (confirmatory, monotonicity) |
| H4 | Does a language-model reading of place names improve the prediction of floods after the model's training cut-off? | M1rc → M3rn (N-read, llmcut split) | secondary |
| H5 | Heterogeneity of the effect (low relief, built-up areas, language areas, lexicon coverage) | strata | exploratory (descriptive only) |

Predictions stated in advance (reported even if they fail): H2 is small (R < 0.05 is likely). H3 declines
monotonically, largest at T0 and smallest at T4. H1 is positive. H4 is small (no larger than N-lex).

## 2. Data

### 2.1 Observed flood labels
- Source: the delineation (DEL) products of Copernicus EMS Rapid Mapping. The inventory is built by
  `scripts/ems_inventory.py`, and the version retrieved on 2026-10-25 is fixed as
  `docs/data/ems_inventory_20261025.json` (activation code, AOI, product URL, version). Activations published after
  that are not used (this is also the upper limit of the time hold-out). The inventory contains no Japanese
  activation (no overlap with Study 1).
- One product per AOI: excluding MONIT (monitoring), preferring wide-area products (OVERVIEW/MAP/PRODUCT), latest
  version (`topsem.world.ems.select_products`).
- Grid: LAEA centred on the AOI centroid, 250 m cells. Cover fractions are computed on a 25 m sub-grid (50 m if the
  AOI exceeds 250,000 cells) (`build_aoi_labels`).
- Evaluable cell: AOI cover ≥ 0.5, observed area (image footprint − Not Analysed) ≥ 0.5, permanent water (the
  product's reference hydrography polygons: hydrographyA / hydrography_poly / Hydrology_A_LAF) < 0.5.
- **Primary label (riverine)**: y = 1 ⇔ cover of flood polygons whose obj_desc is Riverine ≥ 0.5. y = 0 ⇔ flood
  cover of any type < 0.5. Cells flooded only by other types are excluded from the primary analysis.
- Flood types are determined from the **feature attribute obj_desc** (the activation's subCategory is not used).
  Products whose obj_desc does not state the type (e.g. only "Flooded area") follow the rule of §2.4.
- DEM-derived layers such as floodDepthA are not used.

### 2.2 Place names
- GeoNames full dump (retrieval date and sha256 recorded). Primary = feature class P (populated places).
  Sensitivity = P+L+T+S. **H (water bodies) is excluded from all analyses.**
- Strings are `asciiname`. Each cell is assigned the K = 3 nearest P names (with distances).

### 2.3 Covariates (all public versions on AWS Open Data, fixed)
- Terrain: Copernicus GLO-30 DSM, ASF GLO-30 HAND (v1/2021). Reprojected bilinearly to the sub-grid and aggregated
  to 250 m cells (mean, minimum, standard deviation, slope). 1 km mean elevation, relative depression (paper_1
  definition).
- Hydrology: HydroRIVERS v1.0 — distance to all channels and to major channels with upstream area ≥ 1000 km².
- Land cover: ESA CCI Land Cover (1992–2015, v2.0.7cds) / C3S Land Cover (2016–2020, v2.1.1), 300 m annual. **Each
  AOI receives the map of the year before the flood** (the latest, 2020, for floods from 2021), as cover fractions
  of 8 groups: built-up, cropland, tree, shrub, grassland, bare, wetland and water (`topsem.world.landcover`; groups
  defined in CCI_GROUPS in the code). Because it is the previous year's map, no image of the flood year enters.
  Mapped wetland and water also enter the baseline model, and names are tested for information beyond them.
  (The draft planned ESA WorldCover 2021, but for floods of 2012–2020 that would be post-event information, and bare
  ground or channel changes caused by floods could leak label information into the baseline, so it was replaced
  before the freeze.)
- Settlement geometry: distance (log) to the nearest P name, number of P within 1 km and 3 km.
- Timing (v2 section of `docs/data/point_in_time.md`): land cover predates the flood in every case. Channels
  (HydroSHEDS, from the 2000 SRTM) also predate every case. Terrain (GLO-30 acquisitions 2011–2015) can overlap
  floods of 2012–2015, but terrain is treated as nearly invariant to events (sensitivity analysis §7-8 [i.e. item 7
  of §7] excludes events up to 2015). GeoNames is a 2026 dump, after every case, but settlement names are treated as
  long-lived (reported as a limitation).

### 2.4 Inclusion criteria (per AOI, frozen)
EMS products come in 4 schema generations (generation 0, 2012–2013: `AreaOfInterest`, `Crisis_event_A`,
`Hydrology_A_LAF` / generation 1, to 2017: `crisis_information_poly` / generation 2: `observed_event_a` /
generation 3: `observedEventA`). Generations 0 and 1 record neither the flood type nor the image footprint. Flooding
in generation 0 is the `Crisis_event_A` polygons whose SBTYPdes is Flood (all polygons in products without that
field); `intrp`, whose code list is not published, is not used. Permanent water is `Hydrology_A_LAF` only
(`*_Optional_LAF` contains islands). Attribute names are case-insensitive. The generation is recorded as
`schema_generation` in meta.json. **Products that delineate no flood** (no flood layer; 58 of the 999 labelled AOIs)
are kept in the inventory as all-zero labels and are excluded by criterion 1 below (keeping the denominator).
1. Evaluable cells ≥ 400, cells with y = 1 ≥ 20, share of y = 1 ≤ 0.8
2. Riverine obj_desc ≥ 0.5 of the area of the AOI's flood polygons. Products whose untyped flooding (generations 0
   and 1 etc., "other") is ≥ 0.5 of the area are treated as riverine flooding if **inland** (AOI centroid more than
   10 km from the Natural Earth 10 m coastline), and the untyped flooding is also counted as y = 1. Coastal untyped,
   Flash and Coastal products are excluded from the primary analysis and go to a sensitivity analysis (§7-1)
3. At least one P name within the bounding box of the AOI grid (no margin)
4. **At most 10,000 evaluable cells in total per location group (§5).** They are allocated to the AOIs of the group
   in proportion to their evaluable cells (largest remainder) and drawn at random within each AOI with a fixed seed
   (`group_quotas`). Because the primary estimator pools cells, this gives "places" nearly equal weight
   (corresponding to the question of prediction in unseen places). The draft rule (the largest cap at which no
   single AOI exceeds 2% of all cells) gave 298,000 (effectively no cap, 14.92 million cells) on the incidence table
   of 2026-10-25 (500 included AOIs, 195 location groups), and the effective (Kish) number of location groups fell
   to 35 of 195, so it was replaced before the freeze on the basis of a label-only comparison
   (`backend/verification/r_okauchi/v32_cap_choice.py`, AOI caps 1,000–298,000 and group caps 5,000–100,000) (group
   cap 10,000: effective number 165, 1.51 million cells)

## 3. Name channels

| Channel | Definition | Memorisation | Role |
|---|---|---|---|
| N-lex | `lexicon_hydro_v1.yaml` (47 languages, 6 categories) matched to the K = 3 nearest names. 6 indicators of the nearest name + 6 distance-weighted sums over K names + wetside and coverage 2 = 14 dimensions | impossible (the lexicon knows no places) | primary |
| N-learn | Score (1 dimension) of a regularised logistic regression (character 1–4-grams of asciiname, hashing 2^20, L2, C = 1, aggregated and weighted by name) learned only on the training side of the outer fold, by cross-fitting over 5 inner folds | impossible (corpus-independent) | primary |
| N-read | Gemma-4-31B-it (HF snapshot 842da37, vLLM 0.31.0, temperature 0, at most 48 tokens, thinking mode off) judges the **literal meaning** of the elements of each P name into 9 classes (water / wetland / riverine / crossing / lowland / highland / other_environmental / non_environmental / unknown) with a confidence (`topsem.world.nread`; prompt sha256 in §10). Features are the 9-class one-hot, confidence and coverage of the nearest name and the distance-weighted 6 hydrological classes over K names = 17 dimensions. Names whose response cannot be parsed are "not judged" | possible | secondary (H4) |
| placebo | Names are randomly permuted within the set of P names of each AOI (positions and number of settlements unchanged), and N-lex, N-learn and N-read are recomputed by the same procedure (the same permutation for all 3 channels) | — | H2s, H4 |

## 4. Arms (feature sets) and the form of the main comparison

**The main comparison is a nested incremental model** (the same form as the incremental value of adding a new
marker to an existing risk score). The predictions of the base model (the HistGB arms below) are made, on the
training side, by cross-fitting over 5 inner folds (by location group, the same folds as N-learn) and, on the
evaluation side, by a model fitted on the whole training side. The second stage is an L2 logistic regression
(C = 1) on [logit of the base prediction, added features] standardised by the training-side mean and standard
deviation (`world_eval.crossfit_base` / `stage2`).

**The main comparison is the within-place increment**: the reference (M1c) contains the base prediction and the
**regional composition of names** (AOI means, c_*, of the 15 name-feature dimensions — the 14 of N-lex and the
N-learn score), and the names model (M3n) is the reference + the name features themselves. Both share the same base
prediction and the same regional composition, so their difference is "the information that names carry for each
place" (the content of the claim about dangerous place names: this place bears the name of a marsh). The regional
composition is "which names line up in the region"; through language and region it picks up regional flood rates
and works just as well with the permutation placebo (pre-freeze pilot3, §8), so it is not evidence of the meaning
of names.

- Primary arms: M1c / M3n (+ N-lex + N-learn) / M3np (+ placebo) / M3nx, M3nl (channel decomposition) / M1rc,
  M3rn, M3rnp (H4: regional composition = AOI means of the N-read features) / M0c, M0n (H1) / T0c–T3c, T0n–T3n
  (ladder H3; the top rung is M1c → M3n)
- Secondary 1 (total increment including the regional composition): M1s (recalibration of the base prediction
  only) → M3s (+ name features), placebo M3ps; for H4, M1s → M3rs / M3rps. Shows whether names improve prediction
  "as data", but includes proxies for language and region
- Secondary 2 (between flexible models): comparisons that refit the HistGB arms below with the names (M1 → M3 etc.).
  In the pre-freeze pilot2, comparisons between separately fitted flexible models were found to be swamped by model
  fluctuation (up to 29% of the terrain increment), so they were made secondary

Arms of the base and flexible models:

- M0: position (latitude and longitude of the cell centre)
- **M1**: M0 + terrain T4 (1 km elevation, mean elevation, relative depression, minimum elevation, standard
  deviation, slope, HAND mean, HAND minimum) + hydrology (log distance to all channels, log distance to major
  channels) + land cover 8 + settlement geometry 3
- **M3**: M1 + N-lex + N-learn
- M3p: M1 + placebo(N-lex) + placebo(N-learn)
- M3x / M3l: M1 + N-lex only / M1 + N-learn only (channel decomposition, secondary)
- M3r / M3rp: M1 + N-read / M1 + placebo(N-read) (H4)
- M0N: M0 + N-lex + N-learn (H1)
- Ladder (H3): T0 = M0 + hydrology + land cover + settlement geometry, T1 = T0 + 1 km elevation, T2 = T1 + mean
  elevation + relative depression, T3 = T2 + minimum elevation, standard deviation, slope, T4 = T3 + HAND (= M1).
  At each rung ± (N-lex + N-learn)

Learner: `HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1)` (v1 frozen values, otherwise defaults).
Mean of the predicted probabilities over K = 5 seeds (42–46).

## 5. Splits

- **Primary: location groups held out.** The bounding box of each AOI is buffered by 5 km, and AOIs that intersect
  across activations are joined by union-find into connected components = location groups (AOIs of different
  groups are at least 10 km apart). Groups are ordered by continent → number of cells and assigned cyclically to
  10 folds (deterministic). If the same town reappears in a different flood, it does not enter both the training
  and the evaluation side.
- Secondary 1: countries held out (10 folds by GeoNames cc country, assigned by the same rule) — transfer to unseen
  languages. Location groups that straddle a border are not split but go to the country of the majority of their
  cells
- Secondary 2: time held out — evaluation = AOIs with event time ≥ 2024-01-01, training = earlier AOIs **that do
  not belong to the same location group as an evaluation AOI** (earlier floods at the same place are removed from
  both sides).
- H4 split (llmcut): the same rule as secondary 2 with the boundary at 2025-02-01 (the first of the month after the
  pre-training data cut-off of Gemma-4, 2025-01, stated in the model card). The period of the post-training data is
  not published, so the absence of information after the cut-off cannot be guaranteed (reported as a limitation)
- Bootstrap unit: location groups for the primary, time and llmcut splits; countries for the country split

## 6. Estimand, inference and decision

- Primary estimand: the pooled held-out ΔBS(M1c → M3n) = BS(M1c) − BS(M3n) on the primary split (positive =
  improvement by names; the within-place increment of §4).
- Ratio **R = ΔBS(M1c→M3n) / ΔBŜ(M0→M1)** (the within-place increment of place names as a fraction of the
  increment of terrain etc. over position — the point estimate of the flexible models' M0 → M1). Placebo, channels,
  H1, H4 and the ladder take the same form with a common denominator.
- **Fixed-margin method**: the decision is made on the ΔBS scale with the margin fixed at SESOI × ΔBŜ(M0→M1)
  (point estimate) (the fixed-margin method of non-inferiority trials). Interval of R = interval of ΔBS /
  ΔBŜ(M0→M1). The bootstrap coefficient of variation of the terrain increment is large, 0.4–0.55 (pilot3), and the
  bootstrap interval of the ratio is inflated by the fluctuation of the yardstick, so the decision is made with the
  uncertainty of the name increment only. The bootstrap interval of the ratio is given for description.
- Inference: cluster bootstrap with location groups as the unit (B = 4000, fixed seed). R uses the same resamples
  for numerator and denominator. Robustness: two-way cluster bootstrap over activation × country (v1
  implementation).
- **SESOI: R = 0.05**.
- All intervals are percentile intervals.
- Decision (H2, `world_aggregate.decide`; the first that applies, from the top):
  1. undefined-denominator: lower end of the 95% CI of ΔBS(M0→M1) ≤ 0 (R cannot be defined; report ΔBS)
  2. not-meaning: lower end of the 95% CI of R > 0 but **lower end of the 95% CI of R − R_p ≤ 0** (H2s fails; an
     improvement that the placebo matches is attributed not to meaning but to settlement layout or
     permutation-invariant structure)
  3. positive-meaningful: lower end of the 95% CI of R > 0, H2s holds, R̂ ≥ 0.05 (practically meaningful
     improvement)
  4. positive-negligible: as above with R̂ < 0.05 and the 90% CI of R ⊂ (−0.05, +0.05) (larger than 0 but below
     the SESOI = statistically positive, practically negligible)
  5. positive-small: as above with R̂ < 0.05 and equivalence not established (not settled relative to the SESOI)
  6. negative: upper end of the 95% CI of R < 0 (adding names makes prediction worse)
  7. **equivalent-null**: none of the above and the 90% CI of R ⊂ (−0.05, +0.05) (TOST, α = 0.05)
  8. inconclusive: otherwise (report the achieved bound)
  The success of TOST is reported alongside, independently of the label.
- **Interpretation of negative (before the freeze, from the power results)**: with synthetic labels without a name
  effect (δ = 0), 7 of 100 were judged negative (the estimated weight of the names does not transfer to new places
  and slightly degrades prediction). Negative is read not as evidence that "names mislead" but as showing that
  "there is no transferable information", with the success of TOST reported alongside.
- H4: M1rc → M3rn on the llmcut split, decided by the same table (placebo M3rnp).
- Secondary: the same decision table is applied to the total increment (M1s → M3s, placebo M3ps) and to the
  flexible comparison (M1 → M3, placebo M3p) and reported alongside (the primary decision is not changed).
- H3: for the name increments of the ladder Δ_k = ΔBS(T_k c → T_k n) (k = 0..4, T4 = M1), the prior predictions
  are Δ_0 − Δ_4 > 0 (one-sided, bootstrap) and Spearman(k, Δ_k) < 0. The "recovery of terrain information" at T0,
  ΔBS(T0c → T0n) / ΔBS(T0 → T4), is reported.
- H1: ΔBS(M0c → M0n) and the R-like quantity (denominator ΔBS(M0→M1)).
- Multiplicity: the only primary is H2 (H2s is a necessary condition of the primary decision). H3 and H4 are
  secondary, each at α = 0.05 on its own, with Holm-adjusted values also given for interpretation. H5 is
  descriptive only.
- Strata of H5 (`world_aggregate.strata`, `by_continent`): whether the nearest name matches any lexicon category,
  wetside or not, mean in-cell slope at or below the median of evaluable cells (low relief) or not, CCI built-up
  cover > 10% (built-up) or not, and continent. In each stratum ΔBS(M1c → M3n) and its ratio to the within-stratum
  ΔBS(M0 → M1) are described (strata with fewer than 1,000 cells or 20 flooded cells are not reported).
- Secondary metrics: log loss, PR-AUC, capture of flooded cells in the top 10% of area (same definition as v1).

## 7. Sensitivity analyses (listed in advance; the results do not change the primary decision)

On the primary split, only the arms needed for the primary decision (M0, M1, M1c, M3n, M3np) are run, and the same
decision table is applied and reported alongside (options of `world_eval.py`).

1. All flood types (including Flash, Coastal and other; y = 1 ⇔ flood cover ≥ 0.5, without the flood-type
   condition of the inclusion criteria) — `--flood-types all`
2. Weighting of the cap: no location-group cap (all evaluable cells) — `--group-cap 0`; at most 2,000 per AOI
   (AOIs weighted equally) — `--aoi-cap 2000`
3. Place-name classes P+L+T+S for the name channels (settlement geometry stays P) — `--name-class PLTS`
4. Exclude products without an image footprint (generations 0 and 1 = footprint_missing in meta, and products with
   an implicit footprint = footprint_implicit) — `--exclude-no-footprint`
5. Learner K = 1 (seed 42 alone) — the same obligation as A-1 of v1 — `--seeds 42`
6. Location groups with boundary width 10 km → 20 km (buffer 5 → 10 km) — `--group-buffer-km 10`
7. Exclude events up to 2015 (cases overlapping the GLO-30 acquisition period) — `--min-event-year 2016`

The draft's "maximum-extent labels adding MONIT products (union per AOI)" was removed before the freeze. Monitoring
products are images of later dates, and a union mixes flood extents of different dates, changing the meaning of the
label (the extent at one acquisition). It would also require re-acquiring and relabelling products, beyond what
could be done reliably.

## 8. Power (done before the freeze, synthetic labels only)

The real features (terrain, hydrology, land cover, name assignment) and the incidence table (AOIs, groups, cell
counts, flooded fractions) are used as they are, and **only the labels are synthesised**:
y* ~ Bernoulli(σ(α_AOI + f(terrain, hydrology) + δ·g(N-lex) + spatially correlated noise)) (`world_eval.synth_labels`;
the coefficients of f are fixed in advance and not estimated from the real labels). α_AOI is matched to the real
flooded fraction of the AOI (the only use of the real labels, as an AOI-level aggregate). δ is set to the value at
which R = 0.05 by fitting R ≈ cδ² to the R of pilot runs (δ = 0.5, 1, 2, twice each). With δ = 0 and that δ, 100
replicates each (increased to 200 each before the freeze if the decision rate is within ±0.05 of the target), the
procedure of the main analysis (arms, splits, location-group cap, inference) is run in a reduced version (K = 1) to
estimate (a) the probability of declaring a null settled by TOST at δ = 0, and (b) the probability of a positive
decision at R = 0.05. The target is (a) ≥ 0.8 and (b) ≥ 0.8. If not met, the cap of §2.4 and the SESOI are revised
**before the freeze**, with the reasons added to this document. In the synthetic runs too, the arms M0, M1 (flexible,
denominator), M1c, M3n, M3np (primary) and M1s, M3s, M3ps (secondary total increment) are run, and the decision uses
the same `verdict` as the main analysis (fixed margin, placebo condition included) (bootstrap reduced to 1000).
The false-positive rate (frequency of positive at δ = 0) is also reported.

Record of the pilots (the basis of the pre-freeze design changes): in pilot2, which ran the comparison between
flexible models (M1 → M3) twice each at δ = 0.5, 1, 2, ΔBS(M1→M3) fluctuated from −1.2e-3 to +4e-4, down to its
sign, regardless of δ (about 25% of ΔBS(M0→M1) ≈ 4e-3); the placebo fluctuated as much, and the half-width of the
90% CI of R was 0.18 to 0.47. The SESOI (R = 0.05, ΔBS ≈ 2e-4) could be neither detected nor declared equivalent,
so the main comparison was changed to the nested incremental model of §4. In the nested pilot3 (δ = 0, 0.5, 1, 2,
twice each), the improvement from adding names ΔBS(M1s→M3s) was positive, +4.4 to 6.0e-4, even at δ = 0, and the
permutation placebo gave just as much (the regional composition picks up regional flood rates). The difference
between real and placebo responded cleanly to δ (±1e-4 at δ = 0, +1.5 to 3.2e-4 at δ = 2). Also, the bootstrap
coefficient of variation of the terrain increment was 0.4–0.55, and the fluctuation of the yardstick inflated the
interval of the ratio. The main comparison was therefore made the within-place increment (regional composition on
the reference side), and the decision the fixed-margin method.

Calibration of δ (pilot4, within-place increment, fixed margin, δ = 0, 1, 2, 3 twice each, cap 10,000): R was
−0.001 and −0.003 at δ = 0 (half-width of the 90% CI 0.004–0.005, both equivalent-null), 0.039 and 0.027 at δ = 1,
0.196 and 0.184 at δ = 2, and 0.529 and 0.468 at δ = 3. The placebo R was −0.006 to +0.001 in all conditions. A
log-log fit between δ = 1 and 2 (R ∝ δ^2.53) gives δ = 1.18 for R = 0.05 → **δ = 1.2** (rounded to one decimal). A
pilot with a location-group cap of 30,000 (δ = 0 and 2, once each) did not narrow the interval (half-width of the 90%
CI 0.0045), so the cap stays at 10,000. The main computation uses 100 replicates each at δ = 0 and δ = 1.2 (200 if
the decision rate is within ±0.05 of the target).

**Results of the main computation** (runs `power_v3_null`, `power_v3_sesoi`, 100 replicates each, K = 1,
bootstrap 1000): at δ = 0, TOST held in 99% (labels: equivalent-null 91, negative 7, inconclusive 1, not-meaning 1),
positive 0%, mean R −0.0008, median half-width of the 90% CI 0.0035. At δ = 1.2, positive 98%
(positive-meaningful 69, positive-small 29, inconclusive 1, undefined-denominator 1), mean R 0.0567. The targets
(each ≥ 0.8) are met, and the decision rates are more than 0.05 from the targets, so the replicates are not
increased to 200.

## 9. Research integrity

- The 3 Japanese events of v1 are development data whose results have already been seen. All Japanese diagnostics
  from v30 on are reported as exploratory and were used only as the basis of v2 design decisions
  (`2026-10-25_v2_plan.md`).
- The lexicon (N-lex) is compiled by an LLM (Claude Opus 5.5) from general knowledge only and frozen before any
  relation with the global data is seen. If the matching rules or the vocabulary are changed after the freeze, this
  is done by amendment (a new version of the lexicon is used only for secondary analyses).
- All changes after the freeze are recorded in `docs/protocol/amendments/` with B numbers, stating the analyses
  affected.

## 10. Freeze record

| Item | Value |
|---|---|
| Freeze UTC | 2026-10-25T20:58:37Z |
| commit (frozen code and documents) | `e62592444170907104d8375473d40b41cc053097` |
| lexicon_hydro_v1.yaml (version 1.0) sha256 | `6c9999921956ba13f929b18ef921afcd0e05334293ef0989d9357c6a8d59c58d` |
| ems_inventory_20261025.json sha256 | `f5c55eff3e470bce1fee7ae228641d1d3b21361f30b6d5ac9b87eef438028699` |
| GeoNames allCountries.zip sha256, retrieval | `4bc2e106b57fbc34d21f6b66e04fa4d585c74c85f6e5677fd1ec2b14937c3f43`, 2026-10-25 (the NAS copy at 09:38:35Z and the local copy at 08:56:58Z are identical) |
| Natural Earth ne_10m_coastline.zip sha256 | `bfa04cdbcbef07ef90dfca1dabb48062eca29900a113df0f389303e255484017` (version 5.0.0-pre9) |
| commits of the images that built labels, features and land cover | labels: c4672f1 (999 AOIs; the existing 957 match the content fingerprints of the a94ee9a version = v31) / features: a94ee9a (957), c4672f1 (42, same feature code) / land cover: 71611d6 (988), b11f55e (11, differing only by retries) |
| incidence table (incidence.json, sha256 `dafe37008d971d6c…`) | 500/999 AOIs included, 188 activations, 195 location groups, 58 countries, 1,507,486 cells (location-group cap 10,000, Kish effective number 165.4), 2012–2026 |
| N-read | google/gemma-4-31B-it, HF snapshot `842da3794eaa0b77d5f08bae87a17459d91ff475`, prompt sha256 `e9f55e3a65714df6f7412d96a2ac789eaa6f95436f77f77213db7fd08d2f4a8f`, input names.parquet sha256 `42861f3d647af89959481a4d42332e1a1e54c81525f4b975dc32bdab8700d197`, inference Miyabi job 3701504 (repo 9b94b6a, parse success 99.95%), output nread.parquet sha256 `fbd1dbcca2fc4711ab08087b1aefea5ccfd827e0955c3072a7f7ff2fdfc0b8dc` |
| power simulation runs | `power_v3_null` (δ = 0, 100 replicates: TOST held 99%, positive 0%), `power_v3_sesoi` (δ = 1.2, 100 replicates: positive 98%). The calibration pilots pilot2–pilot4 are in §8 |
