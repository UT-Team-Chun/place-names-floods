# B-2 Post hoc robustness analyses in response to the internal review (2026-10-26) — English translation

> English translation (2026-10-27) of `docs/protocol/amendments/B-2_posthoc_review_analyses.md`. The Japanese file is
> the record; where the two differ, the Japanese file governs. Links to internal review reports are omitted (they
> are in the development repository). Translator's notes are in [square brackets].

Status: **registered (2026-10-26, before the runs and computations below)**. The confirmatory verdicts (the labels
of H2, H4 etc. from the decision table of §6) are not changed. Everything below is **post hoc**, added in response to
the internal review held after the confirmatory results had been seen. The paper labels these analyses post hoc and
reports them separately from the confirmatory verdicts. This document fixes what is computed and how it is reported
before the results of the additional runs and computations are seen.
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924).

## Points raised by the review (summary)

1. **Inflation of the yardstick**: in the country, time and llmcut splits, the BS of the position-only flexible model
   M0 is worse than a constant prediction, and ΔBŜ(M0→M1) mixes in "the removal of M0's overfitting". The margin
   (SESOI × yardstick) may become too large
2. **Recalibration mixed into the increment of the regional composition**: M0 → M0c and T0 → T0c include the
   improvement from recalibration in the second stage
3. **Events on both sides**: in the primary split, 54 of the 188 activations span several folds
4. **A single way of splitting**: the variation due to the fold assignment is not reported
5. **No positive control**: it is not shown whether the name classes are linked to terrain (independently of the
   labels)
6. **Two-way version of the inference**: the activation × country two-way bootstrap of the decision (fixed-margin
   TOST) is missing for all splits and for H4
7. **Translation of the effect size**: it is not shown what effect of names the SESOI and the achieved bound
   correspond to

## Additional runs (image 1792984538947461000, commit eaafbd3; with default values the same behaviour as at the freeze)

| run tag | split | arms (HGB / second stage) | purpose (point) |
|---|---|---|---|
| `world_v2_recal` | all 4 splits | M0, T0, M1 / M0s, M0c, M0n, T0s, T0c, T0n, M1s, M1c, M3n, M3np, M3s, M3ps | references with recalibration only, M0s and T0s (2). Reproducibility is also checked by whether M0, M1, M1c and M3n match `world_v2` |
| `world_v2_actgrp` | primary | M0, M1 / M0s, M0c, M1s, M1c, M3n, M3np, M3s, M3ps | `--group-by-activation`: AOIs of the same activation are joined into the same location group (195 → 106 groups) (3) |
| `world_v2_fs1`, `world_v2_fs2` | primary | as above | `--fold-seed 1`, `2`: change the fold assignment while keeping the balance by size band (4) |
| `power_v3_d020`, `d040`, `d060`, `d080` (20 replicates each), `d160` (10) | primary (synthetic labels) | as in §8 (seed 42 only) | a grid of the name effect δ (in addition to the existing δ = 0 and 1.2, 100 each). Random numbers shared across δ within a replicate (7) |

## Post hoc computations (verification scripts v36, v37)

**v36 (robustness of the yardstick and the inference, points 1, 2 and 6).** For each split (and H4), the fixed-margin
TOST of the within-place increment of names ΔBS(M1c → M3n) (H4: M1rc → M3rn) is recomputed with the following
yardsticks Y (margin = 0.05 × Y):

- Y0: as registered, ΔBŜ(M0 → M1) (point estimate)
- Y1: the lower end of the 95% cluster-bootstrap interval of ΔBS(M0 → M1)
- Y2: ΔBŜ(M0s → M1s) (the increment of terrain etc. after recalibrating both)
- Y3: ΔBŜ(constant → M1s). The constant is the training-side flooded share per fold (splits with 10 folds), or the
  evaluation-side flooded share (time splits with 1 fold; favouring the constant and making the yardstick smaller =
  conservative)

Intervals both from the cluster bootstrap (as registered, B = 4000) and from the activation × country two-way
bootstrap (B = 4000). For each yardstick, whether equivalence holds and "the smallest yardstick at which equivalence
holds" (= the maximum absolute value of the 90% interval / 0.05) are reported. The increment of the regional
composition is reported as M0s → M0c, T0s → T0c and M1s → M1c (descriptive only).

**v37 (positive control and description, points 5 and 7)**:

- Positive control (independent of labels): within AOIs, the differences between cells whose nearest P name is
  water-side (`WET_SIDE` of the frozen implementation = water, wet, riverine, crossing, low; feature
  `lex_wetside_n1`) and other cells in the share with HAND ≤ 2 m, the share in the lower half of the AOI's
  elevation, and the log distance to channels, estimated with a cluster bootstrap by AOI (if names remember terrain,
  the surroundings of water-side names are low and close to rivers). As a reverse control, the highland class (high)
  with the same quantities
- ΔBS(M1c → M3n) by name class and by distance band to the nearest name (≤ 0.5 km, 0.5–1 km, > 1 km) (descriptive)
- Overlap of the top 10% predicted cells of M1c and M3n (the whole held-out primary split and within AOIs)
- Translation of the effect size: on the grid of synthetic labels (existing + the additions above), the mean
  ΔBS(M1c → M3n) per δ is obtained, and the δ corresponding to the margin (3.78e-4) and to the 90% upper bound of the
  primary split (2.10e-5) are found by monotone interpolation and reported as odds ratios e^δ (the odds ratio of
  flooding according to whether the nearest settlement name is water-side, with terrain fixed). Stated to be a
  translation under the synthetic form (§8)

## Reporting rules

- The confirmatory verdicts (labels) remain as in §6. If the post hoc analyses disagree with the confirmatory
  verdicts, the disagreement is reported as it is and the claims are narrowed to that scope (the verdicts are not
  reassigned)
- Claims, thresholds and targets are not re-chosen on the basis of the post hoc results. All the quantities above are
  reported (not only the favourable ones)

## Addendum (2026-10-26, same day, before the computations below)

Added after re-reading the full review report (points 3, 4, 9, 10, 12 and 14). Also post hoc, the confirmatory
verdicts are not changed.

- **Y4 added to the yardsticks of v36**: ΔBŜ(M0c → M1c) (a yardstick between second stages, used by the review in its
  trial calculation). The country of the two-way bootstrap is the majority country of the AOI (the same as the unit of
  the country split), B = 4000. Bootstrap intervals of the ratio are given for all R. The H5 strata add lexicon
  coverage (whether the country of the nearest name has a lexicon) and continent, and R is also given for cells
  excluding built-up areas (CCI built-up > 10%). The Holm-adjusted values of H3 and H4 are given numerically. For the
  comparisons that were not settled (H2t, flexible models), the achieved 90% upper bound and an indication of the
  number of location groups needed to exclude the SESOI (current number × (90% half-width / (0.05 − |R̂|))²) are
  given
- **v38 (origin of the labels)**: from the source table of each AOI's product, read the sensor names of the
  post-event images and classify into radar only, optical only and mixed (a one-off check reading the product zips
  on the NAS, `world_submit.py --stage script`). ΔBS(M1c → M3n) and R restricted to radar-only and optical-only AOIs
  are described (from the existing predictions, no refitting). The numbers of AOIs by schema generation and that
  entered by the rule "treat untyped as riverine" are also reported
- **v39 (interpretable effect size)**: a logistic regression with AOI fixed effects, logit p = α_AOI + β·logit(p̂_M1) +
  Σ_c γ_c·[class of nearest name = c] (p̂_M1 the held-out prediction of the primary split), fitted on all cells, giving
  the odds ratio e^γ_c per class with cluster-robust intervals by location group (sandwich estimator). An estimate of
  "with terrain fixed, within the same AOI, how many times higher are the odds of flooding of ground nearest a
  settlement with a water name". Descriptive only
- **Power of the secondary splits**: with synthetic labels (the same form as §8), δ = 0 and δ = 1.2, 20 replicates
  each, for each of the country, time and llmcut splits (tags `power_v3c_*`, `power_v3t_*`, `power_v3l_*`)
- **Early stopping of HistGB**: the learner ran with the scikit-learn defaults (early stopping on 10% validation when
  n > 10,000). The numbers of iterations of M0 and M1 in one fold of the primary split are recorded and reported

## Addendum 2 (2026-10-26, after the second round of review, before the runs and computations below)

A different referee agent checked the revised draft. Points that concern this record itself are written first,
followed by the registration of additional analyses before their runs and computations.

### Deviations from B-2 (already made)

- **Reported without registration**: analyses added while writing v37 and v39 (after the addendum at 12:27 and
  before the first commit at 12:46) — the positive control at the settlement points (cells containing the
  settlement), distance bands of the positive control, distance bands of the flood contrasts, the ≤ 250 m and
  250–500 m divisions of the distance bands of the main comparison and "within 500 m of settlements with water-side
  names", odds ratios by class × distance band and unadjusted odds ratios, and the check of v39's intervals by
  bootstrap. These are reported as **unregistered post hoc analyses** (stated so in the abstract and text too). The
  record cannot confirm which scales and divisions were decided first
- **Registered but not reported**: the increment of the main comparison by class, the odds ratios by class over all
  cells (without distance bands), two of the three cell-level measures of the positive control (lower half of the
  AOI's elevation and log channel distance; both include 0 for water-side names), the smallest yardstick at which
  equivalence holds, and the number of places needed to exclude the SESOI in the flexible comparison. All are
  reported
- **Error in describing a measure**: "HAND ≤ 2 m" of the positive control is the minimum HAND within the cell (any
  25 m pixel within 2 m of the drainage), true of 83% of all cells. It is rewritten as "the share of cells containing
  such ground" rather than "the share of such ground"

### Additional analyses (second round, registered before runs and computations)

1. **Within-place increment that does not depend on the coefficients of the regional composition**: M3d (placebo
   M3dp), which adds to the second stage only the name features (N-lex 14, N-learn) each centred within the AOI
   (feature − AOI mean), compared with the recalibration-only reference M1s. Because the AOI means are not included,
   regional flood rates are not learned. Run on the primary split, the activation-grouped split, the country, time
   and llmcut splits (runs `world_v2_dm`, `world_v2_dm_actgrp`), and described with the same decision table and the
   5 yardsticks as the main comparison
2. **Synthetic labels with the effect confined near settlements**: in the same form as §8, the name effect is put
   only on cells whose nearest settlement name is water-side and that lie within 500 m of the settlement
   (δ_local = 1, 2, 3, 10 replicates each, primary split, seed 42). Shows to how large a local effect the pooled
   within-place increment responds (tags `power_v3_loc_d100`, `d200`, `d300`)
3. **Extension of v36**: the activation-grouped split with 5 yardsticks × 2 bootstraps; H4_total and the regional
   composition of N-read (M1s → M1rc) also with all yardsticks and bootstraps; the registered two-way bootstrap
   (ratio, country of the cell, 1,000 replicates) alongside; strata by language (first language of the country); the
   number of places needed for the flexible comparison; markers for splits and strata whose yardstick is undefined
4. **Extension of v37**: continuous HAND (difference in cell-mean HAND, difference in within-AOI percentile) and the
   reference shares of each group; the positive control without settlements whose coordinates are rounded to whole
   arc-minutes (GeoNames precision, up to ±0.9 km); the increment of the main comparison by class; the flood contrasts
   by distance band with 2 weightings (with numbers of AOIs and places); contrasts within 500 m of settlements
   excluding cells with much built-up (CCI built > 10%) or tree cover (CCI tree > 50%)
5. **Extension of v39**: odds ratios by class over all cells (without distance bands), main effects of distance bands
   (whether flooding is less often recorded near settlements), and odds ratios within 500 m of settlements excluding
   cells with much built-up or tree cover
6. **Coverage of the intervals on synthetic labels (δ = 0)** (90% and 95%)
7. **Presentation of the power of the secondary splits**: the llmcut row is stated to inject through the lexicon (not
   the N-read of H4), with the distribution of decisions, Monte Carlo intervals and the mean R at δ = 1.2

### Reporting rules (additional)

- Concluding sentences (interpretation of effect sizes, summary of robustness) are written after their numbers are
  available. Sentences waiting for numbers are not written in advance
- Robustness is reported symmetrically, including the combinations that fail (yardsticks and splits for which
  equivalence does not hold are written in the text)

## Addendum 3 (2026-10-26 evening, after the third round of review, before the computations below)

A referee agent different from those of rounds 1 and 2 checked the revised draft (Major revision). No new runs are
made; the response consists of computations from existing artefacts and changes of reporting. The following is
registered before the computations.

### Deviations from B-2 (already made)

- **The registered two-way bootstrap (ratio, country of the cell, 1,000 replicates)** was registered to be given
  "alongside", but was shown only for the primary and activation-grouped splits. The country [−0.035, 0.070], time
  [−0.085, 0.116] and llmcut [−0.089, 0.086] splits had not been given. All splits are shown in the table, and those
  outside ±0.05 are written in the text
- **Registered but not reported**: the ≤ 0.5 km band of the increment of the main comparison, the smallest yardstick
  at which equivalence holds (computed but not shown in a table), the bootstrap intervals of the ratio (for all R),
  and the decision-table label of M3d. All are given in the SI
- **The translation of the upper bound of the two-way bootstrap into an effect size was not registered** (what was
  registered was the margin and the 90% upper bound of the cluster bootstrap). Stated as unregistered
- **Two post hoc analyses in the abstract** (the positive control at the settlement points, the odds ratios by class
  × distance band) were unregistered but written only as "post hoc", contrary to the rule of addendum 2 (state
  "unregistered" in the abstract and text too)

### Additional computations (registered before computation; all from existing artefacts)

1. **Translation of the effect size by fold drawing and grouping**: for the two re-assignments of folds of the primary
   split (`world_v2_fs1`, `fs2`), the sensitivity with places merged at 10 km (`world_v2_s6_buffer_b1`) and the
   activation-grouped split (`world_v2_actgrp`), the 90% upper bound of ΔBS(M1c → M3n) from the cluster bootstrap is
   translated into odds ratios with the same δ curves (uniform effect, effect only near settlements). Reported as a
   range rather than taking the value of the confirmatory folds as representative
2. **What the equivalence of the secondary splits excludes**: the margin of each split (0.05 × yardstick, ΔBS) and the
   95% upper bound of ΔBS are compared with the mean ΔBS of the synthetic labels at δ = 1.2 on that split (as a
   ratio). For the time split, whose margin is undefined, only the comparison with the 95% upper bound
3. **Recomputation of coverage on synthetic labels (δ = 0)**: the share of intervals including 0 (90% and 95%), a
   robust spread of R (1.4826 × MAD), and the SD of the standardised deviation (R − mean R)/SE. The reading "the SD
   between replicates is 2.1 times the bootstrap SE" is corrected
4. **Decision-table label of M3d** (the same table as the main comparison) and the real − permuted difference (95%
   interval from the cluster bootstrap)
5. **Odds ratios near settlements by class** (existing values of v39, unregistered) are given in the text for wetland,
   riverine and crossing names, consistent with the example in the abstract. No new combinations of classes are made
   (wetland + riverine, etc.), since they would be combinations made after seeing the values by class
6. **Correspondence of the numbers of countries** (58 in the list of activations, 57 majority countries of AOIs in the
   two-way bootstrap, 53 units of the country split, countries of the nearest settlements of cells) is explained in
   one sentence (descriptive)
7. **Source data of the figures**: the values of each panel are tabulated (generated by `make_figures.py`)

### Reporting rules (additional)

- Post hoc analyses in the abstract and text are labelled "unregistered" if they are. The registered all-cell positive
  control (water-side names as a group do not mark low ground at the cell level) is written beside the result at the
  settlement points
- The odds ratio near settlements is stated as "the value for water-side names as a group", noting that it differs by
  class. The upper end of the 95% interval is not called "the upper bound"
- The translation of the effect size is written conditionally, "under the model of the synthetic labels", with the
  range over fold drawings
- The positive control is summarised by class (riverine and crossing names mark low ground over all cells as well; the
  pooled null comes from water and lowland names)

## Addendum 4 (2026-10-26 evening, after the fifth round of review, before the computations below)

A referee agent different from those of rounds 1–4 checked the revised draft (minor revision). No new runs; the
following is registered before computation.

### Additional computations (from existing artefacts)

1. **95% intervals of the increment of the main comparison by class** (adding 95% intervals to the descriptive
   quantity registered in B-2) and the real − permuted difference by class (the same cluster bootstrap). No
   conclusion is written per class; they are described as scatter in both directions
2. **Operating characteristics against the real margin**: the 90% interval of ΔBS of each replicate of the synthetic
   labels at δ = 1.2 (the interval of R × the synthetic yardstick of that replicate) is compared with the real margin
   of the same split (0.05 × yardstick, ΔBS), giving the share judged equivalent (primary, country and llmcut splits;
   the time split's margin is undefined), with Monte Carlo 95% intervals (Wilson)
3. **An interval including the variation of the fold assignment**: a 90% interval combining the standard deviation of
   R over the confirmatory folds and the two re-assignments with the bootstrap standard error of the confirmatory
   analysis (R̂ ± 1.645 √(SE² + SD²)). Reported as description; the verdict is not changed

### Reporting rules (additional)

- The results by class are written for all 6 classes on the ΔBS scale, with 90% and 95% intervals. The differences in
  the yardstick between classes (up to about 17 times) are noted
- The subject of the Discussion is limited to the prediction that the confirmatory test can falsify (whether names
  improve prediction after terrain is known, within mapped areas). The question about the surroundings of settlements
  is written separately, as unregistered evidence
