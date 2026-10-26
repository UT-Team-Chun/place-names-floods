# B-3 Check with independent flood labels (Copernicus GFM) (2026-10-26) — English translation

> English translation (2026-10-27) of `docs/protocol/amendments/B-3_gfm_independent_labels.md`. The Japanese file is
> the record; where the two differ, the Japanese file governs. Links to internal review reports are omitted.
> Translator's notes are in [square brackets].

Status: **registered (2026-10-26 evening, before the labels below were made and computed)**. The confirmatory
verdicts (the labels of H2, H4 etc. from the decision table of §6) are not changed. Everything below is **post
hoc**, added after the confirmatory results and the post hoc analyses of B-2 had been seen, in response to the
internal review (round 5, proposed change 3) and to the author's decision ("if we should do the ground near
settlements with independent flood data, let's do it"). The paper labels these analyses post hoc and reports them
separately from the confirmatory verdicts. This document fixes, before the labels are made and computed, what is
made, how it is computed and how it is reported.
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924).

## Purpose

All results, including the confirmatory analysis, rest on the Copernicus EMS flood maps (labels) alone. Even with
terrain fixed, the EMS labels mark cells within 500 m of a settlement as flooded less often (odds ratio > 1 km vs
≤ 500 m 1.29, v39 of B-2). If this is due to EMS omissions (missed flooding near buildings), and if the errors near
settlements differ between settlements with water-side names and others, the increment of names could be biased in
either direction. With an independently produced global satellite flood map we check three points.

1. **Check of the labels**: are the EMS labels thinner near settlements than GFM? Does the disagreement between the
   labels differ around settlements with water-side names and around others?
2. **The settlement-scale question**: within 500 m of settlements, are the odds ratios of flooding by name class
   (terrain fixed) the same on the independent labels?
3. **Reproduction of the main comparison**: is the within-place increment of names R (M1c → M3n) inside the margin
   ±0.05 on the independent labels, (a) when the existing held-out predictions are rescored, and (b) when the models
   are refitted on the same cells and folds?

## Independent labels: Copernicus Global Flood Monitoring (GFM)

- The global flood monitoring of CEMS. Every Sentinel-1 (IW GRD) acquisition is classified for flooding by three
  independent methods (LIST, DLR, TUW) combined by majority vote (ensemble observed flood extent, 20 m, Equi7Grid).
  An archive of all acquisitions since 2015 exists. Product definition: GFM Product Definition Document v4.0
  (2025-10-14)
- Access: the public STAC of EODC (`https://stac.eodc.eu/api/v1/collections/GFM`, no authentication) and the
  cloud-optimised GeoTIFFs listed there (`https://data.eodc.eu/collections/...`). Assets used: `ensemble_flood_extent`
  (0 = no flood, 1 = flood, 255 = no data), `exclusion_mask` (0 = classifiable, non-zero = excluded, 255 = no data),
  `reference_water_mask` (0 = normally no water, 1 = permanent water, 2 = seasonal water of that month, 255 = no
  data)
- The exclusion mask covers built-up areas (GHSL built-up probability above 35%), dense vegetation, water look-alike
  surfaces (such as sand), terrain distortion and HAND of 10 m or more, radar shadow, and insufficient Sentinel-1
  observation. GFM can compare only **land that is not excluded** (near settlements, the surrounding land outside
  built-up areas). This is stated as the scope of the analysis
- Degree of independence: in areas where the EMS post-event image is Sentinel-1, the same image classified by a
  different process (automatic majority vote vs. human-assisted interpretation). In areas where the EMS post-event
  image comes from another satellite (COSMO-SkyMed, RADARSAT-2, optical, etc.), a different observation at nearly
  the same time

## What was done before registration (no values bearing on results were seen)

- A feasibility check by catalogue search only (script in the scratchpad, record in this document): of the 500
  analysed AOIs, 348 record the date of the post-event image of the EMS product (schema generations 2 and 3). For
  317 of them there is a GFM acquisition within one day before or after the post-event image date, and for 313 it
  covers at least 90% of the bounding box of the AOI (343 and 340 within three days). Of the 313, 153 include
  Sentinel-1 among the EMS post-event images (149 of them have a GFM acquisition on the same day), 109 have other
  SAR only and 51 optical etc. Catalogue attributes on flood amount such as `flooded_pixels` were not used
- To confirm the sign of the values, only the pixel counts per value of a reduced version of one tile (E048N012T3,
  2023-05-24, northern Italy; including the extent of the analysed EMSR664) were read. Cell-level values, comparisons
  with the EMS labels and relations with names were not looked at
- No GFM label had been made for any AOI

## How the labels are made (the same rules for every AOI)

1. **Target AOIs**: of the 500 analysed AOIs, those whose EMS product records the date of the post-event image (348)
2. **Acquisitions**: all GFM items acquired from the earliest post-event image date − 1 day to the latest + 1 day
   (UTC dates) that intersect the bounding box of the AOI. Items with the same acquisition time (adjacent tiles) are
   grouped as one acquisition
3. **Grid**: the three GFM layers are mapped by nearest neighbour onto the 25 m sub-grid of the AOI's 250 m grid
   (LAEA, the same as the EMS labels; 50 m for AOIs above 250,000 cells, the same rule as the EMS rasterisation). A
   sub-grid pixel is **assessable** when the flood extent is 0 or 1, the exclusion mask is 0 and the reference water
   is 0
4. **Per acquisition**: for each 250 m cell, the share v of assessable sub-grid pixels and the share f of flooded
   pixels among the assessable ones
5. **Cell label**: a cell is assessable by GFM if v ≥ 0.5 on any acquisition; it is flooded (1) if f ≥ 0.5 on any
   acquisition with v ≥ 0.5, otherwise not flooded (0)
6. **Inclusion of AOIs**: at least 400 of the cells evaluable on EMS (the definition of evaluable cells of the
   confirmatory analysis) are also assessable on GFM, at least 20 of them are flooded on GFM, and the GFM flooded
   share is at most 0.8 (the same thresholds as the inclusion criteria of the confirmatory analysis; the riverine-
   share criterion is already met on EMS)
7. **Analysed cells**: the cells of the primary-split sample of the confirmatory analysis (drawn with the cap of
   10,000 cells per location group, world_v2) that lie in included AOIs and are assessable on GFM. Folds and location
   groups stay as in the confirmatory analysis

The labels are made with `topsem.world.gfm` (library) and `scripts/world_gfm.py` (execution), recording the items
used (ID, acquisition time, processing version, URL) for each AOI.

## Computation (verification script v41; refitting with world_eval)

The bootstrap is, as in the confirmatory analysis, a cluster bootstrap with location groups as the unit (B = 4,000,
seed 20261027). Within-AOI contrasts are, as in v35 and v37, averages over AOIs with at least 20 cells in each group,
weighted by the number of cells of the smaller group.

**A. Check of the labels** (analysed cells)
- A0 (descriptive): flooded shares of EMS and GFM, agreement, Cohen's κ, P(GFM = 1 | EMS = 0) and
  P(EMS = 1 | GFM = 0). Overall and by the kind of EMS post-event image (includes Sentinel-1 / other SAR only /
  optical etc.)
- A1: within AOIs, the difference in [EMS flooded share − GFM flooded share] between cells within 500 m of the
  nearest settlement and cells beyond 1 km
- A2: within AOIs, the difference in [EMS flooded share − GFM flooded share] between cells whose nearest settlement
  name is water-side and other cells (within 500 m and all cells)
- A3: in the same within-AOI logistic regression as v39 (AOI intercepts, distance-band indicators, logit of the
  held-out prediction of the terrain model M1), the odds ratio > 1 km vs ≤ 500 m, on the GFM labels and on the EMS
  labels of the same cells

**B. The settlement scale** (analysed cells within 500 m of the nearest settlement, and all cells)
- With the same regression as v39 (AOI intercepts, name-class indicators, and in the adjusted model the logit of the
  held-out prediction of M1; standard errors clustered by location group with the G/(G − 1) correction), the odds
  ratios of the water-side group and of the 6 classes, on the GFM labels and on the EMS labels of the same cells. AOIs
  without any flood carry no information and are removed (as in v39). Holm-corrected p values of the classes are also
  given
- The M1 used for adjustment is the held-out prediction of the confirmatory analysis (trained on EMS). As a
  sensitivity, also the M1 trained on GFM in C2

**C. Reproduction of the main comparison**
- C1 (rescoring): the held-out predictions of the world_v2 primary split (M0, M1, M1c, M3n, M3np) are scored against
  the GFM labels on the analysed cells. ΔBS(M1c → M3n), yardstick ΔBS(M0 → M1), R, 90% and 95% intervals, fixed margin
  (0.05 × the point estimate of the yardstick on GFM), difference from the placebo, and the label of the decision
  table of §6 (`world_aggregate.decide`). Also with the EMS labels of the same cells
- C2 (refitting): `world_eval --labels gfm --tag world_v2_gfm --split primary` (default arms, 5 seeds, 10 folds;
  after building the same table and folds as the confirmatory analysis, restrict to the analysed cells and replace y
  by the GFM labels) and `world_aggregate --tag world_v2_gfm --split primary`. The §6 label of the main comparison and
  R with its intervals
- Sensitivity (B and C1): only the AOIs whose EMS post-event images do not include Sentinel-1 (a different
  observation)

## Reporting rules

- All the quantities above are reported (a new SI Note and table). The text updates the results paragraph (marked
  post hoc and registered) and the limitations of the Discussion (label errors), and the abstract mentions the C2
  result in one phrase, whatever its direction
- Interpretation rules (fixed before computation):
  - If the 95% interval of A1 lies below 0, write "the EMS labels are thinner than GFM near settlements"; if it
    includes 0, write "the thinness near settlements is not specific to EMS (the same on GFM)". A3 is presented in the
    same direction
  - If the 95% interval of A2 excludes 0, write "the disagreement between the labels differs around settlements with
    water-side names" and state in which direction the increment of names could be biased. If it includes 0, write
    "no difference by name was detected" and give the width of the interval
  - B: if the 95% interval of the adjusted odds ratio of the water-side group on the GFM labels lies above 1, write
    "on the independent labels, with terrain fixed, the surroundings of settlements with water-side names flooded
    more"; below 1, the opposite; if it includes 1, write "no association was detected on the independent labels
    either" and give the upper bound. The values by class are all listed as description, and conclusions by class
    are limited to classes whose Holm-corrected 95% intervals exclude 1. The settlement-scale question is moved from
    "unresolved" only if the interval of the group odds ratio excludes 1
  - The label of C2 is written as the verdict of the reproduction on independent labels. C1 is auxiliary. If either
    disagrees with the confirmatory label, this is written as it is in the abstract and Discussion (the confirmatory
    label is not reassigned)
- Because of GFM's exclusion mask (built-up areas, dense vegetation, HAND of 10 m or more, etc.), the comparison is
  limited to land that is not excluded; this is written in the same paragraph as the results
- After seeing results, the window (± 1 day), thresholds, inclusion criteria, target AOIs and order of analysis are
  not changed. If a change becomes necessary, an addendum is added to this document and the results before the
  change are also reported

## Addendum 1 (2026-10-26 evening, while the labels were being made and before any relation with names was computed)

During the making of the labels (cluster, 6 shards), for one of the 14 items in the window of EMSR265_AOI24
(`ENSEMBLE_FLOOD_20180202T180221_VV_NA020M_E036N075T3`, a tile of the North American continental grid), the reference
water file returned HTTP 403 (file absent), and the whole AOI failed. The rules above did not cover this case, so the
following is added. Relations with names and comparisons with the EMS labels had not yet been computed for any AOI
(only the success of the jobs and their run times had been seen).

- **Items with an unreadable layer are not used**: an item for which any of the three layers cannot be read (HTTP
  errors, etc.) is removed from its acquisition, and its ID is recorded in the provenance (`skipped_unreadable` in
  `gfm_meta.json`). The other items of the acquisition are used as they are. The number of removed items and the
  number of AOIs that lost an acquisition through this are reported
- The AOIs that succeeded in the first run had no unreadable items (otherwise they would have failed), so adding the
  rule does not change their results. Only the failed AOIs are rebuilt with the code that includes this rule
- The output was changed so that an interruption leaves no partial table: the table replaces a temporary file after
  the provenance has been written (same behaviour)
- Supplement (same evening, found with EMSR861 and EMSR871): items that lack one of the three assets in the catalogue
  (e.g. `ENSEMBLE_FLOOD_20260130T061030_VV_AF020M_E030N090T3` has no exclusion mask) are also removed and recorded,
  like unreadable items
- Supplement 2 (same evening, found in the rebuilding run): removing items without an asset worked, but the job
  stopped where the item URLs are written to the provenance, which accessed the same asset (5 AOIs including
  EMSR265). The provenance now writes None for an absent layer (the label rules are the same). Only the stopped AOIs
  are rebuilt

## Addendum 2 (early 2026-10-27, a check added after seeing the result of B; post hoc)

The result of B (water-side group within 500 m, adjusted for terrain: GFM OR 0.55 [0.37, 0.82], EMS on the same cells
0.79 [0.58, 1.07]) reads, by the registered interpretation rule, "on the independent labels, the surroundings of
settlements with water-side names flooded less than the surroundings of other settlements". Because this is opposite
to the hypothesis, we check one point only: whether it can be explained by the way GFM classifies. It is stated to be
**a check added after seeing the result of B**, and the registered interpretation of B is not changed.

- **Partial masking**: GFM removes permanent-water and excluded pixels within a cell from classification. Settlements
  with water-side names lie near water, so the cells around them are more likely to be partly masked, and the
  remaining pixels may be biased towards the dry side. Restricting to cells whose GFM-assessable share v (maximum over
  acquisitions) is at least 0.9, and to cells with v = 1 (all pixels assessable), the same regression as B (within
  500 m and all cells, GFM and EMS) is computed (v45)
- Reporting rule: if, with v ≥ 0.9, the 95% interval of the GFM group odds ratio still lies below 1, write "it does
  not depend on partially masked cells"; if it comes to include 1, write "it weakens when partially masked cells are
  removed, so masking may have contributed". All values are listed
