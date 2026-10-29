# place-names-floods

Analysis code, frozen protocol and lexicon for the study
**"Landscape terms in settlement names, assigned to the nearest 250 m cells, do not improve within-area predictions of
mapped river flood extent: a pre-specified test in 500 Copernicus EMS areas"** (R. Okauchi and P.-J. Chun, The University
of Tokyo).

The pre-specified test assigns each name reading to all ground nearest its settlement and scores how well the names place
flooding within mapped areas, adjusting for which names occur in each area. An effect confined to the ground near
settlements could not have exceeded its margin, which the authors found only after the confirmatory result; the claim
that names record the flooding of settlements is not tested here.

The study tests whether settlement names improve the prediction of satellite-observed river flooding
(Copernicus Emergency Management Service delineations; 500 areas of interest, 188 activations, 58
countries, 2012–2026) beyond terrain, hydrology, land cover and settlement layout, in places a model has
never seen. The protocol was frozen before any model was fitted to the cell-level flood labels; see
[`PROVENANCE.md`](PROVENANCE.md) for the freeze record, what it can and cannot show, and how to verify it.
Everything outside the frozen protocol is exploratory; most exploratory analyses were fixed in a dated amendment
before they were computed, and the paper lists the exceptions and the planning status of every group of estimates.

## Layout

| Path | Contents |
|---|---|
| `docs/protocol/preregistration_v2.md` | Frozen analysis protocol (hypotheses, estimand, splits, decision rules, smallest effect size of interest, sensitivity analyses, power analysis); in Japanese, the frozen record |
| `docs/protocol/en/` | English translations of the protocol and of amendments B-1 to B-11, drafted with a language model after the analysis and not yet checked by the authors (the Japanese files govern where they differ) |
| `docs/protocol/v2/lexicon_hydro_v1.yaml` | Frozen multilingual lexicon of water, wetland, riverine, crossing, lowland and highland word elements |
| `docs/protocol/amendments/B-*.md` | Post-freeze amendments: corrections of implementation errors with both runs reported (B-1: one sensitivity analysis; B-4: the lexicon was not applied to Norwegian names because the YAML parser read the country code `NO` as false); post hoc analyses fixed before they were run (B-2: after critiques by language-model agents with access to the project records; B-3: independent flood labels from Copernicus Global Flood Monitoring; B-5 to B-11: after seven rounds of critiques by language-model agents given only the manuscript and supplement). These critiques are not peer review; the paper (Appendix A) lists the analyses that were planned or withdrawn after related results had been seen |
| `docs/research/2026-10-26_v2_interpretation_plan.md` | Interpretation plan written after the freeze and before any confirmatory result |
| `docs/data/` | Copernicus EMS product inventory (frozen snapshot) and the point-in-time register of all inputs |
| `backend/packages/topsem/` | Analysis package (`topsem.world`: EMS parsing, rasterisation, terrain, hydrology, land cover, names, lexicon and language-model channels, GFM labels, the post hoc modules `b6`–`b10`) and scripts (`scripts/world_*.py`, `nread_*.py`) with unit tests |
| `backend/packages/gateways/` | An unused template package of the original project scaffold, kept only because the locked workspace lists it; the analysis does not import it |
| `CLAUDE.md` | The analysis rules (pre-specification, leakage, run contract, traceability); the scripts also use this file to find the repository root, and the coding assistant used during development reads it |
| `backend/verification/r_okauchi/v31`–`v68` | One-off verification scripts: label identity, cell cap, freeze record, channel agreement, exploratory contrasts (v31–v35); post hoc analyses of B-2 (v36–v40); GFM labels and their checks (v41, v45; B-3); B-5 analyses from existing artefacts (v42), etymology checks against Wikidata (v43) and Wikipedia (v44), the size of the B-4 defect (v46) and the aggregation of the B-3 to B-5 cluster runs (v47); B-6 analyses from existing predictions (v48: the ceiling of a settlement-scale effect, AOI weights, rounded coordinates, interval coverage), surface water around settlements (v49), diagnostics of the GFM labels and the variant with seasonal water (v50), documentation and strict reading of the lexicon (v51), aggregation of the B-6 cluster runs including the plasmode bounds (v52), the surface-water positive control (v53) and the descriptions of B-6 addendum 1 (v54); B-7 analyses from existing artefacts (v55: floodplain reading, partly flooded cells, event magnitude, built-up and tree adjustment, GFM discrimination, pairs bootstrap, joint tests, digraphs, recall by language and further descriptions) and the aggregation of the B-7 cluster runs (v56: permutations with recomputed composition, plasmode with refitted models, the activation split with the corrected code, the recurrent-flooding control, the floodplain reading and the run without partly flooded cells); B-8 analyses
from existing artefacts (v57: achieved bound, groups needed, permutation test, the two label sets on the same cells,
estimates near settlements by class and variant, heights above the main rivers, calibration and skill, audit of
matched forms, a reading by river and lake names), a second coder of the etymology excerpts (v57b), the cleaned lexicon
(`v57_lexicon_clean`), surface water of each cell (v58), the aggregation of the B-8 runs (v59: bounds with Monte Carlo
intervals, the second stage refitted in the bootstrap, operating characteristics at the margin, surface water of each
cell, the cleaned lexicon), where the flooded cells near water-side names lie (v60), estimates without Pakistan and
Bangladesh (v61) and the precision of the strict reading by language (v62); B-9 analyses from existing artefacts (v63:
Murphy decomposition, within-AOI evaluation, flood-map skill, activations left out, Benjamini–Hochberg values) and the
aggregation of the B-9 runs (v64: bounds read three ways, protective and heterogeneous injected effects, simulated
power with refitted models, one permanent-water definition, subsampling seeds, the bare-earth terrain model), the
lowest rung of the terrain ladder (v65); B-10 analyses from existing artefacts (v66: the islands-and-meadows group at
settlement sites, the floodplain reading by country and form, detectable effects and three kinds of interval near
settlements, misreading over a grid of rates, validation by country, sensors and product generations, the
activation-merged split, the within-AOI evaluation on GFM labels with the placebo of addendum 2, the
near-settlement estimator on injected labels,
one definition of the HAND-only map, the between-area relation within language areas, cells read in another
language, the class distribution of the language-model reading) and the aggregation of the B-10 runs (v67: distance
and continent interactions, a tuned penalty, country indicators, spatially clustered injection, JRC flood hazard maps
as part of the baseline); B-11 analyses (v68: the change of the activation-merged split between code versions,
odds ratios near settlements on every evaluable cell with both label sets, bounds on the estimator's scale, the GFM
domain, the denominator of the labels, corrected and combined intervals, a census of exploratory estimates, feature codes
and rounded coordinates, the pilots, the toponymic audit sample, digraphs, and the aggregation of the country-split
permutations with rebuilt composition) |
| `docs/verification/` | Their outputs, including the coded Wikipedia etymology excerpts with sources and quotations (`v44_etym_codes.csv`) and the lexicon tables (`v51_lexicon_table.csv`, `v51_top_forms.csv`, `v51_examples.csv`) |
| `docs/verification/v51_audit_sample.csv`, `v51_audit_guide.md` | A stratified sample of 100 settlement names and a coding guide for an independent coding of their recorded etymologies (amendment B-6; superseded by the sample below) |
| `docs/verification/v68_audit_sample.csv`, `v68_audit_guide.md` | A sample of 150 settlement names stratified by matched form and country, weighted by each form's contribution to the results, and a guide for toponymists that asks for standard dictionaries, homographs and whether an element names the site or a nearby feature; the study's readings of the sample are withheld until the check is complete |
| `docs/verification/v68_estimates.csv` | Machine-readable census of the exploratory interval estimates (amendment B-11, item B8): source file, JSON path, kind, estimate, 95 % interval, whether it excludes the no-effect value, approximate two-sided p, amendment, family (name effect on flooding, site description or nuisance coefficient) and the rules applied; the planning status of each group of estimates is in the Supplement |
| `docs/paper/paper_2_nhess/figures/` | `make_figures.py`, which generates every number, table and figure of the paper from the run artefacts, and its outputs |

## Environment

Python 3.12 with [uv](https://docs.astral.sh/uv/). From `backend/`:

```bash
uv sync --frozen --all-extras      # the locked environment used for the analysis
uv run pytest packages/topsem/tests -q
```

The confirmatory runs used a Linux amd64 container built from this code. Results on other platforms can
differ: on an arm64 machine the evaluations reproduced the primary comparison to the fifth decimal of R, but the
simulations that refit the character model and the second stage on injected labels differed from the cluster's by
up to 1.4 Monte Carlo standard errors, whereas a rerun on the cluster reproduced them (paper, Appendix A). Run the
simulations in the amd64 container.

## Reproducing the analysis

Inputs are public (see Code and data availability in the paper). The pipeline stages are scripts in
`backend/packages/topsem/scripts/`, reading and writing under `TOPSEM_DATA_DIR` (default `./data`):

1. `ems_inventory.py` → `world_labels.py` — inventory and rasterise the EMS delineations (250 m cells)
2. `world_features.py`, `world_landcover.py`, `build_rivers_global.py` — terrain (GLO-30, HAND), land cover
   (ESA CCI / C3S, year before each event, 2020 for later events), rivers (HydroRIVERS), place names (GeoNames)
3. `world_incidence.py` — label-only incidence table and inclusion (the only label use before the freeze)
4. `nread_names.py`, `nread_infer.py` — language-model reading of names (Gemma-4-31B-it, pinned revision, vLLM)
5. `world_eval.py --synthetic {null,sesoi} --delta <δ> --reps a:b` — power analysis on synthetic labels
   (aggregated by `world_aggregate.py --tag <run> --power`)
6. `world_eval.py --tag <run> --split {primary,country,temporal,llmcut} --fold <k>` — one held-out fold per
   call (10 folds; one for the time and language-model splits); sensitivity options are command-line flags
7. `world_aggregate.py --tag <run> --split <split>` — cluster bootstrap, decision table, diagnostics
8. Post hoc: `world_gfm.py` builds the GFM labels (B-3) and `world_eval.py --labels gfm` refits on them;
   `world_eval.py` options `--perms`, `--save-coefs`, `--oracle`, `--label-variant`, `--fold-seed` and
   `--group-by-activation` serve the analyses of B-2 and B-5, and `--plasmode`, `--oracle-label`, `--gsw`,
   `--strict-names`, `--gn-filter`, `--exclude-no-hydro` and `--exclude-untyped` those of B-6 (`v49` builds the
   surface-water table that `--gsw` reads), and `--perms-recomp` with `--perm-first`, `--plasmode-refit`,
   `--recurrence`, `--floodplain` and `--exclude-partial` those of B-7; `--plasmode-deltas`, `--exclude-gsw-perm`,
   `--table-seed` and `--terrain fabdem` (with `world_fabdem.py`) those of B-9; the B-10 stacks of `--stacks`,
   `--plasmode-forms AC` with `--plasmode-rho`, and `--hazard jrc` (with `world_jrc.py`) those of B-10; the scripts
   `v36`–`v67` compute the rest
9. `docs/paper/paper_2_nhess/figures/make_figures.py` — all numbers, tables and figures

The derived data (cell tables, name readings, held-out predictions, run summaries and the GeoNames snapshot used)
will be archived on Zenodo before publication, with a manifest that lists every file with its sha256.

## Licence

MIT (see `LICENSE`). Input datasets keep their own licences (Copernicus EMS and GFM, Copernicus DEM, ESA CCI /
C3S, HydroSHEDS, JRC Global Surface Water, JRC global river flood hazard maps ODbL 1.0, FABDEM CC BY-NC-SA 4.0,
GeoNames CC BY 4.0, Natural Earth public domain, Wikidata CC0). The
quotations in `v44_etym_codes.csv` and `v51_examples.csv` are excerpts of Wikipedia articles (CC BY-SA 4.0), each
attributed by the URL of its article in the same row; they are shared under that licence.
