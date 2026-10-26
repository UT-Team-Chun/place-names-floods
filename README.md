# place-names-floods

Analysis code, frozen protocol and lexicon for the study
**"Do place names remember floods? A pre-specified global test against satellite-observed inundation"**
(R. Okauchi and P.-J. Chun, The University of Tokyo).

The study tests whether settlement names improve the prediction of satellite-observed river flooding
(Copernicus Emergency Management Service delineations; 500 areas of interest, 188 activations, 58
countries, 2012–2026) beyond terrain, hydrology, land cover and settlement layout, in places a model has
never seen. The protocol was frozen before any model was fitted to the observed flood labels; see
[`PROVENANCE.md`](PROVENANCE.md) for the freeze record and how to verify it.

## Layout

| Path | Contents |
|---|---|
| `docs/protocol/preregistration_v2.md` | Frozen analysis protocol (hypotheses, estimand, splits, decision rules, smallest effect size of interest, sensitivity analyses, power analysis); in Japanese |
| `docs/protocol/v2/lexicon_hydro_v1.yaml` | Frozen multilingual lexicon of water, wetland, riverine, crossing, lowland and highland morphemes |
| `docs/protocol/amendments/B-1_*.md` | The single post-freeze amendment (correction of an implementation defect in one sensitivity analysis) |
| `docs/research/2026-10-26_v2_interpretation_plan.md` | Interpretation plan written after the freeze and before any confirmatory result |
| `docs/data/` | Copernicus EMS product inventory (frozen snapshot) and the point-in-time register of all inputs |
| `backend/packages/topsem/` | Analysis package (`topsem.world`: EMS parsing, rasterisation, terrain, hydrology, land cover, names, lexicon and language-model channels) and scripts (`scripts/world_*.py`, `nread_*.py`, `sim_power_v3.py`) with unit tests |
| `backend/verification/r_okauchi/v31`–`v35` | One-off verification scripts (label identity, cell cap choice, freeze record, channel agreement, exploratory name-class contrast) |
| `docs/verification/` | Their outputs |
| `docs/paper/paper_2_global/figures/` | `make_figures.py`, which generates every number, table and figure of the paper from the run artefacts, and its outputs |
| `docs/paper/paper_2_global/submission/` | Manifest (with sha256) of the derived data archived on Zenodo |

## Environment

Python 3.12 with [uv](https://docs.astral.sh/uv/). From `backend/`:

```bash
uv sync --frozen --all-extras      # the locked environment used for the analysis
uv run pytest packages/topsem/tests -q
```

The confirmatory runs used a Linux amd64 container built from this code; results on other platforms
may differ in the last floating-point digits.

## Reproducing the analysis

Inputs are public (see Data availability in the paper). The pipeline stages are scripts in
`backend/packages/topsem/scripts/`, reading and writing under `TOPSEM_DATA_DIR` (default `./data`):

1. `ems_inventory.py` → `world_labels.py` — inventory and rasterise the EMS delineations (250 m cells)
2. `world_features.py`, `world_landcover.py`, `build_rivers_global.py` — terrain (GLO-30, HAND), land cover
   (ESA CCI / C3S, year before each event), rivers (HydroRIVERS), place names (GeoNames)
3. `world_incidence.py` — label-only incidence table and inclusion (the only label use before the freeze)
4. `nread_names.py`, `nread_infer.py` — language-model reading of names (Gemma-4-31B-it, pinned revision, vLLM)
5. `sim_power_v3.py` — power analysis on synthetic labels
6. `world_eval.py --tag <run> --split {primary,country,temporal,llmcut} --fold <k>` — one held-out fold per
   call (10 folds; one for the time and language-model splits); sensitivity options are command-line flags
7. `world_aggregate.py --tag <run> --split <split>` — cluster bootstrap, decision table, diagnostics
8. `docs/paper/paper_2_global/figures/make_figures.py` — all numbers, tables and figures

The derived data (cell tables, name readings, held-out predictions and run summaries) are archived on
Zenodo; `docs/paper/paper_2_global/submission/zenodo_data_manifest.json` lists every file with its sha256.

## Licence

MIT (see `LICENSE`). Input datasets keep their own licences (Copernicus EMS, Copernicus DEM, ESA CCI /
C3S, HydroSHEDS, GeoNames CC BY 4.0, Natural Earth public domain).
