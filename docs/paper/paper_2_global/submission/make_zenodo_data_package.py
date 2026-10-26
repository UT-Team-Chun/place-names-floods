# %% paper_2 の派生データを Zenodo の Web 画面から上げられる形（zip 数個 + 台帳 + 説明）にまとめる
"""手元で zip を作るだけで、どこにも送らない。

出力（既定: リポジトリの外 ../zenodo_upload/place-names-floods-data-v1/）:
- cell_tables.zip   data/global/labels/<AOI>/{cells,features,landcover}.parquet と *meta.json
- runs.zip          data/runs/world_v2/、感度 world_v2_s*/、事後（B-2）の world_v2_*（集計・予測・来歴）、検出力の集計
- records.zip       発生表・言語モデルの読み（N-read）・EMS 製品台帳・検証結果（v31–v40）と v37・v38 の中間表
- MANIFEST.json     zip ごとの中身（パス・大きさ・sha256）と zip 自体の sha256
- README.md         内容・コードとの対応・入力データのライセンスと帰属表示
Zenodo は 1 レコードのファイル数に上限があるので、個々のファイルは zip にまとめる。
"""

from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT.parent / "zenodo_upload" / "place-names-floods-data-v1"
STORED = {".parquet", ".npz", ".zip"}  # 既に圧縮済みの形式は再圧縮しない

ARCHIVES = {
    "cell_tables.zip": ["data/global/labels/*/*.parquet", "data/global/labels/*/*meta.json"],
    # 確証（world_v2）・感度（world_v2_s*）・事後（world_v2_recal・actgrp・fs1・fs2・dm・dm_actgrp、amendment B-2）
    "runs.zip": ["data/runs/world_v2*/*.json", "data/runs/world_v2*/*.png",
                 "data/runs/world_v2*/predictions/*/*",
                 "data/runs/power_v3*/power_summary.json"],
    "records.zip": ["data/global/incidence.json", "data/global/nread/*",
                    "data/global/geonames/aoi_subset.parquet",
                    "docs/data/ems_inventory_20261025.json",
                    "docs/verification/v3[1-9]_*.json", "docs/verification/v40_*.json",
                    "data/verification/v38_ems_sensors_raw.json",
                    "data/verification/v37_table.parquet"],
}

README = """# place-names-floods — derived data

Derived data for the study "Do place names help predict where floods go? A pre-specified global test against
satellite-observed inundation" (R. Okauchi and P.-J. Chun). Code, frozen protocol and lexicon:
https://github.com/UT-Team-Chun/place-names-floods (archived on Zenodo as a separate record).
Paths inside the archives are relative to that repository's root, so unpacking the archives there
reproduces the layout the scripts expect (`TOPSEM_DATA_DIR=./data`).

## Contents

- `cell_tables.zip` — one directory per Copernicus EMS area of interest (AOI), one row per 250 m cell:
  `cells.parquet` — grid indices and coordinates (`ix`, `iy`, `x`, `y`, `lon`, `lat`), cover fractions of
  the AOI (`f_aoi`), of the analysed image area (`f_observed`) and of permanent water (`f_perm_water`), and
  flooded fractions by flood type (`f_riverine`, `f_flash`, `f_coastal`, `f_other`, `f_trace`,
  `f_flood_any`); `features.parquet` — terrain statistics (`elev_*`, `slope_mean`, `hand_*`, `depress`,
  `dem_valid_frac`), distances to rivers (`dist_river_all`, `dist_river_major`), the three nearest populated
  places with distances and the number within 1 and 3 km (`P_*`), the same for GeoNames feature classes
  P, L, T and S (`PLTS_*`), and ESA WorldCover 2021 fractions (`lc_*`) computed before the analysis switched
  to the year-before-event maps and **not used by the analysis**; `landcover.parquet` — ESA CCI / C3S
  land-cover fractions for the year before each event (`cci_*`, used by the analysis); and the metadata
  JSON files (product, schema generation, content digest of the cells).
- `runs.zip` — the confirmatory run `world_v2` (splits primary, country, temporal, llmcut), the
  sensitivity runs `world_v2_s*`, and the post hoc runs of amendment B-2 (`world_v2_recal`: recalibration-only
  reference models; `world_v2_actgrp`: areas of one flood event grouped; `world_v2_fs1`, `world_v2_fs2`: folds
  reassigned; `world_v2_dm`, `world_v2_dm_actgrp`: name features centred within each AOI): held-out cell-level predictions of every arm (`predictions/<split>/fold*.npz`), fold provenance
  (`fold*.json`), aggregated metrics with bootstrap intervals and decisions (`metrics_<split>.json`), run
  summaries and diagnostics; and the summaries of the power-analysis runs (`power_v3_*`, including the post hoc
  grid of name effects and the secondary splits).
- `records.zip` — the label-only incidence table (`incidence.json`), the language-model name readings
  (`nread/`), the GeoNames populated places within the AOIs (`geonames/aoi_subset.parquet`, derived from the
  GeoNames dump, CC BY 4.0), the frozen EMS product inventory, the outputs of verifications v31–v40, the
  sensors read from each product's source table (`v38_ems_sensors_raw.json`) and the analysis table of the
  primary split aligned with its held-out predictions (`v37_table.parquet`, used by v37 and v39).
- `MANIFEST.json` — every file in every archive with its size and sha256, and the sha256 of each archive.

## Licences and attribution of the inputs

These data are derived from: Copernicus Emergency Management Service Rapid Mapping products (contains
modified Copernicus Emergency Management Service information, 2012–2026); the Copernicus DEM GLO-30
(produced using Copernicus WorldDEM-30 © DLR e.V. 2010–2014 and © Airbus Defence and Space GmbH 2014–2018,
provided under COPERNICUS by the European Union and ESA) and the ASF GLO-30 height above nearest
drainage; ESA CCI / Copernicus Climate Change Service land cover; HydroRIVERS (HydroSHEDS); GeoNames
(CC BY 4.0); Natural Earth (public domain). Name readings were produced with the open-weight model
Gemma-4-31B-it. Users of these data should also acknowledge the original sources.
"""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest: dict = {"archives": {}}
    for name, patterns in ARCHIVES.items():
        files = sorted({p for pat in patterns for p in ROOT.glob(pat) if p.is_file()})
        rows = []
        with zipfile.ZipFile(OUT / name, "w", allowZip64=True) as z:
            for p in files:
                rel = str(p.relative_to(ROOT))
                comp = zipfile.ZIP_STORED if p.suffix in STORED else zipfile.ZIP_DEFLATED
                z.write(p, rel, compress_type=comp)
                rows.append({"path": rel, "bytes": p.stat().st_size, "sha256": sha256(p)})
        manifest["archives"][name] = {"sha256": sha256(OUT / name), "bytes": (OUT / name).stat().st_size,
                                      "n_files": len(rows), "files": rows}
        print(f"{name}: {len(rows)} files, {(OUT / name).stat().st_size / 2**30:.2f} GiB")
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    (OUT / "README.md").write_text(README, encoding="utf-8")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
