"""amendment B-10 の C5: JRC の全球河川洪水ハザード図 v2.1.2 から、AOI ごとのセルの特徴を作る。

    python world_jrc.py --prepare            # 解析する AOI にかかる 10° タイル（RP10・RP100・RP500・恒久水域）を NAS へ
    python world_jrc.py --aois 0:600         # 採用 AOI（名前の順）の範囲の特徴を作る
    python world_jrc.py --summary            # 解析するセルのうちハザード図のタイルがあるセルの割合（実行の条件）

出力: 各 AOI のラベルの場所に ``features_jrc.parquet``（cells.parquet と同じ行の順。hz_f10・hz_f100・hz_f500・
hz_d100・hz_cover）と ``features_jrc_meta.json``。タイルは ``data/global/jrc_hazard/``（配布元と同じ相対パス）、
タイルの一覧は ``tile_extents.geojson``。配布: Baugh et al. (2024)、ODbL 1.0。
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))

import world_eval as ev

from topsem.world.b10 import JRC_BASE, JRC_RPS, jrc_features, jrc_file, jrc_tiles, warp_jrc
from topsem.world.terrain import grid_lonlat_bounds

ROOT = ev.ROOT
CACHE = ROOT / "jrc_hazard"
INDEX = CACHE / "tile_extents.geojson"
OUT_NAME = "features_jrc.parquet"
WORKERS = 8
RETRIES = 6
COVER_RUN = 0.95  # 実行の条件: 解析するセルのうちハザード図のタイルがあるセルの割合（B-10 の C5）
_HEADERS = {"User-Agent": "toponym-semantics research"}


def _gn() -> pd.DataFrame:
    return pd.read_parquet(ROOT / "geonames" / "aoi_subset.parquet")  # world_eval と同じ


def included_dirs() -> list[Path]:
    dirs, _, _ = ev.select_aois(_gn())
    return dirs


def _get(url: str, dest: Path) -> int:
    for k in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers=_HEADERS)  # noqa: S310
            with urllib.request.urlopen(req, timeout=600) as resp:  # noqa: S310
                data = resp.read()
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return len(data)
        except urllib.error.HTTPError as e:
            if e.code == 404 or k == RETRIES - 1:  # noqa: PLR2004
                raise
        except (urllib.error.URLError, TimeoutError, OSError):
            if k == RETRIES - 1:
                raise
        time.sleep(2 ** k)
    return 0


def _index() -> dict:
    if not INDEX.is_file():
        INDEX.parent.mkdir(parents=True, exist_ok=True)
        _get(f"{JRC_BASE}/tile_extents.geojson", INDEX)
    return json.loads(INDEX.read_text(encoding="utf-8"))


def needed_tiles(dirs: list[Path]) -> list[tuple[int, str]]:
    idx = _index()
    tiles: set[tuple[int, str]] = set()
    for d in dirs:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        tiles |= set(jrc_tiles(idx, grid_lonlat_bounds(meta)))
    return sorted(tiles)


def prepare(workers: int) -> dict:
    tiles = needed_tiles(included_dirs())
    todo = [jrc_file(t, rp) for t in tiles for rp in (*JRC_RPS, None)]
    todo = [f for f in todo if not (CACHE / f).is_file()]
    print(f"{len(tiles)} tiles, {len(todo)} files to fetch", flush=True)
    got, nbytes, err = 0, 0, []

    def one(rel: str) -> int:
        dest = CACHE / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        return _get(f"{JRC_BASE}/{rel}", dest)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(one, f): f for f in todo}
        for fu in as_completed(futs):
            try:
                nbytes += fu.result()
                got += 1
            except Exception as e:  # noqa: BLE001
                err.append({"file": futs[fu], "error": f"{type(e).__name__}: {e}"})
    res = {"tiles": len(tiles), "fetched": got, "bytes": nbytes, "errors": err,
           "tile_names": [n for _, n in tiles]}
    (CACHE / "prepare.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k != "tile_names"}), flush=True)
    return res


def one(d: Path, idx: dict, *, force: bool = False) -> dict:
    out = d / OUT_NAME
    if not force and out.is_file():
        return {"aoi": d.name, "status": "cached"}
    t0 = time.time()
    try:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        cells = pq.read_table(d / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        g = meta["grid"]
        flat = (cells["iy"].to_numpy() * g["nx"] + cells["ix"].to_numpy()).astype(np.int64)
        tiles = jrc_tiles(idx, grid_lonlat_bounds(meta))
        files = {rp: [CACHE / jrc_file(t, rp) for t in tiles if (CACHE / jrc_file(t, rp)).is_file()]
                 for rp in (*JRC_RPS, None)}
        warped = warp_jrc(meta, files)
        feat = jrc_features(warped, meta)
        cols = {k: np.asarray(v)[flat] for k, v in feat.items()}
        pq.write_table(pa.table(cols), out, compression="zstd")
        info = {"aoi": d.name, "status": "ok", "tiles": [n for _, n in tiles],
                "files_missing": {str(rp): len(tiles) - len(v) for rp, v in files.items()},
                "n_cells": len(cells), "cover_mean": float(np.mean(cols["hz_cover"])),
                "share_f100_pos": float(np.mean(np.nan_to_num(cols["hz_f100"]) > 0)),
                "sec": round(time.time() - t0, 1),
                "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        (d / "features_jrc_meta.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        return {"aoi": d.name, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-1500:]}
    else:
        return info


def summary() -> dict:
    """実行の条件（B-10 の C5）: ハザード図のタイルがあるセルの割合。

    解析するセル（確証と同じ上限つきの表）のうち、100 年のタイルのある画素が半分以上のセルの割合。
    """
    dirs, _, infos = ev.select_aois(_gn())
    df = ev.build_table(dirs, seed=20261025, infos=infos, extra_cols=("ix", "iy"))
    ok, n = 0, 0
    per = {}
    for aoi, x in df.groupby("aoi"):
        f = ROOT / "labels" / aoi / OUT_NAME
        if not f.is_file():
            per[aoi] = None
            n += len(x)
            continue
        cells = pq.read_table(ROOT / "labels" / aoi / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        fe = pq.read_table(f, columns=["hz_cover", "hz_f100"]).to_pandas()
        fe["ix"], fe["iy"] = cells["ix"].to_numpy(), cells["iy"].to_numpy()
        m = x[["ix", "iy"]].merge(fe, on=["ix", "iy"], how="left")
        v = np.nan_to_num(m["hz_cover"].to_numpy(np.float64)) >= 0.5  # noqa: PLR2004
        per[aoi] = float(v.mean())
        ok += int(v.sum())
        n += len(x)
    res = {"rule": "amendment B-10 item C5: run only if the hazard map covers at least 95 % of the analysed cells",
           "n_cells": n, "n_covered": ok, "share_covered": ok / max(n, 1),
           "condition_met": ok / max(n, 1) >= COVER_RUN,
           "aois_missing": sorted(a for a, v in per.items() if v is None),
           "aois_below_90pct": sorted(a for a, v in per.items() if v is not None and v < 0.9)}  # noqa: PLR2004
    (CACHE / "coverage.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--aois", default="", help="採用 AOI（名前の順）の範囲 a:b")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if args.prepare:
        prepare(args.workers)
    if args.aois:
        a, b = (int(x) for x in args.aois.split(":"))
        idx = _index()
        for d in included_dirs()[a:b]:
            r = one(d, idx, force=args.force)
            print(json.dumps({k: v for k, v in r.items() if k != "trace"}), flush=True)
            if r["status"] == "error":
                print(r["trace"], flush=True)
    if args.summary:
        print(json.dumps(summary(), indent=1), flush=True)


if __name__ == "__main__":
    main()
