"""KSJ 河川（一級水系全数）→ 全国河道 bitmap ビルド（A-7・r4' F1′）。

geoshape の水系 index から全一級水系の stream.json を集め、500 m グローバル
セル（基準緯度固定の等距離格子）に量子化して重複除去した点集合を保存する。
dist_river（全国版）はこの点群への最近傍距離として計算される — 水系別
マッピングの誤り（吉田=名取・大牟田=本川遠隔）を構造的に起こさない。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/build_rivers_global.py

生成物: data/covariates/ksj_rivers_global/cells.npz (ix, iy int32)・
取得ログ cells_manifest.json（水系コードと点件数）。クラスタからは NAS 経由
（push_nas_data で covariates ごと同期）。
"""

from __future__ import annotations

import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.config import COVARIATES_DIR
from topsem.covariates.rivers import _densify, fetch_basin

OUT = COVARIATES_DIR / "ksj_rivers_global"
CELL_M = 500.0
_MIN_LINE_PTS = 2  # 単独頂点は線分でないのでスキップ
_DEG_PER_M = 1.0 / 111_320.0
_REF_LAT = 36.5  # 全国単一の固定基準（quantize 専用・距離には使わない）
_INDEX = "https://geoshape.ex.nii.ac.jp/river/resource/"


def all_basin_codes() -> list[str]:
    """resource index 直下の href="NNNNNN/" から一級水系（8 _prefix 109 個）を列挙。"""
    h = urllib.request.urlopen(_INDEX, timeout=90).read().decode()
    return sorted({c for c in re.findall(r'href="(\d{6})/"', h) if c.startswith("8")})


def quantize_points(paths: list[Path]) -> np.ndarray:
    out = []
    for p in paths:
        d = json.loads(p.read_text(encoding="utf-8"))
        for f in d["features"]:
            g = f["geometry"]
            if not g:
                continue
            lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
            for line in lines:
                c = np.asarray(line, dtype=float)
                if len(c) < _MIN_LINE_PTS:
                    continue
                pts = _densify(c, CELL_M)
                cl = CELL_M * _DEG_PER_M / np.cos(np.radians(_REF_LAT))
                gx = np.floor(pts[:, 0] / cl).astype(np.int32)
                gy = np.floor(pts[:, 1] / (_CELL_DEG_LAT)).astype(np.int32)
                out.append(np.column_stack([gx, gy]))
    if not out:
        return np.zeros((0, 2), dtype=np.int32)
    return np.unique(np.vstack(out), axis=0)


_CELL_DEG_LAT = CELL_M * _DEG_PER_M


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    codes = all_basin_codes()
    print(f"basins: {len(codes)}")
    allcells, manifest = [], {}
    for i, bc in enumerate(codes):
        try:
            paths = fetch_basin(bc)
        except (OSError, ValueError) as e:
            print(f"  {bc}: fetch failed ({e}) — skip (記録は manifest)")
            manifest[bc] = "FAILED"
            continue
        cells = quantize_points(paths)
        if len(cells):
            allcells.append(cells)
        manifest[bc] = len(cells)
        print(f"  [{i + 1}/{len(codes)}] {bc}: {len(cells)} cells", flush=True)
    if not allcells:
        msg = "全水系の取得に失敗"
        raise RuntimeError(msg)
    pts = np.unique(np.vstack(allcells), axis=0)
    np.savez_compressed(OUT / "cells.npz", ix=pts[:, 0].astype(np.int32), iy=pts[:, 1].astype(np.int32))
    (OUT / "cells_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"global river bitmap: {len(pts)} cells (500 m) -> {OUT / 'cells.npz'}")


if __name__ == "__main__":
    main()
