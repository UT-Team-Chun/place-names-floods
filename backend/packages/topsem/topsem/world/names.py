"""GeoNames 全件ダンプ → AOI 格子セルへの地名割当（最近傍 K 名 + 集落幾何）。

- 主: 地物クラス P（集落）。感度: P+L+T+S（地域名・地形名・地点名）
- **H（水域名）は除外**: 川・湖・沼そのものの名前は水体の位置を直接ラベルしており、
  「人が地名に残した記憶」ではなく水文データの再表現になる（漏洩）
- 文字列は GeoNames の ``asciiname``（ASCII 化された綴り）を正とし、``name``
  （現地表記）を併記する

セルごとに保存するもの: 最近傍 K=3 の geonameid と距離（m）、1 km / 3 km 内の件数。
意味特徴（辞書照合・n-gram 得点）はここでは作らない（`topsem.world.lexicon` 側）。
"""

from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree

COLS = {0: "gid", 1: "name", 2: "asciiname", 4: "lat", 5: "lon", 6: "fclass", 7: "fcode",
        8: "cc", 14: "population"}
CLASS_SETS = {"P": ("P",), "PLTS": ("P", "L", "T", "S")}
K = 3


def geonames_to_parquet(zip_path: Path, out: Path) -> int:
    """allCountries.zip → 必要列だけの parquet（水域 H も含めて保存し、使用時に絞る）。"""
    with zipfile.ZipFile(zip_path) as zf, zf.open("allCountries.txt") as raw:
        txt = io.TextIOWrapper(raw, encoding="utf-8")
        df = pd.read_csv(
            txt, sep="\t", header=None, usecols=list(COLS), names=None,
            quoting=csv.QUOTE_NONE, dtype={1: "string", 2: "string", 6: "string",
                                           7: "string", 8: "string"},
            keep_default_na=False, na_values=[], low_memory=False,
        )
    df = df.rename(columns=COLS)
    df["population"] = pd.to_numeric(df["population"], errors="coerce").fillna(0).astype("int64")
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, compression="zstd", index=False)
    return len(df)


def subset_for_bounds(gn: pd.DataFrame, bounds: list[tuple[float, float, float, float]],
                      pad_deg: float = 0.3) -> pd.DataFrame:
    """AOI 群の外接範囲（+余白）に入る地名だけを残す（ワーカーのメモリ節約）。"""
    keep = np.zeros(len(gn), dtype=bool)
    lon = gn["lon"].to_numpy(dtype=float)
    lat = gn["lat"].to_numpy(dtype=float)
    for w, s, e, n in bounds:
        keep |= (lon >= w - pad_deg) & (lon <= e + pad_deg) & (lat >= s - pad_deg) & (
            lat <= n + pad_deg)
    return gn.loc[keep].reset_index(drop=True)


def assign_names(cells_xy: tuple[np.ndarray, np.ndarray], crs: CRS, gn: pd.DataFrame,
                 bounds: tuple[float, float, float, float], pad_deg: float = 0.25,
                 ) -> dict[str, np.ndarray]:
    """セル中心（LAEA m）ごとに最近傍 K 名と近傍件数（クラス集合ごと）。"""
    w, s, e, n = bounds
    sub = gn.loc[
        gn["lon"].between(w - pad_deg, e + pad_deg) & gn["lat"].between(s - pad_deg, n + pad_deg)
    ]
    tr = Transformer.from_crs(4326, crs, always_xy=True)
    x, y = cells_xy
    q = np.column_stack([x, y])
    out: dict[str, np.ndarray] = {}
    for tag, classes in CLASS_SETS.items():
        pts = sub.loc[sub["fclass"].isin(classes)]
        if len(pts) == 0:
            for j in range(1, K + 1):
                out[f"{tag}_gid{j}"] = np.full(len(x), -1, dtype=np.int64)
                out[f"{tag}_d{j}"] = np.full(len(x), np.inf, dtype=np.float32)
            out[f"{tag}_n1km"] = np.zeros(len(x), dtype=np.int32)
            out[f"{tag}_n3km"] = np.zeros(len(x), dtype=np.int32)
            continue
        px, py = tr.transform(pts["lon"].to_numpy(dtype=float), pts["lat"].to_numpy(dtype=float))
        tree = cKDTree(np.column_stack([px, py]))
        k = min(K, len(pts))
        d, idx = tree.query(q, k=k)
        d = d.reshape(len(x), k)
        idx = idx.reshape(len(x), k)
        gids = pts["gid"].to_numpy(dtype=np.int64)
        for j in range(K):
            if j < k:
                out[f"{tag}_gid{j + 1}"] = gids[idx[:, j]]
                out[f"{tag}_d{j + 1}"] = d[:, j].astype(np.float32)
            else:
                out[f"{tag}_gid{j + 1}"] = np.full(len(x), -1, dtype=np.int64)
                out[f"{tag}_d{j + 1}"] = np.full(len(x), np.inf, dtype=np.float32)
        out[f"{tag}_n1km"] = np.array([len(v) for v in tree.query_ball_point(q, 1000.0)],
                                      dtype=np.int32)
        out[f"{tag}_n3km"] = np.array([len(v) for v in tree.query_ball_point(q, 3000.0)],
                                      dtype=np.int32)
    return out
