"""KSJ W05 河川線形（Geoshape 経由）→ 河川からの距離（M1 候補共変量）。

取得構造: `river/resource/{6桁水系}/` のページに 10 桁河川コード一覧、
`river/resource/{6桁}/{10桁}/stream.json` にその河川の LineString 群。
水系ごとに全河川ファイルをキャッシュし、線分を SAMPLE_M で sampled した
点群を cKDTree にして最近傍距離を出す（区間近似誤差 ≤ SAMPLE_M/2）。

時点注意: KSJ W05 は 2006–09 年版 → 2010 年以降のイベントに対し
河道形状は事前情報（point_in_time.md 参照）。改修で動いた河道は誤差側。
**凍結済み M1（elev, depress）には未接続** — M1 共変量リスト確定（protocol §7
TBD）で amd 後に対応づけるときに build_xy へ加える。
"""

from __future__ import annotations

import json
import math
import re
import urllib.request
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from topsem.config import COVARIATES_DIR

RIVERS_DIR = COVARIATES_DIR / "ksj_rivers"
SAMPLE_M = 25.0
_DEG_PER_M = 1.0 / 111_320.0
_BASE = "https://geoshape.ex.nii.ac.jp/river/resource/"

# ラベルメタの river 表記 → Geoshape 水系 6 桁コード（2026-09-20 検索ページで確認）
BASIN_CODE = {
    "久慈川水系（久慈川）": "830301",
    "信濃川水系（千曲川）": "840403",
    "吉田川水系（吉田川）": "820203",  # 鳴瀬川水系（宮城・山本町の吉田川）。
    # 2026-09-30 r4' F1: 820202（名取）は誤り（wet セントロイドから河道 23.8 km）。
    # geoshape 8202030004 の W05_004=吉田川・bbox 一致で実測確認（v15 監査）。
    "荒川水系（入間川・越辺川・都幾川）": "830304",
    "那珂川水系（那珂川）": "830302",
    "阿武隈川水系（阿武隈川）": "820201",
    "rinkaku_mogamigawa01": "820211",  # 最上川
    # 2020/2024 イベント水系（A-4 発効で event-loo の訓練プールが
    # 2020+2024 だけの fold が全 NaN 化して sklearn binning がクラッシュした
    # 件 v14 で特定。コードは geoshape resource index を 2026-09-26 実照合）
    "最上川水系（最上川）": "820211",
    "筑後川水系（筑後川）": "890906",
    "球磨川水系（球磨川）": "890911",
}


def basin_river_codes(basin: str) -> list[str]:
    """水系ページから 10 桁河川コード一覧を取る（キャッシュ {basin}_codes.json）。"""
    cache = RIVERS_DIR / f"{basin}_codes.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    url = _BASE + f"{basin}/"
    html = urllib.request.urlopen(url, timeout=60).read().decode()  # noqa: S310
    codes = sorted(set(re.findall(r'href="(\d{6}\d{4})/"', html)))
    RIVERS_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(codes), encoding="utf-8")
    return codes


def fetch_basin(basin: str, *, force: bool = False) -> list[Path]:
    """水系の全 stream.json を data/covariates/ksj_rivers/ に集めてパス返す。"""
    RIVERS_DIR.mkdir(parents=True, exist_ok=True)
    out = []
    for rc in basin_river_codes(basin):
        p = RIVERS_DIR / f"{rc}.json"
        if force or not p.exists():
            url = _BASE + f"{basin}/{rc}/stream.json"
            p.write_bytes(urllib.request.urlopen(url, timeout=120).read())  # noqa: S310
        out.append(p)
    return out


_MIN_LINE_PTS = 2


def _densify(coords: np.ndarray, cell_m: float) -> np.ndarray:
    """折れ線頂点 (lon,lat) を cell_m 刻みで補間した点群 (lon,lat) に。"""
    dx = np.diff(coords[:, 0]) * math.cos(math.radians(coords[:, 1].mean())) / _DEG_PER_M
    dy = np.diff(coords[:, 1]) / _DEG_PER_M
    seg = np.hypot(dx, dy)
    pts = [coords[0]]
    for i, ln in enumerate(seg):
        if ln <= 0:
            continue
        k = int(np.ceil(ln / cell_m))
        t = np.linspace(0, 1, k + 1, endpoint=False)[1:]
        pts.append(coords[:-1][i] + np.outer(t, np.diff(coords, axis=0)[i]))
    return np.vstack(pts) if len(pts) > 1 else coords


def basin_points(paths: list[Path]) -> np.ndarray:
    """stream.json 群 → sampled 点群（lon,lat）。

    窓外の線分も点群に含めてよい — 距離は幾何の情報源であって
    訓練側 fit 処理はなく、窓外情報で「学習」する対象がない（point-in-time
    台帳上 KSJ W05 は 2006–09 時点で公開）。
    """
    parts = []
    for p in paths:
        d = json.loads(p.read_text(encoding="utf-8"))
        for f in d["features"]:
            g = f["geometry"]
            if not g:
                continue
            lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
            for line in lines:
                c = np.asarray(line, dtype=float)
                if len(c) >= _MIN_LINE_PTS:
                    parts.append(_densify(c, SAMPLE_M))
    return np.vstack(parts) if parts else np.zeros((0, 2))


def river_distance_m(
    lon: np.ndarray,
    lat: np.ndarray,
    pts: np.ndarray,
    lat_ref: float,
) -> np.ndarray:
    """各格子点から最も近い河道線分までの距離（m、平面近似・基準緯度補正）。"""
    if len(pts) == 0:
        return np.full(len(lon), np.nan)
    sx = math.cos(math.radians(lat_ref)) / _DEG_PER_M
    sy = 1.0 / _DEG_PER_M
    tree = cKDTree(np.column_stack([pts[:, 0] * sx, pts[:, 1] * sy]))
    q = np.column_stack([lon * sx, lat * sy])
    d, _ = tree.query(q, k=1)
    return d


_GNPZ = COVARIATES_DIR / "ksj_rivers_global" / "cells.npz"
_GLOBAL_CELL_M = 500.0
_GLOBAL_REF_LAT = 36.5
_global_cache: dict = {}


def global_river_points() -> np.ndarray:
    """全国河道 bitmap（build_rivers_global.py 生成）→ (lon,lat) 点群。無ければ空。"""
    if "pts" not in _global_cache:
        if _GNPZ.exists():
            z = np.load(_GNPZ)
            cl = _GLOBAL_CELL_M * _DEG_PER_M / math.cos(math.radians(_GLOBAL_REF_LAT))
            lon = (z["ix"].astype(np.float64) + 0.5) * cl
            lat = (z["iy"].astype(np.float64) + 0.5) * _GLOBAL_CELL_M * _DEG_PER_M
            _global_cache["pts"] = np.column_stack([lon, lat])
        else:
            _global_cache["pts"] = np.zeros((0, 2))
    return _global_cache["pts"]
