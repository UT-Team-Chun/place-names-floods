"""KSJ 河川距離共変量の幾何部分（ネットワーク不要・一時ファイルのみ）。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from topsem.covariates.rivers import (
    _DEG_PER_M,
    BASIN_CODE,
    _densify,
    basin_points,
    river_distance_m,
)


def test_basin_code_mapping_sanity() -> None:
    """r4'-F1 回帰テスト: 吉田川水系は鳴瀬川 820203（名取 820202 は同姓異川の罠）。

    地理的監査（水系 bbox 中心 vs河道最近傍 <5 km）は
    backend/verification/r_okauchi/v15_river_mapping_audit.py で実施する。
    """
    assert BASIN_CODE["吉田川水系（吉田川）"] == "820203"
    assert "820202" not in BASIN_CODE.values()  # 名取川はラベル水系に無い
    # 共有を許すのは最上の alias キーのみ（rinkaku_mogamigawa01 ≡ 最上川水系）
    codes = [c for k, c in BASIN_CODE.items() if k != "rinkaku_mogamigawa01"]
    assert len(set(codes)) == len(codes)


def _fake_stream(tmp: Path, lines: list) -> Path:
    feats = [
        {"geometry": {"type": "LineString", "coordinates": ln}, "properties": {}}
        for ln in lines
    ]
    p = tmp / "stream.json"
    p.write_text(json.dumps({"type": "FeatureCollection", "features": feats}), encoding="utf-8")
    return p


def test_densify_spacing() -> None:
    # 東西に約 1 km の直線 → 25 m 刻みなら 40 点前後
    lon0, lat0 = 140.0, 36.0
    length_deg = 1000.0 * _DEG_PER_M / np.cos(np.radians(lat0))
    coords = np.array([[lon0, lat0], [lon0 + length_deg, lat0]])
    pts = _densify(coords, 25.0)
    assert 35 <= len(pts) <= 45
    assert np.allclose(pts[0], coords[0])


def test_basin_points_and_distance(tmp_path: Path) -> None:
    lat0, lon0 = 36.0, 140.0
    line = [[lon0 - 0.01, lat0], [lon0 + 0.01, lat0]]  # 東西の河道
    p = _fake_stream(tmp_path, [line])
    pts = basin_points([p])
    assert len(pts) > 10
    # 真北 500 m の点 → 距離 ≈ 500 m（平面近似・同緯度基準で誤差小）
    lat_n = lat0 + 500.0 * _DEG_PER_M
    d = river_distance_m(np.array([lon0]), np.array([lat_n]), pts, lat_ref=lat0)
    assert abs(d[0] - 500.0) < 5.0


def test_distance_empty_is_nan() -> None:
    d = river_distance_m(np.array([140.0]), np.array([36.0]), np.zeros((0, 2)), lat_ref=36.0)
    assert np.isnan(d).all()
