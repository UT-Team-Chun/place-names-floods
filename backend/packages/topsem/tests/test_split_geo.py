"""topsem.splits.geo — 幾何排他マスク（A-7 / r4' F3 回帰テスト）。"""

from __future__ import annotations

import math

import pandas as pd

from topsem.splits.geo import cell_meters, exclude_near_holdout

_DEG = 1.0 / 111_320.0


def _unit(ix0: int, iy0: int, n: int, lat_ref: float) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ix": [ix0 + i for i in range(n)],
            "iy": [iy0 + i for i in range(n)],
            "lat_ref": [lat_ref] * n,
        }
    )


def test_cell_meters_is_lat_ref_invariant() -> None:
    """同一物理点は unit の lat_ref が違っても同じメートル座標に戻る。

    v4 まで漏れていた穴：ix,iy は lat_ref 依存のローカル指標で、
    隣接窓の key 一致チェックでは重複を検出できなかった。
    """
    lon, lat = 140.5, 37.2
    for ref in (36.0, 37.0, 38.5):
        cl = 250.0 * _DEG / math.cos(math.radians(ref))
        ix, iy = int(lon / cl), int(lat / (_DEG * 250.0))
        xm, ym = cell_meters([ix], [iy], [ref])
        # 量化誤差 ≤1 cell（250 m）以内で復元できる
        assert abs(xm[0] - lon * 111_320.0 * math.cos(math.radians(36.5))) < 300.0
        assert abs(ym[0] - lat * 111_320.0) < 200.0  # 半セル（125 m）量子化


def test_exclude_near_holdout_removes_shared_and_adjacent() -> None:
    hold = _unit(1000, 1000, 3, 36.5)
    # 同一位置セル（同格納値）と遠方セル
    tr = pd.concat([_unit(1000, 1000, 3, 36.5), _unit(1400, 1000, 3, 36.5)], ignore_index=True)
    kept, n = exclude_near_holdout(tr, hold, buffer_m=250.0)
    assert n == 3  # 同一位置の 3 cell が消える
    assert len(kept) == 3
    assert set(kept["ix"]) == {1400, 1401, 1402}


def test_cross_lat_ref_same_place_is_excluded() -> None:
    """lat_ref の違う unit が同一物理地面を共有する場合を捕まえる（F3 の核心）。"""
    # ref=36.5 の ix=1000 番地の lon を、ref=37.5 で同じ ix として書くと別場所。
    # 正しくは lon から再量子化 → ix がずれる。その「ずれた同一点」を構成する。
    ref_a, ref_b = 36.5, 37.5
    lon, lat = 140.0, 37.0
    qa = (int(lon // (250.0 * _DEG / math.cos(math.radians(ref_a)))),
          int(lat // (250.0 * _DEG)))
    qb = (int(lon // (250.0 * _DEG / math.cos(math.radians(ref_b)))),
          int(lat // (250.0 * _DEG)))
    assert qa != qb  # キー一致チェックでは絶対に合わない（v4 の穴そのもの）
    hold = pd.DataFrame({"ix": [qa[0]], "iy": [qa[1]], "lat_ref": [ref_a]})
    tr = pd.DataFrame({"ix": [qb[0], qb[0] + 900], "iy": [qb[1], qb[1]], "lat_ref": [ref_b] * 2})
    kept, n = exclude_near_holdout(tr, hold, buffer_m=250.0)
    assert n == 1  # 近接（同一地点）だけ除去され遠方は残る
    assert len(kept) == 1


def test_zero_buffer_is_noop() -> None:
    hold = _unit(1000, 1000, 2, 36.5)
    tr = _unit(1000, 1000, 2, 36.5)
    kept, n = exclude_near_holdout(tr, hold, buffer_m=0.0)
    assert n == 0
    assert len(kept) == 2
