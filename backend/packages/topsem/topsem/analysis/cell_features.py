"""(event × watershed) 窓格子 → 学習 table の共通部品（pilot/本体解析で共有）。

特徴は訓練側だけで計算される（この関数群に fit 処理はない。埋め込みは
ハッシュ固定次元で fit 不要、route_a も辞書参照のみ）。
"""

from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd

from topsem.config import LABELS_DIR
from topsem.covariates.dem import relief_depression, sample_grid
from topsem.covariates.rivers import (
    BASIN_CODE,
    basin_points,
    fetch_basin,
    global_river_points,
    river_distance_m,
)
from topsem.features.route_a import ElementIndex, name_features
from topsem.features.route_b import DIM, embed_name

WET_TH = 0.5
_DEG_PER_M = 1.0 / 111_320.0
_VDATE_LEN = 8  # source.version_date の YYYYMMDD 長


def discover_label_units() -> list[dict]:
    """data/labels/*/derived/*/metadata.json を巡回し unit 記述を返す。

    unit ごとの event_date（版日付）を返す。イベント時刻が必要な呼び手側では
    **最小版日付（= 発生開始）**を使うこと — 複数の水系を跨ぐイベント
    （2020 年 7 月豪雨など）で temporal の訓練/評価境界が dict 反復順で
    揺れないため（`scripts/freeze_splits.py` 参照）。
    """
    units = []
    for event_dir in sorted(p for p in LABELS_DIR.iterdir() if p.is_dir()):
        derived = event_dir / "derived"
        if not derived.is_dir():
            continue
        for d in sorted(derived.iterdir()):
            meta = json.loads((d / "metadata.json").read_text(encoding="utf-8"))
            vdate = str(meta["source"].get("version_date", ""))
            date = f"{vdate[0:4]}-{vdate[4:6]}-{vdate[6:8]}" if len(vdate) == _VDATE_LEN else ""
            units.append(
                {
                    "event": meta["event"],
                    "event_date": date,
                    "watershed": meta["river"],
                    # 運用者が build 時に渡した top-level role が正（sources 側の既定値は
                    # 常時 robustness_check_only で、2024_mogami の snapshot_only が
                    # 覆されていた — レビュー C-2・2026-09-23）
                    "role": meta.get("role") or meta["source"].get("role", ""),
                    "meta": meta,
                    "grid": pd.read_csv(d / "grid_window.csv"),
                }
            )
    return units


def river_cells(
    unit: dict,
    pop: pd.DataFrame,
    index: ElementIndex,
    feat_cache: dict,
    emb_cache: dict,
) -> pd.DataFrame:
    """1 (event, watershed) の窓グリッドに住所名特徴とラベルを接合する。"""
    m, grid = unit["meta"], unit["grid"]
    size, lat_ref = m["grid_size_m"], m["lat_ref"]
    minx, miny, maxx, maxy = m["bbox"]
    pad = 0.03
    sub = pop.loc[
        pop["lat"].between(miny - pad, maxy + pad) & pop["lon"].between(minx - pad, maxx + pad)
    ]
    cell_lon = size * _DEG_PER_M / math.cos(math.radians(lat_ref))
    ix = np.floor(sub["lon"].to_numpy(dtype=float) / cell_lon).astype(int)
    iy = np.floor(sub["lat"].to_numpy(dtype=float) / (size * _DEG_PER_M)).astype(int)

    rows = []
    for name, kana, x, y in zip(sub["name"], sub["kana"], ix, iy, strict=True):
        n = "" if pd.isna(name) else str(name)
        k = None if pd.isna(kana) else str(kana)
        key = (n, k)
        if key not in feat_cache:
            f = name_features(n, k, index)
            feat_cache[key] = (f["a_flood_spec_sum"], f["a_any_hit"], f["a_n_substituted"])
        rows.append((x, y, *feat_cache[key]))
    g = _aggregate_to_grid(grid, rows, sub, emb_cache)
    g["y"] = (g["wet_frac"] >= WET_TH).astype(int)
    g["event"] = unit["event"]
    g["watershed"] = unit["watershed"]
    g["role"] = unit["role"]
    # ix,iy は unit の lat_ref 依存のローカル指標。幾何比較（split v5・A-7）は
    # この列から共通メートル座標へ復元して行う（topsem.splits.geo）
    g["lat_ref"] = float(lat_ref)
    g = add_dem(g, m)
    if os.environ.get("TOPSEM_M1_RIVER", "").lower() in ("1", "global"):
        g = add_river_distance(g, m)
    return g


def population_tag() -> str:
    """load_population が選択したファイルのタグ（来歴記録用・r3-M2）。"""
    from topsem.config import DATA_DIR  # noqa: PLC0415

    cands = sorted((DATA_DIR / "population").glob("*/population.csv.gz"))
    return cands[-1].parent.name if cands else "unknown"


def _aggregate_to_grid(
    grid: pd.DataFrame,
    rows: list,
    sub: pd.DataFrame,
    emb_cache: dict,
) -> pd.DataFrame:
    if not rows:
        g = grid.copy()
        for cname in ["a_flood", "a_any", "a_sub", "n_addr"]:
            g[cname] = 0.0
        g["emb"] = [np.zeros(DIM, dtype=np.float32)] * len(g)
        return g
    rdf = pd.DataFrame(rows, columns=["ix", "iy", "a_flood", "a_any", "a_sub"])
    agg = (
        rdf.groupby(["ix", "iy"])
        .agg(
            a_flood=("a_flood", "max"),
            a_any=("a_any", "sum"),
            a_sub=("a_sub", "sum"),
            n_addr=("a_flood", "size"),
        )
        .reset_index()
    )
    names_in_order = ["" if pd.isna(v) else str(v) for v in sub["name"]]
    emb_rows = np.vstack(
        [
            emb_cache.setdefault(nm, np.asarray(embed_name(nm), dtype=np.float32))
            for nm in names_in_order
        ]
    )
    codes, keys = pd.factorize(
        pd.Series(list(zip(rdf["ix"], rdf["iy"], strict=True)), dtype=object)
    )
    codes = codes.astype(int)
    cell_emb = np.zeros((len(keys), DIM), dtype=np.float32)
    np.add.at(cell_emb, codes, emb_rows)
    counts = np.bincount(codes, minlength=len(keys)).astype(np.float32)
    cell_emb /= counts[:, None]
    kmap = {kk: i for i, kk in enumerate(keys)}
    g = grid.merge(agg, on=["ix", "iy"], how="left")
    for cname in ["a_flood", "a_any", "a_sub", "n_addr"]:
        g[cname] = g[cname].fillna(0.0)
    embs_col = np.zeros((len(g), DIM), dtype=np.float32)
    ix_arr = g["ix"].to_numpy()
    iy_arr = g["iy"].to_numpy()
    hit = ~np.isnan(ix_arr.astype(np.float64))
    for j in np.flatnonzero(hit):
        ki = kmap.get((int(ix_arr[j]), int(iy_arr[j])))
        if ki is not None:
            embs_col[j] = cell_emb[ki]
    g["emb"] = list(embs_col)
    return g


def add_dem(g: pd.DataFrame, m: dict) -> pd.DataFrame:
    """窓格子の標高・相対低地度を接合する。"""
    size = m["grid_size_m"]
    cell_lat = size * _DEG_PER_M
    cell_lon = size * _DEG_PER_M / math.cos(math.radians(m["lat_ref"]))
    ixs = g["ix"].to_numpy(dtype=int)
    iys = g["iy"].to_numpy(dtype=int)
    ix0, ix1 = ixs.min(), ixs.max()
    iy0, iy1 = iys.min(), iys.max()
    gx, gy = np.meshgrid(np.arange(ix0, ix1 + 1), np.arange(iy0, iy1 + 1))
    elev = sample_grid(gx.ravel(), gy.ravel(), cell_lon, cell_lat)
    elev2d = elev.reshape(gy.shape)
    g["elev"] = elev2d[iys - iy0, ixs - ix0]
    g["depress"] = relief_depression(elev2d)[iys - iy0, ixs - ix0]
    return g


def build_xy(  # noqa: PLR0911
    df: pd.DataFrame, arm: str, *, drop_names: bool = False
) -> np.ndarray:
    """腕名 → 設計行列。drop_names は消融去（M3 → 名前特征なし）用。

    A-4 発効（dist_river 列あり = env TOPSEM_M1_RIVER=1）時、M1/M2/M3 は
    河道距離 log1p を地形と併せ、消融去腕 M1r は base+距離のみ、M1rN は
    base+距離+地名（地形なし；r3-F4 の分解腕）を設計行列とする。
    """
    base = np.column_stack([df["ix"], df["iy"]]).astype(np.float32)
    if arm == "M0":
        return base
    terr = np.column_stack([df["elev"], df["depress"]]).astype(np.float32)
    riv = (
        np.log1p(df["dist_river"].to_numpy(dtype=np.float32))[:, None]
        if "dist_river" in df.columns
        else np.zeros((len(df), 0), dtype=np.float32)
    )
    if arm == "M0":
        return base
    if drop_names:  # M3 の消融去 = 名前特征を落とすと M1 と同一設計
        return np.column_stack([base, terr, riv])
    emb = np.vstack(df["emb"].to_numpy())
    match arm:
        case "M1":
            return np.column_stack([base, terr, riv])
        case "M1r":
            return np.column_stack([base, riv])
        case "M1rN":
            return np.column_stack([base, riv, emb])
        case "M2":
            return np.column_stack(
                [base, terr, riv,
                 df[["a_flood", "a_any", "a_sub"]].to_numpy(dtype=np.float32)]
            )
        case _:  # M3
            return np.column_stack([base, terr, riv, emb])


def load_population() -> pd.DataFrame:
    """最新 population.csv.gz を正規化列つきで返す。"""
    import gzip  # noqa: PLC0415

    from topsem.config import DATA_DIR  # noqa: PLC0415

    cands = sorted((DATA_DIR / "population").glob("*/population.csv.gz"))
    if not cands:
        msg = "data/population/<tag>/population.csv.gz がない。make population を先に実行"
        raise FileNotFoundError(msg)
    with gzip.open(cands[-1], "rt", encoding="utf-8") as f:
        pop = pd.read_csv(f)
    pop["name"] = pop["大字町丁目名"].astype("string")
    pop["kana"] = pop["大字町丁目名カナ"].astype("string")
    pop["lat"] = pd.to_numeric(pop["緯度"], errors="coerce")
    pop["lon"] = pd.to_numeric(pop["経度"], errors="coerce")
    return pop


_RP: dict[str, np.ndarray] = {}


def add_river_distance(g: pd.DataFrame, m: dict) -> pd.DataFrame:
    """窓全セルに KSJ 河道からの距離（m）を付す。A-4（M1+dist_river）用。

    水系名が BASIN_CODE に無い unit（2020 の sub-basin 等）は NaN →
    build_xy 側で elev と同じく教訓的に扱うのではなく、**调用側で
    母集団から除外する**（silent partial coverage を防ぐ）。
    """
    if os.environ.get("TOPSEM_M1_RIVER", "").lower() == "global":
        # A-7: 全国河道 bitmap 最近傍（水系マッピング誤りに頑健）。v7 既定。
        size = m["grid_size_m"]
        cell_lat = size * _DEG_PER_M
        cell_lon = size * _DEG_PER_M / math.cos(math.radians(m["lat_ref"]))
        glon = (g["ix"].to_numpy(dtype=float) + 0.5) * cell_lon
        glat = (g["iy"].to_numpy(dtype=float) + 0.5) * cell_lat
        pts = global_river_points()
        if len(pts) == 0:
            msg = (
                "TOPSEM_M1_RIVER=global なのに ksj_rivers_global/cells.npz が無い"
                "（scripts/build_rivers_global.py を実行）"
            )
            raise RuntimeError(msg)
        g["dist_river"] = river_distance_m(glon, glat, pts, m["lat_ref"])
        return g
    code = BASIN_CODE.get(m["river"])
    size = m["grid_size_m"]
    cell_lat = size * _DEG_PER_M
    cell_lon = size * _DEG_PER_M / math.cos(math.radians(m["lat_ref"]))
    lon = (g["ix"].to_numpy(dtype=float) + 0.5) * cell_lon
    lat = (g["iy"].to_numpy(dtype=float) + 0.5) * cell_lat
    if code is None:
        #  unmapped水系の silent NaN は退化列として sklearn binning を
        # クラッシュさせる（v14・2026-09-26）。拡張時はここへ登録すること
        msg = (
            f"BASIN_CODE 未登録の river キー: {m['river']!r}"
            "（geoshape resource index で照合して追加）"
        )
        raise KeyError(msg)
    if code not in _RP:
        _RP[code] = basin_points(fetch_basin(code))
    g["dist_river"] = river_distance_m(lon, lat, _RP[code], m["lat_ref"])
    return g
