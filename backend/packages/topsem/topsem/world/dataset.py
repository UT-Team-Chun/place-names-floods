"""v2 全球検定の解析テーブル: AOI 採用・評価セル抽出・場所グループ・fold 割当。

すべて preregistration_v2.md §2.4・§5 の規則の実装。ここで使うのはラベル単独の
集計（評価セル数・浸水セル数・種別の面積比）と位置だけで、特徴との関係は見ない。
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

EVAL_MIN_AOI = 0.5
EVAL_MIN_OBSERVED = 0.5
EVAL_MAX_PERM_WATER = 0.5
WET_TH = 0.5
MIN_EVAL_CELLS = 400
MIN_WET_CELLS = 20
MAX_WET_FRAC = 0.8
MIN_RIVERINE_SHARE = 0.5
INLAND_KM = 10.0  # 種別未記載の製品を河川とみなす条件: AOI 重心が海岸線から 10 km 超
GROUP_BUFFER_KM = 5.0
N_FOLDS = 10
# §2.4-4: 場所グループあたりの評価セルの上限（グループ内の AOI へ評価セル数に比例配分）。
# 草案の「AOI あたり 2% 以下の最大値」（265,000）は実効グループ数を 34/194 に落とすので、
# ラベル単独の発生表で比べて凍結前に差し替えた（v32: Kish 実効数 164/194、148 万セル）
GROUP_CAP = 10_000


@dataclass(frozen=True)
class AoiInfo:
    tag: str
    activation: str
    countries: tuple[str, ...]
    event_time: str
    bounds: tuple[float, float, float, float]
    n_eval: int
    n_wet: int
    riverine_share: float
    footprint_implicit: bool
    unspecified_share: float = 0.0
    coast_km: float = float("inf")
    footprint_missing: bool = False  # 第 0・第 1 世代（撮像 footprint 層なし）
    all_types: bool = False  # 感度 §7-1: 全浸水種別のラベル

    @property
    def unspecified_as_riverine(self) -> bool:
        """種別未記載（第 1 世代など）が過半で、かつ内陸 → 河川浸水とみなす（§2.4）。"""
        return (self.riverine_share < MIN_RIVERINE_SHARE
                and self.unspecified_share >= MIN_RIVERINE_SHARE
                and self.coast_km > INLAND_KM)


def evaluable(df: pd.DataFrame) -> np.ndarray:
    return (
        (df["f_aoi"].to_numpy() >= EVAL_MIN_AOI)
        & (df["f_observed"].to_numpy() >= EVAL_MIN_OBSERVED)
        & (df["f_perm_water"].to_numpy() < EVAL_MAX_PERM_WATER)
    )


def riverine_labels(df: pd.DataFrame, *, unspecified_as_riverine: bool = False,
                    all_types: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """主ラベル: y=1 河川浸水被覆 ≥ 0.5、y=0 全種別の浸水被覆 < 0.5、それ以外は除外。

    unspecified_as_riverine（内陸の種別未記載製品）では未記載の浸水も河川に数える。
    all_types（感度 §7-1）では種別を問わず y=1 ⇔ 浸水被覆 ≥ 0.5（除外なし）。
    """
    if all_types:
        y1 = df["f_flood_any"].to_numpy() >= WET_TH
        return y1.astype(np.int8), np.ones(len(df), dtype=bool)
    riv = df["f_riverine"].to_numpy()
    if unspecified_as_riverine:
        riv = riv + df["f_other"].to_numpy()
    y1 = riv >= WET_TH
    y0 = df["f_flood_any"].to_numpy() < WET_TH
    return y1.astype(np.int8), (y1 | y0)


def aoi_info(d: Path, coast: "CoastIndex | None" = None, *,  # noqa: UP037
             all_types: bool = False) -> AoiInfo:
    """ラベル単独の AOI 要約（採用判定と発生表に使う）。all_types は感度 §7-1。"""
    meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
    cells = pq.read_table(d / "cells.parquet").to_pandas()
    ev = evaluable(cells)
    flood_area = float(cells.loc[ev, "f_flood_any"].sum())
    riv_area = float(cells.loc[ev, "f_riverine"].sum())
    oth_area = float(cells.loc[ev, "f_other"].sum())
    lon, lat = cells["lon"].to_numpy(), cells["lat"].to_numpy()
    coast_km = coast.distance_km(float(np.mean(lon)), float(np.mean(lat))) if coast else np.inf
    riv_share = riv_area / flood_area if flood_area > 0 else 0.0
    oth_share = oth_area / flood_area if flood_area > 0 else 0.0
    as_riv = (riv_share < MIN_RIVERINE_SHARE and oth_share >= MIN_RIVERINE_SHARE
              and coast_km > INLAND_KM)
    y, keep = riverine_labels(cells, unspecified_as_riverine=as_riv, all_types=all_types)
    sel = ev & keep
    act = meta.get("activation", {})
    return AoiInfo(
        tag=d.name,
        activation=str(act.get("code", d.name.split("_")[0])),
        countries=tuple(act.get("countries") or ()),
        event_time=str(act.get("event_time") or act.get("activation_time") or ""),
        bounds=(float(lon.min()), float(lat.min()), float(lon.max()), float(lat.max())),
        n_eval=int(sel.sum()),
        n_wet=int(y[sel].sum()),
        riverine_share=riv_share,
        footprint_implicit=bool(meta.get("footprint_implicit", False)),
        unspecified_share=oth_share,
        coast_km=float(coast_km),
        footprint_missing=bool(meta.get("footprint_missing", False)),
        all_types=all_types,
    )


def include(info: AoiInfo, *, n_names_inside: int) -> tuple[bool, str]:
    """§2.4 の採用基準。理由文字列を返す（除外理由の集計に使う）。"""
    if info.n_eval < MIN_EVAL_CELLS:
        return False, "eval<400"
    if info.n_wet < MIN_WET_CELLS:
        return False, "wet<20"
    if info.n_wet / info.n_eval > MAX_WET_FRAC:
        return False, "wetfrac>0.8"
    if (not info.all_types and info.riverine_share < MIN_RIVERINE_SHARE
            and not info.unspecified_as_riverine):
        return False, "not-riverine"
    if n_names_inside < 1:
        return False, "no_P_name"
    return True, "ok"


def location_groups(bounds: dict[str, tuple[float, float, float, float]],
                    buffer_km: float = GROUP_BUFFER_KM) -> dict[str, int]:
    """緩衝した外接範囲が交差する AOI を union-find で束ねる（activation を跨ぐ）。"""
    tags = sorted(bounds)
    parent = list(range(len(tags)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    boxes = []
    for t in tags:
        w, s, e, n = bounds[t]
        dlat = buffer_km / 111.32
        dlon = buffer_km / (111.32 * max(math.cos(math.radians((s + n) / 2)), 0.05))
        boxes.append((w - dlon, s - dlat, e + dlon, n + dlat))
    arr = np.array(boxes)
    for i in range(len(tags)):
        hit = np.flatnonzero(
            (arr[:, 0] <= arr[i, 2]) & (arr[:, 2] >= arr[i, 0])
            & (arr[:, 1] <= arr[i, 3]) & (arr[:, 3] >= arr[i, 1])
        )
        for j in hit[hit > i]:
            ri, rj = find(i), find(int(j))
            if ri != rj:
                parent[rj] = ri
    roots = {find(i) for i in range(len(tags))}
    gid = {r: k for k, r in enumerate(sorted(roots))}
    return {t: gid[find(i)] for i, t in enumerate(tags)}


def assign_folds(unit_key: dict[str, str], unit_sort: dict[str, tuple],
                 n_folds: int = N_FOLDS) -> dict[str, int]:
    """単位（場所グループ・国）を並べ替えて n_folds へ巡回割当（決定的）。

    unit_sort[u] は並べ替えキー（例: (大陸, −セル数, 単位 ID)）。
    """
    units = sorted(set(unit_key.values()), key=lambda u: unit_sort[u])
    fold_of_unit = {u: k % n_folds for k, u in enumerate(units)}
    return {t: fold_of_unit[u] for t, u in unit_key.items()}


def continent_map(country_info: Path) -> dict[str, str]:
    """GeoNames countryInfo.txt → {ISO2: 大陸コード}。"""
    out = {}
    for line in country_info.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        out[f[0]] = f[8]
    return out


def group_quotas(n_eval: dict[str, int], group: dict[str, int],
                 cap: int = GROUP_CAP) -> dict[str, int]:
    """場所グループごとに合計 cap 個まで、AOI へ評価セル数に比例して配分（最大剰余法、決定的）。"""
    members: dict[int, list[str]] = {}
    for a in sorted(n_eval):
        members.setdefault(group[a], []).append(a)
    out: dict[str, int] = {}
    for aois in members.values():
        tot = sum(n_eval[a] for a in aois)
        if tot <= cap:
            out |= {a: n_eval[a] for a in aois}
            continue
        exact = {a: cap * n_eval[a] / tot for a in aois}
        q = {a: math.floor(v) for a, v in exact.items()}
        rest = cap - sum(q.values())
        for a in sorted(aois, key=lambda a: (-(exact[a] - q[a]), a))[:rest]:
            q[a] += 1
        out |= q
    return out


def load_aoi_table(d: Path, *, cap: int, seed: int, unspecified_as_riverine: bool = False,
                   all_types: bool = False) -> pd.DataFrame:
    """1 AOI の評価セル（主ラベル付き）+ 特徴。cap を超えれば固定シードで抽出。"""
    cells = pq.read_table(d / "cells.parquet").to_pandas()
    feats = pq.read_table(d / "features.parquet").to_pandas()
    lc = pq.read_table(d / "landcover.parquet").to_pandas()  # 洪水前年の CCI（§2.3）
    df = pd.concat([cells.reset_index(drop=True), feats.reset_index(drop=True),
                    lc.reset_index(drop=True)], axis=1)
    y, keep = riverine_labels(df, unspecified_as_riverine=unspecified_as_riverine,
                              all_types=all_types)
    df["y"] = y
    df = df.loc[evaluable(df) & keep]
    if len(df) > cap:
        rng = np.random.default_rng([seed, _stable_hash(d.name)])
        df = df.iloc[np.sort(rng.choice(len(df), size=cap, replace=False))]
    df = df.reset_index(drop=True)
    df["aoi"] = d.name
    return df


def _stable_hash(s: str) -> int:
    h = 1469598103934665603
    for ch in s.encode("utf-8"):
        h = ((h ^ ch) * 1099511628211) % (1 << 63)
    return h


class CoastIndex:
    """Natural Earth 10 m 海岸線の頂点 KD 木（AOI 重心の海岸距離、§2.4 の内陸判定用）。"""

    def __init__(self, zip_path: Path) -> None:
        import shapely  # noqa: PLC0415
        from pyogrio import raw  # noqa: PLC0415
        from scipy.spatial import cKDTree  # noqa: PLC0415

        _, _, g, _ = raw.read(f"/vsizip/{zip_path}/ne_10m_coastline.shp")
        coords = shapely.get_coordinates(shapely.from_wkb(g))
        self._tree = cKDTree(_unit(coords[:, 0], coords[:, 1]))

    def distance_km(self, lon: float, lat: float) -> float:
        chord, _ = self._tree.query(_unit(np.array([lon]), np.array([lat]))[0])
        return float(2 * 6371.0 * np.arcsin(min(chord / 2, 1.0)))


def _unit(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lo, la = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
