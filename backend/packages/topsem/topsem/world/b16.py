"""amendment B-16（初見の査読 第 12 回、事後）の実装。

- C1: 最寄りの河道の上流の面積。HydroRIVERS の全ての河道の頂点（AOI の範囲 ± 0.2°）のうち、
  セルの中心に最も近い頂点の UPLAND_SKM の log10 と、その頂点までの距離（B-9 の C7 の河道からの
  高さと同じ球面の弦の距離）
- C5: EMS の詳細図（DETAIL の製品）を、概観図から作ったラベルの AOI の格子に同じ規則で落とす
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree

from topsem.world.b9 import EARTH_M, _unit_xyz
from topsem.world.ems import (
    FLOOD_CLASSES,
    GEN0_UNTYPED,
    LAYER_KEYS,
    Grid,
    _find_layer,
    _to_local,
    cover_fraction,
    flood_class,
    read_layer,
    schema_generation,
)


def in_box(lon: np.ndarray, lat: np.ndarray, box: tuple[float, float, float, float]) -> np.ndarray:
    """経緯度が範囲（w, s, e, n、両端を含む）に入るか。"""
    w, s, e, n = box
    return (lon >= w) & (lon <= e) & (lat >= s) & (lat <= n)


def nearest_upland(lon: np.ndarray, lat: np.ndarray, vlon: np.ndarray, vlat: np.ndarray,
                   vup: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """各点（セルの中心）に最も近い頂点の上流の面積の log10 と距離（m）。頂点が無ければ NaN。"""
    if len(vlon) == 0:
        nan = np.full(len(lon), np.nan, dtype=np.float32)
        return nan, nan.copy()
    tree = cKDTree(_unit_xyz(np.asarray(vlon, np.float64), np.asarray(vlat, np.float64)))
    d, j = tree.query(_unit_xyz(np.asarray(lon, np.float64), np.asarray(lat, np.float64)))
    up = np.asarray(vup, np.float64)[j]
    lupl = np.where(up > 0, np.log10(np.where(up > 0, up, 1.0)), np.nan)
    return lupl.astype(np.float32), (d * EARTH_M).astype(np.float32)


def merge_nearest(best_d: np.ndarray, best_v: np.ndarray, d: np.ndarray, v: np.ndarray) -> None:
    """二つの候補の集合（大陸ごとの河道）の最寄りを合わせる。

    best_* をその場で更新する。NaN の距離は候補が無いものとする。
    """
    take = np.isfinite(d) & (~np.isfinite(best_d) | (d < best_d))
    best_d[take] = d[take]
    best_v[take] = v[take]


# %% C5: EMS の詳細図を AOI の格子に落とす
def aoi_grid(meta: dict) -> Grid:
    """ラベルの meta.json の格子（LAEA、250 m、副格子）を作り直す。"""
    g = meta["grid"]
    sub = round(float(g["cell_m"]) / float(g["fine_res_m"]))
    return Grid(crs=CRS.from_proj4(meta["crs_proj4"]), x0=float(g["x0"]), y1=float(g["y1"]),
                nx=int(g["nx"]), ny=int(g["ny"]), sub=sub)


def _event_parts(z: Path, member: str, gen: int | None, tr: Transformer, acc: dict) -> None:
    """浸水の層を種別・痕跡・Not Analysed に分けて acc に足す（build_aoi_labels と同じ読み方）。"""
    ev = read_layer(z, member)
    n_ev = len(ev.geoms)
    al = {k.lower(): v for k, v in ev.attrs.items()}
    untyped = [GEN0_UNTYPED if gen == 0 else ""] * n_ev
    etype = [str(v) for v in
             (al.get("event_type") or al.get("subtype") or al.get("sbtypdes") or untyped)]
    desc = [str(v) for v in al.get("obj_desc", [""] * n_ev)]
    nota = [str(v) for v in (al.get("notation") or al.get("interpret") or [""] * n_ev)]
    hole = np.array(["not analys" in n.lower() for n in nota], dtype=bool)
    is_flood = np.array(
        [("flood" in t.lower() or t.startswith("5") or "flood" in n.lower()) and not h
         for t, n, h in zip(etype, nota, hole, strict=True)], dtype=bool)
    cls = np.array([flood_class(d) for d in desc], dtype=object)
    trace = np.array(["trace" in n.lower() for n in nota], dtype=bool)
    ev_l = _to_local(ev.geoms, tr) if n_ev else ev.geoms
    for c in FLOOD_CLASSES:
        acc[f"cls_{c}"] += list(ev_l[is_flood & (cls == c)])
    acc["trace"] += list(ev_l[is_flood & trace])
    acc["any"] += list(ev_l[is_flood])
    acc["hole"] += list(ev_l[hole])


def _footprint_parts(z: Path, member: str | None, aoi_l: list, tr: Transformer,
                     acc: dict) -> None:
    """撮像 footprint と Not Analysed。footprint の多角形の無い製品は AOI を観測域とする。"""
    img_any = False
    if member:
        fp = read_layer(z, member)
        fl = {k.lower(): v for k, v in fp.attrs.items()}
        otype = [str(v) for v in fl.get("obj_type", [""] * len(fp.geoms))]
        fp_l = _to_local(fp.geoms, tr)
        na = np.array(["not analys" in o.lower() for o in otype], dtype=bool)
        img = np.array(["footprint" in o.lower() for o in otype], dtype=bool)
        img_any = bool(img.any())
        acc["seen"] += list(fp_l[img])
        acc["hole"] += list(fp_l[na])
    if not img_any:
        acc["seen"] += aoi_l


def detail_fractions(zips: list[Path], grid: Grid) -> dict[str, np.ndarray]:
    """詳細図の製品（複数）の層を合わせて一つの製品として、AOI の格子に被覆率で落とす。

    層の読み方・浸水の種別・痕跡・Not Analysed・撮像 footprint・恒久水は build_aoi_labels と
    同じ規則。多角形は製品をまたいで合わせる（重なりは一つと数える）。戻り値は格子全体
    （ny*nx、行 0 = 北）と、製品の書式の世代。
    """
    tr = Transformer.from_crs(4326, grid.crs, always_xy=True)
    acc: dict[str, list] = {k: [] for k in ("aoi", "seen", "hole", "perm", "trace", "any",
                                            *(f"cls_{c}" for c in FLOOD_CLASSES))}
    gens = []
    for z in zips:
        with zipfile.ZipFile(z) as zf:
            names = zf.namelist()
        found = {k: _find_layer(names, k) for k in LAYER_KEYS}
        if not found["aoi"]:
            continue
        gen = schema_generation(found)
        gens.append(-1 if gen is None else gen)
        aoi_l = list(_to_local(read_layer(z, found["aoi"]).geoms, tr))
        acc["aoi"] += aoi_l
        if found["event"]:
            _event_parts(z, found["event"], gen, tr, acc)
        _footprint_parts(z, found["footprint"], aoi_l, tr, acc)
        if found["hydro"]:
            acc["perm"] += list(_to_local(read_layer(z, found["hydro"]).geoms, tr))
    seen = cover_fraction(acc["seen"], grid)
    out = {"f_aoi": cover_fraction(acc["aoi"], grid),
           "f_observed": np.clip(seen - cover_fraction(acc["hole"], grid), 0.0, 1.0),
           "f_perm_water": cover_fraction(acc["perm"], grid),
           "f_trace": cover_fraction(acc["trace"], grid),
           "f_flood_any": cover_fraction(acc["any"], grid)}
    out |= {f"f_{c}": cover_fraction(acc[f"cls_{c}"], grid) for c in FLOOD_CLASSES}
    out["generations"] = np.array(gens, dtype=np.int64)
    return out
