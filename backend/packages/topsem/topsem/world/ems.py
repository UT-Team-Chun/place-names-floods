"""Copernicus EMS delineation 製品 → AOI ローカル 250 m 格子の観測浸水ラベル。

ラベルの定義（protocol v2 で凍結する前提の実装。値の閾値はここでは適用せず、
セルごとの被覆率として保存する）:

- 評価領域 = AOI ∩ 撮像 footprint −「Not Analysed」（雲・欠測）
- 浸水 = `observedEventA` のうち event_type が洪水のもの。種別は**地物属性の
  obj_desc**（Riverine / Flash / Coastal …）で分ける。activation 単位の
  subCategory は誤りがある（EMSR772 は内陸河川なのに "Coastal flood"）ので使わない
- 恒常水域 = `hydrographyA`（参照データの面状水域）
- floodDepthA 等の DEM 由来層は**使わない**（地形特徴との循環）

格子は AOI 重心を中心とする LAEA（m 単位）で 250 m。被覆率は 25 m（大 AOI は
50 m）の副格子でラスタ化してブロック平均する。
"""

from __future__ import annotations

import hashlib
import json
import re
import urllib.request
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow.parquet as pq
import shapely
from affine import Affine
from pyogrio import raw
from pyproj import CRS, Transformer
from rasterio import features
from rasterio.transform import from_origin

CELL_M = 250.0
MAX_FINE_CELLS = 250_000  # これを超える AOI は副格子 50 m（メモリ上限）
UA = {"User-Agent": "toponym-semantics research (ryoutaokauchi@g.ecc.u-tokyo.ac.jp)"}

# 層名は 4 世代で綴りが違う（2026-10-25 に全世代の実物で確認）:
#   第 0 世代（2012–2013）: *_AreaOfInterest.shp / *_Crisis_event_A.shp（SBTYPdes
#                        "EM009 - Flood"、種別の記載なし。intrp は符号表が公開されていないので
#                        使わない）/ *_Hydrology_A_LAF.shp（*_Hydrology_A_Optional_LAF には島が
#                        入るので使わない）、撮像 footprint 層なし
#   第 1 世代（〜2017）: *_area_of_interest.shp / *_crisis_information_poly.shp（洪水は
#                        subtype "EM009 - Flood"・interpret "Flooded Area"、種別の記載なし）/
#                        *_hydrography_poly.shp、撮像 footprint 層なし
#   第 2 世代（2018–2021 頃）: *_area_of_interest_a.shp / *_observed_event_a.shp /
#                        *_sensor_metadata_a
#   第 3 世代（新系）: *_areaOfInterestA_v1.shp / *_observedEventA_v2.shp / *_imageFootprintA
LAYER_KEYS: dict[str, tuple[str, ...]] = {
    "aoi": ("areaofinteresta", "area_of_interest", "areaofinterest."),
    "event": ("observedeventa", "observed_event_a", "crisis_information_poly", "crisis_event_a."),
    "footprint": ("imagefootprinta", "image_footprint_a", "sensor_metadata_a"),
    "hydro": ("hydrographya", "hydrography_a", "hydrography_poly", "hydrology_a_laf."),
}
# 第 0 世代で種別欄（SBTYPdes）を持たない浸水層の多角形（洪水 activation の Crisis_event_A）
GEN0_UNTYPED = "flood (generation-0 crisis event, untyped)"


def schema_generation(found: dict[str, str | None]) -> int | None:
    """AOI 層の綴りから書式の世代（0–3）を判定（来歴と感度解析用）。"""
    a = (found.get("aoi") or "").lower()
    if "areaofinteresta" in a:
        return 3
    if "area_of_interest_a" in a:
        return 2
    if "area_of_interest" in a:
        return 1
    if "areaofinterest." in a:
        return 0
    return None
FLOOD_CLASSES = ("riverine", "flash", "coastal", "other")


def flood_class(obj_desc: str) -> str:
    """observedEventA.obj_desc → 種別。"""
    s = (obj_desc or "").lower()
    if "river" in s or "fluvial" in s:
        return "riverine"
    if "flash" in s or "torrent" in s:
        return "flash"
    if "coast" in s or "storm surge" in s or "tidal" in s:
        return "coastal"
    return "other"


def download(url: str, dst: Path) -> str:
    """zip を取得してキャッシュ（新系は署名付き S3 への 302 を urllib が追従）。sha256 を返す。

    壊れた/zip でない応答（エラーページ等）はキャッシュに残さず 1 回だけ取り直す。
    それでも zip でなければ ValueError（台帳にあるが配信されていない製品）。
    """
    if dst.is_file() and not zipfile.is_zipfile(dst):
        dst.unlink()
    for attempt in range(2):
        _fetch(url, dst)
        if zipfile.is_zipfile(dst):
            break
        dst.unlink()
        if attempt == 1:
            msg = f"product not a zip (unavailable): {url}"
            raise ValueError(msg)
    return hashlib.sha256(dst.read_bytes()).hexdigest()


def _fetch(url: str, dst: Path) -> None:
    if not dst.is_file():
        dst.parent.mkdir(parents=True, exist_ok=True)
        headers = {**UA}
        if "rapidmapping.emergency.copernicus.eu" in url:
            headers |= {
                "Origin": "https://mapping.emergency.copernicus.eu",
                "Referer": "https://mapping.emergency.copernicus.eu/",
            }
        req = urllib.request.Request(url, headers=headers)  # noqa: S310
        tmp = dst.with_suffix(".part")
        with urllib.request.urlopen(req, timeout=300) as r, tmp.open("wb") as f:  # noqa: S310
            while chunk := r.read(1 << 20):
                f.write(chunk)
        tmp.rename(dst)


def _find_layer(names: list[str], key: str) -> str | None:
    """zip 内の .shp から層キーに合う 1 本（無ければ .json）。"""
    pats = LAYER_KEYS[key]
    for ext in (".shp", ".json"):
        for n in names:
            low = n.lower()
            if low.endswith(ext) and any(p in low for p in pats):
                return n
    return None


@dataclass
class Layer:
    geoms: np.ndarray  # shapely geometry array (EPSG:4326)
    attrs: dict[str, list[Any]] = field(default_factory=dict)


def read_layer(zip_path: Path, member: str) -> Layer:
    """zip 内の 1 層を WGS84 の shapely 配列 + 属性として読む。"""
    vsi = f"/vsizip/{zip_path}/{member}"
    meta, _, geom, fields = raw.read(vsi)
    geoms = shapely.from_wkb(geom)
    crs = meta.get("crs")
    if crs and CRS.from_user_input(crs) != CRS.from_epsg(4326):
        tr = Transformer.from_crs(crs, 4326, always_xy=True)
        geoms = shapely.transform(
            geoms, lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1])))
    attrs = {str(k): list(v) for k, v in zip(meta["fields"], fields, strict=True)}
    return Layer(geoms=shapely.make_valid(geoms), attrs=attrs)


def laea_proj4(lon0: float, lat0: float) -> str:
    return f"+proj=laea +lat_0={lat0:.6f} +lon_0={lon0:.6f} +datum=WGS84 +units=m +no_defs"


def laea(lon0: float, lat0: float) -> CRS:
    return CRS.from_proj4(laea_proj4(lon0, lat0))


def _to_local(geoms: np.ndarray, tr: Transformer) -> np.ndarray:
    return shapely.transform(geoms, lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1])))


@dataclass
class Grid:
    crs: CRS
    x0: float  # 左上 x（m）
    y1: float  # 左上 y（m）
    nx: int
    ny: int
    sub: int  # 1 セルあたりの副格子数（一辺）

    @property
    def fine_res(self) -> float:
        return CELL_M / self.sub

    def fine_transform(self) -> Affine:
        return from_origin(self.x0, self.y1, self.fine_res, self.fine_res)

    def cell_centres_local(self) -> tuple[np.ndarray, np.ndarray]:
        ix, iy = np.meshgrid(np.arange(self.nx), np.arange(self.ny))
        x = self.x0 + (ix + 0.5) * CELL_M
        y = self.y1 - (iy + 0.5) * CELL_M
        return x.ravel(), y.ravel()


def make_grid(aoi_local: shapely.Geometry, crs: CRS) -> Grid:
    minx, miny, maxx, maxy = aoi_local.bounds
    x0 = np.floor(minx / CELL_M) * CELL_M
    y1 = np.ceil(maxy / CELL_M) * CELL_M
    nx = int(np.ceil((maxx - x0) / CELL_M))
    ny = int(np.ceil((y1 - miny) / CELL_M))
    sub = 10 if nx * ny <= MAX_FINE_CELLS else 5
    return Grid(crs=crs, x0=float(x0), y1=float(y1), nx=nx, ny=ny, sub=sub)


def cover_fraction(geoms: list, grid: Grid) -> np.ndarray:
    """多角形群の被覆率を 250 m セルごとに（[ny*nx]、行 0 = 北）。"""
    shape = (grid.ny * grid.sub, grid.nx * grid.sub)
    if not geoms:
        return np.zeros(grid.ny * grid.nx, dtype=np.float32)
    fine = features.rasterize(
        ((g, 1) for g in geoms if g is not None and not g.is_empty),
        out_shape=shape,
        transform=grid.fine_transform(),
        fill=0,
        dtype="uint8",
    )
    blk = fine.reshape(grid.ny, grid.sub, grid.nx, grid.sub).sum(axis=(1, 3), dtype=np.uint16)
    return (blk.astype(np.float32) / float(grid.sub * grid.sub)).ravel()


def build_aoi_labels(  # noqa: PLR0915 — 層ごとの手順を 1 か所に並べて読める形を優先
    zip_path: Path,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """1 製品 zip → セル列（dict of arrays）と来歴メタデータ。"""
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
    found = {k: _find_layer(names, k) for k in LAYER_KEYS}
    if not found["aoi"]:
        msg = f"AOI 層が無い: {zip_path.name} {found}"
        raise ValueError(msg)
    # observedEvent 層の無い製品 = 浸水を描画しなかった AOI（撮像時に水が引いていた等）。
    # 2015–2018 の旧世代に 300 件超（2026-10-25 実測）。全被覆 0 のラベルとして残し、
    # 採用基準（浸水セル ≥ 20）で除外される経路を台帳に残す（分母を誠実に保つ）
    aoi = read_layer(zip_path, found["aoi"])
    aoi_geom = shapely.union_all(aoi.geoms)
    c = aoi_geom.centroid
    crs = laea(c.x, c.y)
    tr = Transformer.from_crs(4326, crs, always_xy=True)
    aoi_l = shapely.union_all(_to_local(aoi.geoms, tr))
    grid = make_grid(aoi_l, crs)

    if found["event"]:
        ev = read_layer(zip_path, found["event"])
    else:
        ev = Layer(geoms=np.array([], dtype=object), attrs={})
    n_ev = len(ev.geoms)
    gen = schema_generation(found)
    # 属性名の大文字小文字は製品で揺れる（第 0 世代 EMSR009 は大文字）ので小文字で引く。
    # 第 1 世代は event_type→subtype、notation→interpret、第 0 世代は SBTYPdes。obj_desc が
    # 無い世代は種別不明 = "other"（§2.4 の内陸規則）
    al = {k.lower(): v for k, v in ev.attrs.items()}
    untyped = [GEN0_UNTYPED if gen == 0 else ""] * n_ev
    etype = [str(v) for v in
             (al.get("event_type") or al.get("subtype") or al.get("sbtypdes") or untyped)]
    desc = [str(v) for v in al.get("obj_desc", [""] * n_ev)]
    nota = [str(v) for v in (al.get("notation") or al.get("interpret") or [""] * n_ev)]
    ev_hole = np.array(["not analys" in n.lower() for n in nota], dtype=bool)
    is_flood = np.array(
        [("flood" in t.lower() or t.startswith("5") or "flood" in n.lower()) and not h
         for t, n, h in zip(etype, nota, ev_hole, strict=True)], dtype=bool)
    cls = np.array([flood_class(d) for d in desc], dtype=object)
    trace = np.array(["trace" in n.lower() for n in nota], dtype=bool)
    ev_l = _to_local(ev.geoms, tr) if len(ev.geoms) else ev.geoms

    cols: dict[str, np.ndarray] = {}
    for k in FLOOD_CLASSES:
        cols[f"f_{k}"] = cover_fraction(list(ev_l[is_flood & (cls == k)]), grid)
    cols["f_trace"] = cover_fraction(list(ev_l[is_flood & trace]), grid)
    cols["f_flood_any"] = cover_fraction(list(ev_l[is_flood]), grid)
    cols["f_aoi"] = cover_fraction([aoi_l], grid)

    fp_meta: dict[str, Any] = {"footprint_layer": found["footprint"]}
    if found["footprint"]:
        fp = read_layer(zip_path, found["footprint"])
        fl = {k.lower(): v for k, v in fp.attrs.items()}
        otype = [str(v) for v in fl.get("obj_type", [""] * len(fp.geoms))]
        fp_l = _to_local(fp.geoms, tr)
        na = np.array(["not analys" in o.lower() for o in otype])
        img = np.array(["footprint" in o.lower() for o in otype])
        # 撮像 footprint 多角形を持たない製品がある（EMSR772 AOI02: 'Not Analysed' のみ）。
        # その場合は「AOI − Not Analysed」を観測域とし、flag を残す（感度分析で除外可能）
        implicit = not bool(img.any())
        seen = cols["f_aoi"].copy() if implicit else cover_fraction(list(fp_l[img]), grid)
        hole = cover_fraction(list(fp_l[na]), grid)
        cols["f_observed"] = np.clip(seen - hole, 0.0, 1.0)
        fp_meta |= {"n_footprint": int(img.sum()), "n_not_analysed": int(na.sum()),
                    "footprint_implicit": implicit}
    else:
        cols["f_observed"] = cols["f_aoi"].copy()
        fp_meta |= {"footprint_missing": True}
    if ev_hole.any():  # 第 1 世代は Not Analysed を浸水層の interpret に持つことがある
        cols["f_observed"] = np.clip(
            cols["f_observed"] - cover_fraction(list(ev_l[ev_hole]), grid), 0.0, 1.0)
    if found["hydro"]:
        hy = read_layer(zip_path, found["hydro"])
        cols["f_perm_water"] = cover_fraction(list(_to_local(hy.geoms, tr)), grid)
    else:
        cols["f_perm_water"] = np.zeros(grid.nx * grid.ny, dtype=np.float32)

    x, y = grid.cell_centres_local()
    inv = Transformer.from_crs(crs, 4326, always_xy=True)
    lon, lat = inv.transform(x, y)
    ix, iy = np.meshgrid(np.arange(grid.nx), np.arange(grid.ny))
    cols |= {
        "ix": ix.ravel().astype(np.int32),
        "iy": iy.ravel().astype(np.int32),
        "x": x.astype(np.float32),
        "y": y.astype(np.float32),
        "lon": np.asarray(lon, dtype=np.float64),
        "lat": np.asarray(lat, dtype=np.float64),
    }
    keep = cols["f_aoi"] > 0
    cols = {k: v[keep] for k, v in cols.items()}
    counts = {
        "obj_desc": dict(zip(*map(list, np.unique(desc, return_counts=True)), strict=True))
        if desc else {},
        "notation": dict(zip(*map(list, np.unique(nota, return_counts=True)), strict=True))
        if nota else {},
        "event_type": {str(k): int(v) for k, v in
                       zip(*np.unique(etype, return_counts=True), strict=True)} if etype else {},
        "n_event_polygons": len(ev.geoms),
        "observed_event_layer": bool(found["event"]),
    }
    counts["obj_desc"] = {str(k): int(v) for k, v in counts["obj_desc"].items()}
    counts["notation"] = {str(k): int(v) for k, v in counts["notation"].items()}
    meta = {
        "layers": found,
        "schema_generation": gen,
        "crs_proj4": laea_proj4(c.x, c.y),
        "grid": {"x0": grid.x0, "y1": grid.y1, "nx": grid.nx, "ny": grid.ny,
                 "cell_m": CELL_M, "fine_res_m": grid.fine_res},
        "aoi_area_km2": float(aoi_l.area / 1e6),
        "n_cells": int(keep.sum()),
        **counts,
        **fp_meta,
    }
    return cols, meta


def cells_digest(path: Path) -> str:
    """cells.parquet の内容の指紋（列名順に列名と値のバイト列を sha256）。

    parquet の書き込みメタデータに依らず、ラベル再構築で中身が変わったかを比べるために使う。
    """
    tab = pq.read_table(path)
    h = hashlib.sha256()
    for name in sorted(tab.column_names):
        h.update(name.encode("utf-8"))
        h.update(np.ascontiguousarray(tab.column(name).to_numpy()).tobytes())
    return h.hexdigest()


def select_products(del_products: list[dict]) -> list[dict]:
    """AOI ごとに 1 製品: MONIT 以外・広域（OVERVIEW/MAP/PRODUCT 優先、DETAIL は代替）・最新版。"""
    rank = {"PRODUCT": 0, "OVERVIEW": 0, "MAP": 0, "DETAIL": 1}
    by_aoi: dict[str, list[dict]] = {}
    for p in del_products:
        if p.get("kind") != "PRODUCT":
            continue
        by_aoi.setdefault(p["aoi"], []).append(p)
    out = []
    for aoi, ps in sorted(by_aoi.items()):
        best = min(
            ps,
            key=lambda p: (rank.get(p.get("scale", "PRODUCT"), 1), -int(p.get("version", 0)),
                           _detail_no(p["url"])),
        )
        out.append({**best, "aoi": aoi})
    return out


def _detail_no(url: str) -> int:
    m = re.search(r"DETAIL(\d+)", url)
    return int(m.group(1)) if m else 0


def write_meta(path: Path, meta: dict) -> None:
    path.write_text(json.dumps(meta, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
