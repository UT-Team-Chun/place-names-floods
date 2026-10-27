"""Copernicus GFM（全球洪水監視）の観測浸水域を AOI の 250 m セルのラベルに（事後、B-3）。

GFM は Sentinel-1（IW GRD）の全撮像を 3 つの手法（LIST・DLR・TUW）で浸水判定し、多数決で
合わせた観測浸水域（20 m、Equi7Grid、2015 年以降）。EMS とは独立に作られた全球の浸水図として、
EMS のラベルの点検と、名前の増分の再現に使う。出所は EODC の公開 STAC（collection ``GFM``）と、
そこに記載された COG（認証なし）。値は GFM Product Definition Document v4.0 のとおり:
浸水域 0 = 浸水なし・1 = 浸水、
除外域 0 = 判定可能・それ以外 = 除外（市街地・密な植生・水に似た面・地形の歪みと HAND 10 m 以上・
レーダーの影・観測不足）、参照水域 0 = 通常は水なし・1 = 恒常水域・2 = その月の季節水域。いずれも
255 = データなし。

ラベルの規則（B-3 で計算の前に固定）:

- 撮像: EMS の事後画像の日付の最小 − 1 日から最大 + 1 日まで（UTC の日付）に取得され、
  AOI に掛かる全項目。同じ取得時刻の項目（隣り合うタイル）は 1 つの撮像
- EMS と同じ副格子（25 m、大きな AOI は 50 m）に最近傍で写し、浸水域 ∈ {0, 1}・除外域 = 0・
  参照水域 = 0 の画素を判定可能とする
- 撮像ごとに、セルの判定可能な割合 v と、判定可能な部分に占める浸水の割合 f
- どれかの撮像で v ≥ 0.5 なら評価可能、v ≥ 0.5 の撮像のどれかで f ≥ 0.5 なら浸水
- AOI の採用: EMS で評価可能なセルのうち GFM でも評価可能なものが 400 以上、そのうち浸水が
  20 以上、浸水率が 0.8 以下（確証解析の採用基準と同じ閾値）
"""

from __future__ import annotations

import datetime as dt
import json
import random
import time
import urllib.parse
import urllib.request
from typing import Any

import numpy as np
import rasterio
from pyproj import CRS
from rasterio.enums import Resampling
from rasterio.errors import RasterioIOError
from rasterio.warp import reproject

from topsem.world.landcover import GDAL_ENV
from topsem.world.terrain import _block, grid_lonlat_bounds

STAC_ITEMS = "https://stac.eodc.eu/api/v1/collections/GFM/items"
LAYERS = ("ensemble_flood_extent", "exclusion_mask", "reference_water_mask")
NODATA = 255
WINDOW_DAYS = 1
V_MIN = 0.5  # セルの判定可能な割合の下限
F_MIN = 0.5  # 判定可能な部分の浸水の割合の下限（EMS の WET_TH と同じ）
MIN_EVAL, MIN_WET, MAX_WET_FRAC = 400, 20, 0.8  # 確証解析の採用基準と同じ閾値
DATE_FORMATS = ("%d/%m/%Y", "%Y/%m/%d")
RETRIES = 8
MAX_PAGES = 50
# EMS の出典表に現れる SAR 衛星の名（小文字の部分一致）。Sentinel-1 は別に数える
SAR_NAMES = ("radarsat", "rcm", "cosmo", "terrasar", "tandem", "paz", "iceye", "alos", "kompsat-5",
             "saocom", "capella", "risat", "gaofen-3", "umbra")


def sensor_group(sources: list[dict]) -> str:
    """EMS の事後画像の種類: s1（Sentinel-1 を含む）・sar（他の SAR だけ）・optical・none。"""
    names = [str(s.get("name", "")).lower() for s in sources
             if "post" in str(s.get("phase", "")).lower()]
    if not names:
        return "none"
    if any("sentinel-1" in n or "sentinel 1" in n for n in names):
        return "s1"
    if any(k in n for n in names for k in SAR_NAMES):
        return "sar"
    return "optical"


def post_event_window(sources: list[dict], days: int = WINDOW_DAYS,
                      ) -> tuple[dt.date, dt.date] | None:
    """EMS 製品の出典表（v38）の事後画像の日付から、GFM を探す日付の範囲（両端を含む）。"""
    dates = []
    for s in sources:
        if "post" not in str(s.get("phase", "")).lower():
            continue
        for fmt in DATE_FORMATS:
            try:
                text = str(s.get("date", "")).strip()
                d = dt.datetime.strptime(text, fmt).replace(tzinfo=dt.UTC)
            except ValueError:
                continue
            dates.append(d.date())
            break
    if not dates:
        return None
    return min(dates) - dt.timedelta(days=days), max(dates) + dt.timedelta(days=days)


def _get_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "toponym-semantics research"})  # noqa: S310 — 固定の https
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
                return json.loads(r.read().decode("utf-8"))
        except OSError:
            if attempt == RETRIES - 1:
                raise
            time.sleep(min(120.0, 4.0 * 2**attempt) + random.uniform(0, 3))  # noqa: S311 — 待ち時間の揺らぎのみ
    msg = "unreachable"
    raise RuntimeError(msg)


def search_items(bounds: tuple[float, float, float, float], d0: dt.date, d1: dt.date,
                 ) -> list[dict]:
    """AOI の lon/lat 範囲に掛かり、d0 から d1 まで（UTC の日付、両端を含む）に取得された項目。"""
    q = urllib.parse.urlencode({
        "bbox": ",".join(f"{v:.6f}" for v in bounds),
        "datetime": f"{d0.isoformat()}T00:00:00Z/{d1.isoformat()}T23:59:59Z",
        "limit": 200,
    })
    url, feats = f"{STAC_ITEMS}?{q}", []
    for _ in range(MAX_PAGES):
        page = _get_json(url)
        feats += page.get("features", [])
        nxt = [lk["href"] for lk in page.get("links", []) if lk.get("rel") == "next"]
        if not nxt:
            break
        url = nxt[0]
    return sorted(feats, key=lambda f: (f["properties"]["datetime"], f["id"]))


def group_scenes(items: list[dict]) -> dict[str, list[dict]]:
    """同じ取得時刻の項目（隣り合うタイル）を 1 つの撮像にまとめる。"""
    out: dict[str, list[dict]] = {}
    for f in items:
        out.setdefault(f["properties"]["datetime"], []).append(f)
    return out


def valid_and_wet(flood: np.ndarray, excl: np.ndarray, refw: np.ndarray,
                  ) -> tuple[np.ndarray, np.ndarray]:
    """副格子の画素: 判定可能（浸水域 0/1・除外域 0・参照水域 0）と、判定可能かつ浸水。"""
    valid = ((flood == 0) | (flood == 1)) & (excl == 0) & (refw == 0)
    return valid, valid & (flood == 1)


def cell_fractions(valid: np.ndarray, wet: np.ndarray, sub: int, ny: int, nx: int,
                   ) -> tuple[np.ndarray, np.ndarray]:
    """副格子 → 250 m セル（行 0 = 北）: 判定可能な割合 v と、判定可能な部分の浸水の割合 f。"""
    nv = _block(valid, sub, ny, nx).sum(axis=2)
    nw = _block(wet, sub, ny, nx).sum(axis=2)
    v = nv / float(sub * sub)
    f = np.divide(nw, nv, out=np.zeros(nv.shape, dtype=np.float64), where=nv > 0)
    return v.astype(np.float32).ravel(), f.astype(np.float32).ravel()


def combine_scenes(vs: list[np.ndarray], fs: list[np.ndarray],
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """撮像ごとの (v, f) → セルの最大の v・評価可能・ラベル（評価不能は −1）・評価可能な撮像の数。

    評価可能 = どれかの撮像で v ≥ 0.5。浸水 = v ≥ 0.5 の撮像のどれかで f ≥ 0.5。
    """
    v = np.stack(vs)
    ok = v >= V_MIN
    ev = ok.any(axis=0)
    wet = (ok & (np.stack(fs) >= F_MIN)).any(axis=0)
    y = np.where(ev, wet.astype(np.int8), np.int8(-1)).astype(np.int8)
    return v.max(axis=0), ev, y, ok.sum(axis=0).astype(np.int8)


def _item_layers(item: dict, shape: tuple[int, int], dst_tr: rasterio.Affine, dst_crs: CRS,
                 ) -> dict[str, np.ndarray]:
    """1 項目（タイル）の 3 層を副格子へ最近傍で。

    どれかの層が読めなければ例外（呼び出し側で項目ごと外す）。
    """
    out = {}
    for k in LAYERS:
        with rasterio.open(f"/vsicurl/{item['assets'][k]['href']}") as src:
            tmp = np.full(shape, NODATA, dtype=np.uint8)
            reproject(source=rasterio.band(src, 1), destination=tmp, dst_transform=dst_tr,
                      dst_crs=dst_crs, src_nodata=NODATA, dst_nodata=NODATA,
                      resampling=Resampling.nearest)
            out[k] = tmp
    return out


def scene_fractions(items: list[dict], meta: dict,
                    ) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """1 つの撮像（タイルの組）の 3 層を AOI の副格子へ最近傍で写し、セルの v と f。

    層が読めない項目（例: 別の大陸格子のタイルで参照水域のファイルが無く HTTP 403）は使わず、
    その ID を返す（B-3 追記 1）。
    """
    g = meta["grid"]
    ny, nx = g["ny"], g["nx"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    shape = (ny * sub, nx * sub)
    dst_tr = rasterio.transform.from_origin(g["x0"], g["y1"], g["fine_res_m"], g["fine_res_m"])
    dst_crs = CRS.from_proj4(meta["crs_proj4"])
    layers = {k: np.full(shape, NODATA, dtype=np.uint8) for k in LAYERS}
    skipped = []
    with rasterio.Env(**GDAL_ENV):
        for item in items:
            try:
                got = _item_layers(item, shape, dst_tr, dst_crs)
            except (RasterioIOError, KeyError):  # 読めない層・目録に無い資産（B-3 追記 1）
                skipped.append(item["id"])
                continue
            for k in LAYERS:
                fill = (layers[k] == NODATA) & (got[k] != NODATA)
                layers[k][fill] = got[k][fill]
    valid, wet = valid_and_wet(*(layers[k] for k in LAYERS))
    v, f = cell_fractions(valid, wet, sub, ny, nx)
    return v, f, skipped


def item_provenance(item: dict) -> dict[str, Any]:
    """来歴に残す項目の ID・版・各層の URL（目録に無い層は None。B-3 追記 1 の補足）。"""
    return {"id": item["id"], "version": item["properties"].get("processing:version"),
            "hrefs": {k: item["assets"].get(k, {}).get("href") for k in LAYERS}}


def gfm_labels(meta: dict, sources: list[dict], ix: np.ndarray, iy: np.ndarray,
               ) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """1 AOI のセル（格子の ix, iy）の GFM ラベルと来歴。撮像が無ければ全セル評価不能。"""
    n = len(ix)
    cols = {"gfm_v": np.zeros(n, dtype=np.float32), "gfm_eval": np.zeros(n, dtype=bool),
            "gfm_y": np.full(n, -1, dtype=np.int8), "gfm_n_scenes": np.zeros(n, dtype=np.int8)}
    win = post_event_window(sources)
    if win is None:
        return cols, {"reason": "no post-event date", "window": None, "scenes": []}
    prov: dict[str, Any] = {"window": [win[0].isoformat(), win[1].isoformat()], "scenes": []}
    scenes = group_scenes(search_items(grid_lonlat_bounds(meta), *win))
    if not scenes:
        return cols, prov | {"reason": "no GFM item in window"}
    vs, fs = [], []
    for when, items in scenes.items():
        v, f, skipped = scene_fractions(items, meta)
        vs.append(v)
        fs.append(f)
        prov["scenes"].append({"datetime": when, "skipped_unreadable": skipped,
                               "items": [item_provenance(it) for it in items]})
    vmax, ev, y, nsc = combine_scenes(vs, fs)
    idx = iy.astype(np.int64) * int(meta["grid"]["nx"]) + ix.astype(np.int64)
    cols = {"gfm_v": vmax[idx], "gfm_eval": ev[idx], "gfm_y": y[idx], "gfm_n_scenes": nsc[idx]}
    return cols, prov | {"reason": "ok", "n_scenes": len(scenes)}


# %% 事後 B-6 の 4: 層の成分と、季節水域を判定可能に戻した変種 V1
LAYERS_V1 = (*LAYERS, "ensemble_water_extent")
REFW_PERMANENT, REFW_SEASONAL = 1, 2  # 参照水域の値
COMPONENTS = ("excl", "perm", "seas", "nodata")


def valid_and_wet_v1(water: np.ndarray, excl: np.ndarray, refw: np.ndarray,
                     ) -> tuple[np.ndarray, np.ndarray]:
    """変種 V1 の副格子の画素。

    判定可能 = 観測水域 0/1・除外域 0・参照水域 ≠ 1（恒常水域でない。季節水域 2 は判定可能に）。
    浸水 = 判定可能かつ観測水域 = 1（季節水域でも水が観測されれば浸水と数える）。
    """
    valid = (((water == 0) | (water == 1)) & (excl == 0) & (refw != REFW_PERMANENT)
             & (refw != NODATA))
    return valid, valid & (water == 1)


def component_fractions(excl: np.ndarray, refw: np.ndarray, sub: int, ny: int, nx: int,
                        ) -> dict[str, np.ndarray]:
    """副格子 → セル: 除外域・恒常水域・季節水域・データなし（その撮像が覆わない）の画素の割合。

    重なりは除外域を優先する。浸水域の層は恒常水域で値を持たない（255）ので、データなしは除外域の
    層で決める。
    """
    nodata = excl == NODATA
    ex = (excl != 0) & ~nodata
    perm = (refw == REFW_PERMANENT) & ~ex & ~nodata
    seas = (refw == REFW_SEASONAL) & ~ex & ~nodata
    tot = float(sub * sub)
    return {k: (_block(m, sub, ny, nx).sum(axis=2) / tot).astype(np.float32).ravel()
            for k, m in zip(COMPONENTS, (ex, perm, seas, nodata), strict=True)}


def scene_fractions_v1(items: list[dict], meta: dict) -> dict[str, Any]:
    """1 撮像: B-3 の (v, f)、変種 V1 の (v1, f1)、層の成分の割合。

    読めない項目は外す（B-3 追記 1 と同じ）。
    """
    g = meta["grid"]
    ny, nx = g["ny"], g["nx"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    shape = (ny * sub, nx * sub)
    dst_tr = rasterio.transform.from_origin(g["x0"], g["y1"], g["fine_res_m"], g["fine_res_m"])
    dst_crs = CRS.from_proj4(meta["crs_proj4"])
    layers = {k: np.full(shape, NODATA, dtype=np.uint8) for k in LAYERS_V1}
    skipped = []
    with rasterio.Env(**GDAL_ENV):
        for item in items:
            try:
                got = {}
                for k in LAYERS_V1:
                    with rasterio.open(f"/vsicurl/{item['assets'][k]['href']}") as src:
                        tmp = np.full(shape, NODATA, dtype=np.uint8)
                        reproject(source=rasterio.band(src, 1), destination=tmp,
                                  dst_transform=dst_tr, dst_crs=dst_crs, src_nodata=NODATA,
                                  dst_nodata=NODATA, resampling=Resampling.nearest)
                        got[k] = tmp
            except (RasterioIOError, KeyError):
                skipped.append(item["id"])
                continue
            for k in LAYERS_V1:
                fill = (layers[k] == NODATA) & (got[k] != NODATA)
                layers[k][fill] = got[k][fill]
    flood, excl, refw, water = (layers[k] for k in LAYERS_V1)
    valid, wet = valid_and_wet(flood, excl, refw)
    v, f = cell_fractions(valid, wet, sub, ny, nx)
    valid1, wet1 = valid_and_wet_v1(water, excl, refw)
    v1, f1 = cell_fractions(valid1, wet1, sub, ny, nx)
    comp = component_fractions(excl, refw, sub, ny, nx)
    return {"v": v, "f": f, "v1": v1, "f1": f1, "skipped": skipped,
            **{f"c_{k}": a for k, a in comp.items()}}


def gfm_labels_v1(meta: dict, sources: list[dict], ix: np.ndarray, iy: np.ndarray,
                  ) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """1 AOI のセルの B-3 のラベル（作り直しの照合用）・変種 V1・層の成分。

    層の成分は、B-3 の判定可能な割合 v が最大の撮像の値。
    """
    n = len(ix)
    cols: dict[str, np.ndarray] = {
        "gfm_eval": np.zeros(n, dtype=bool), "gfm_y": np.full(n, -1, dtype=np.int8),
        "v1_v": np.zeros(n, dtype=np.float32), "v1_eval": np.zeros(n, dtype=bool),
        "v1_y": np.full(n, -1, dtype=np.int8),
        **{f"c_{k}": np.full(n, np.nan, dtype=np.float32) for k in COMPONENTS}}
    win = post_event_window(sources)
    if win is None:
        return cols, {"reason": "no post-event date", "window": None, "scenes": []}
    prov: dict[str, Any] = {"window": [win[0].isoformat(), win[1].isoformat()], "scenes": []}
    scenes = group_scenes(search_items(grid_lonlat_bounds(meta), *win))
    if not scenes:
        return cols, prov | {"reason": "no GFM item in window"}
    res = []
    for when, items in scenes.items():
        r = scene_fractions_v1(items, meta)
        res.append(r)
        prov["scenes"].append({"datetime": when, "skipped_unreadable": r["skipped"],
                               "items": [item_provenance(it) for it in items]})
    _, ev, y, _ = combine_scenes([r["v"] for r in res], [r["f"] for r in res])
    v1max, ev1, y1, _ = combine_scenes([r["v1"] for r in res], [r["f1"] for r in res])
    best = np.stack([r["v"] for r in res]).argmax(axis=0)
    comp = {k: np.stack([r[f"c_{k}"] for r in res])[best, np.arange(best.size)]
            for k in COMPONENTS}
    idx = iy.astype(np.int64) * int(meta["grid"]["nx"]) + ix.astype(np.int64)
    cols = {"gfm_eval": ev[idx], "gfm_y": y[idx], "v1_v": v1max[idx], "v1_eval": ev1[idx],
            "v1_y": y1[idx], **{f"c_{k}": comp[k][idx] for k in COMPONENTS}}
    return cols, prov | {"reason": "ok", "n_scenes": len(scenes)}


def gfm_include(ems_eval: np.ndarray, gfm_eval: np.ndarray, gfm_y: np.ndarray,
                ) -> tuple[bool, dict[str, Any]]:
    """AOI の採用（B-3）: EMS と GFM の両方で評価可能なセルが 400 以上、浸水 20 以上。

    浸水率は 0.8 以下。閾値は確証解析の採用基準と同じ。
    """
    m = ems_eval & gfm_eval
    n = int(m.sum())
    w = int((gfm_y[m] == 1).sum())
    frac = w / n if n else 0.0
    if n < MIN_EVAL:
        why = "few-evaluable"
    elif w < MIN_WET:
        why = "few-flooded"
    elif frac > MAX_WET_FRAC:
        why = "too-flooded"
    else:
        why = "ok"
    return why == "ok", {"n_eval": n, "n_wet": w, "wet_frac": frac, "why": why}
