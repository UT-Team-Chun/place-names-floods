"""年次土地被覆（ESA CCI / C3S、300 m）を洪水の**前年**版で AOI の副格子へ（時点管理の準拠）。

WorldCover は 2021 年版しかなく、2012–2020 年の洪水には事後の情報になる（洪水で生じた裸地や
流路の変化がラベル情報を基準モデルへ漏らし得る）。ESA CCI Land Cover（1992–2015, v2.0.7cds）と
C3S Land Cover（2016–2020, v2.1.1）は毎年の版があるので、各 AOI に洪水の前年の版を当てる
（2021 年以降の洪水は最新の 2020 年版）。どの事例でも洪水より前の版になる。

前年の版なので、洪水年の撮像が混入する懸念（WorldCover で水域・湿地を外した理由）は無い。
地図化された湿地・水域も基準モデルに入れ、名前はそれを超える情報を持つかで試される。

出所: Microsoft Planetary Computer の COG（collection ``esa-cci-lc``、45° タイル）。
資産 URL は匿名 SAS トークンで署名して GDAL で読む。
"""

from __future__ import annotations

import calendar
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
from rasterio.warp import reproject

from topsem.world.terrain import _block, grid_lonlat_bounds

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS_ENDPOINT = "https://planetarycomputer.microsoft.com/api/sas/v1/token/esa-cci-lc"
COLLECTION = "esa-cci-lc"
FIRST_YEAR, LAST_YEAR = 1992, 2020
# LCCS の区分 → 8 群（0 = 欠測）。40（自然植生 > 50% のモザイク）は灌木に、110（草本 > 50% の
# モザイク）は草地に、まばらな植生（150–153）と雪氷（220）は裸地に寄せる
CCI_GROUPS: dict[str, tuple[int, ...]] = {
    "crop": (10, 11, 12, 20, 30),
    "tree": (50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100),
    "shrub": (40, 120, 121, 122),
    "grass": (110, 130, 140),
    "bare": (150, 151, 152, 153, 200, 201, 202, 220),
    "wetland": (160, 170, 180),
    "built": (190,),
    "water": (210,),
}
GDAL_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "GDAL_HTTP_MULTIRANGE": "YES",
    "GDAL_HTTP_MERGE_CONSECUTIVE_RANGES": "YES",
    "VSI_CACHE": "TRUE",
    "GDAL_HTTP_MAX_RETRY": "5",
    "GDAL_HTTP_RETRY_DELAY": "2",
}
_TOKEN: dict[str, Any] = {"value": "", "expiry": 0.0}
RETRIES = 8


def lc_year(event_time: str) -> int:
    """洪水の前年（最新の 2020 年版で頭打ち）。"""
    return max(FIRST_YEAR, min(int(event_time[:4]) - 1, LAST_YEAR))


def _get_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "toponym-semantics research"})  # noqa: S310 — 固定の https
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310
                return json.loads(r.read().decode("utf-8"))
        except OSError as e:
            if attempt == RETRIES - 1:
                raise
            # 回数制限（429）は Retry-After に従い、無ければ指数的に待つ。
            # 並列ワーカーの同時再試行を散らすため小さな揺らぎを足す
            hdrs = getattr(e, "headers", None)
            after = hdrs.get("Retry-After") if hdrs else None
            backoff = min(120.0, 4.0 * 2**attempt)
            wait = float(after) if after and after.isdigit() else backoff
            time.sleep(wait + random.uniform(0, 3))  # noqa: S311 — 待ち時間の揺らぎのみ
    msg = "unreachable"
    raise RuntimeError(msg)


def _token() -> str:
    """匿名 SAS トークン（期限の 5 分前に取り直す）。"""
    if time.time() > _TOKEN["expiry"] - 300:
        d = _get_json(SAS_ENDPOINT)
        _TOKEN["value"] = d["token"]
        _TOKEN["expiry"] = calendar.timegm(time.strptime(d["msft:expiry"][:19],
                                                         "%Y-%m-%dT%H:%M:%S"))
    return str(_TOKEN["value"])


def item_hrefs(bounds: tuple[float, float, float, float], year: int) -> list[str]:
    """AOI の lon/lat 範囲に掛かるタイルの lccs_class COG（署名前の URL）。"""
    q = urllib.parse.urlencode({"collections": COLLECTION, "bbox": ",".join(map(str, bounds)),
                                "datetime": f"{year}-01-01/{year}-12-31", "limit": 20})
    feats = _get_json(f"{STAC}?{q}").get("features", [])
    return sorted(f["assets"]["lccs_class"]["href"] for f in feats)


def group_lut() -> np.ndarray:
    """LCCS 区分（uint8）→ 群番号（1..8、0 = 欠測・未知）。"""
    lut = np.zeros(256, dtype=np.uint8)
    for gi, codes in enumerate(CCI_GROUPS.values(), start=1):
        lut[list(codes)] = gi
    return lut


def cci_features(meta: dict, event_time: str) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """1 AOI の全セル（格子全体 ny*nx、行 0 = 北）の CCI 土地被覆群の被覆率と来歴。"""
    g = meta["grid"]
    ny, nx = g["ny"], g["nx"]
    sub = round(g["cell_m"] / g["fine_res_m"])
    shape = (ny * sub, nx * sub)
    year = lc_year(event_time)
    hrefs = item_hrefs(grid_lonlat_bounds(meta), year)
    cls = np.zeros(shape, dtype=np.uint8)
    dst_tr = rasterio.transform.from_origin(g["x0"], g["y1"], g["fine_res_m"], g["fine_res_m"])
    dst_crs = CRS.from_proj4(meta["crs_proj4"])
    with rasterio.Env(**GDAL_ENV):
        for href in hrefs:
            with rasterio.open(f"/vsicurl/{href}?{_token()}") as src:
                tmp = np.zeros(shape, dtype=np.uint8)
                reproject(source=rasterio.band(src, 1), destination=tmp, dst_transform=dst_tr,
                          dst_crs=dst_crs, src_nodata=0, dst_nodata=0,
                          resampling=Resampling.nearest)
                cls[tmp != 0] = tmp[tmp != 0]
    grp = group_lut()[cls]
    gb = _block(grp, sub, ny, nx)
    out = {f"cci_{k}": np.mean(gb == gi, axis=2).astype(np.float32)
           for gi, k in enumerate(CCI_GROUPS, start=1)}
    out["cci_valid_frac"] = np.mean(gb != 0, axis=2).astype(np.float32)
    prov = {"collection": COLLECTION, "year": year, "event_time": event_time,
            "hrefs": hrefs, "groups": {k: list(v) for k, v in CCI_GROUPS.items()}}
    return {k: v.ravel() for k, v in out.items()}, prov
