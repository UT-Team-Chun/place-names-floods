# %% v49（事後・amendment B-6 の 3(b)(c)）: 集落ごとの繰り返す水（JRC Global Surface Water v1.4）
"""解析した AOI の集落（GeoNames P 類）ごとに、点から 500 m 以内の 30 m 画素の JRC GSW occurrence（水があった
時間の割合、%）と seasonality（水のあった月の数）の平均を出す。恒常水域（occurrence ≥ 90 かつ seasonality = 12）
を除いた平均も。浸水のラベルは使わない。GSW のタイル（10° 四方、EPSG:4326、0.00025°）は公開の GeoTIFF を
範囲読み出しで AOI ごとに読む（Pekel et al. 2016, Nature 540; v1.4、1984–2021）。

出力: ``data/global/gsw/settlement_gsw.parquet``（gid, occ, seas, occ_np, seas_np, n_px, n_px_np）と
``settlement_gsw_meta.json``（タイル・取得日時・規則）。クラスタの script 段で実行する（ネットワークに出る）。
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())  # コンテナには CLAUDE.md が無い
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import rasterio
import world_eval as ev
from rasterio.windows import from_bounds

from topsem.config import DATA_DIR

URL = "https://storage.googleapis.com/global-surface-water/downloads2021/{layer}/{layer}_{lon}_{lat}v1_4_2021.tif"
OUT = DATA_DIR / "global" / "gsw"
RADIUS_M = 500.0
PAD_DEG = 0.06  # 集落の点の取り込み（AOI の範囲 ± 0.05°）と 500 m の余白
NODATA = 255
PERM_OCC, PERM_SEAS = 90, 12
ENV = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR", "GDAL_HTTP_MAX_RETRY": "8", "GDAL_HTTP_RETRY_DELAY": "2",
       "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif"}


def tile_name(lon0: int, lat1: int) -> tuple[str, str]:
    """タイルの左上の角（経度は 10 の倍数、緯度は 10 の倍数）→ 名前の部品（例 10E, 50N）。"""
    return (f"{abs(lon0)}{'E' if lon0 >= 0 else 'W'}", f"{abs(lat1)}{'N' if lat1 >= 0 else 'S'}")


def read_layer(layer: str, w: float, s: float, e: float, n: float) -> tuple[np.ndarray, tuple[float, float]]:
    """範囲 [w, e] × [s, n] の GSW の層を、0.00025° の格子で（タイルをまたげば継ぎ合わせて）読む。"""
    res = 0.00025
    nx = math.ceil((e - w) / res)
    ny = math.ceil((n - s) / res)
    out = np.full((ny, nx), NODATA, dtype=np.uint8)
    for lon0 in range(int(math.floor(w / 10) * 10), int(math.floor(e / 10) * 10) + 1, 10):
        for lat1 in range(int(math.ceil(s / 10) * 10), int(math.ceil(n / 10) * 10) + 1, 10):
            lo, la = tile_name(lon0, lat1)
            url = URL.format(layer=layer, lon=lo, lat=la)
            tw, te, ts, tn = max(w, lon0), min(e, lon0 + 10), max(s, lat1 - 10), min(n, lat1)
            if tw >= te or ts >= tn:
                continue
            try:
                with rasterio.Env(**ENV), rasterio.open(f"/vsicurl/{url}") as src:
                    win = from_bounds(tw, ts, te, tn, src.transform).round_offsets().round_lengths()
                    arr = src.read(1, window=win, boundless=True, fill_value=NODATA)
            except rasterio.errors.RasterioIOError:  # 陸の無いタイル（海）は公開されていない
                continue
            c0 = round((tw - w) / res)
            r0 = round((n - tn) / res)
            h, wd = min(arr.shape[0], ny - r0), min(arr.shape[1], nx - c0)
            out[r0:r0 + h, c0:c0 + wd] = arr[:h, :wd]
    return out, (w, n)


def settlement_means(occ: np.ndarray, seas: np.ndarray, origin: tuple[float, float],
                     pts: pd.DataFrame) -> pd.DataFrame:
    res = 0.00025
    w, n = origin
    rows = []
    for gid, lon, lat in zip(pts["gid"].to_numpy(), pts["lon"].to_numpy(float), pts["lat"].to_numpy(float),
                             strict=True):
        dlat = RADIUS_M / 110_540.0
        dlon = RADIUS_M / (111_320.0 * max(math.cos(math.radians(lat)), 1e-3))
        c0, c1 = int((lon - dlon - w) / res), int((lon + dlon - w) / res) + 1
        r0, r1 = int((n - (lat + dlat)) / res), int((n - (lat - dlat)) / res) + 1
        c0, r0 = max(c0, 0), max(r0, 0)
        o = occ[r0:r1, c0:c1].astype(np.float64)
        sz = seas[r0:r1, c0:c1].astype(np.float64)
        if o.size == 0:
            rows.append((gid, np.nan, np.nan, np.nan, np.nan, 0, 0))
            continue
        cx = w + (np.arange(c0, c0 + o.shape[1]) + 0.5) * res
        cy = n - (np.arange(r0, r0 + o.shape[0]) + 0.5) * res
        dx = (cx[None, :] - lon) * 111_320.0 * math.cos(math.radians(lat))
        dy = (cy[:, None] - lat) * 110_540.0
        inside = (dx ** 2 + dy ** 2) <= RADIUS_M ** 2
        ok = inside & (o != NODATA) & (sz != NODATA)
        perm = (o >= PERM_OCC) & (sz == PERM_SEAS)
        ok_np = ok & ~perm
        rows.append((gid, float(o[ok].mean()) if ok.any() else np.nan,
                     float(sz[ok].mean()) if ok.any() else np.nan,
                     float(o[ok_np].mean()) if ok_np.any() else np.nan,
                     float(sz[ok_np].mean()) if ok_np.any() else np.nan, int(ok.sum()), int(ok_np.sum())))
    return pd.DataFrame(rows, columns=["gid", "occ", "seas", "occ_np", "seas_np", "n_px", "n_px_np"])


def main() -> None:
    t0 = time.time()
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    dirs, _, infos = ev.select_aois(gn, ev.Sensitivity())
    p = gn.loc[gn["fclass"] == "P"]
    parts = []
    for k, d in enumerate(dirs):
        w, s, e, n = infos[d.name].bounds
        pts = p.loc[p["lon"].between(w - 0.05, e + 0.05) & p["lat"].between(s - 0.05, n + 0.05)]
        if len(pts) == 0:
            continue
        box = (w - PAD_DEG, s - PAD_DEG, e + PAD_DEG, n + PAD_DEG)
        occ, origin = read_layer("occurrence", *box)
        seas, _ = read_layer("seasonality", *box)
        parts.append(settlement_means(occ, seas, origin, pts))
        if k % 25 == 0:
            print(f"{k + 1}/{len(dirs)} {d.name}: {len(pts)} settlements ({time.time() - t0:.0f}s)", flush=True)
    out = pd.concat(parts, ignore_index=True)
    out = out.groupby("gid", as_index=False).first()  # AOI が重なれば同じ集落が二度出る（値は同じ）
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT / "settlement_gsw.parquet", index=False)
    meta = {"rule": "amendment B-6 item 3(b)", "source": "JRC Global Surface Water v1.4 (1984-2021)",
            "url_pattern": URL, "radius_m": RADIUS_M, "permanent": {"occ_ge": PERM_OCC, "seas_eq": PERM_SEAS},
            "n_settlements": len(out), "share_with_pixels": float((out["n_px"] > 0).mean()),
            "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "sec": round(time.time() - t0)}
    (OUT / "settlement_gsw_meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
