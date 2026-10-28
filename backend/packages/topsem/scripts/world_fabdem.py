"""amendment B-9 の C7: FABDEM v1-2 の地形の特徴（建物・樹冠を除いた標高）と河道からの高さを AOI ごとに作る。

    python world_fabdem.py --prepare            # 解析する AOI に要る 1° タイルを NAS へ取り出す
    python world_fabdem.py --aois 0:50          # 採用 AOI（名前の順）の範囲の特徴を作る
    python world_fabdem.py --summary            # 解析するセルのうち FABDEM の標高が得られた割合（実行の条件）

出力: 各 AOI のラベルの場所に ``features_fabdem.parquet``（cells.parquet と同じ行の順。elev_1km・elev_mean・depress・
elev_min・elev_std・slope_mean・dem_valid_frac・h_above_100・h_above_1000・dist_100_m・dist_1000_m）と
``features_fabdem_meta.json``。タイルは ``data/global/fabdem/tiles/``（海域で zip に無いタイルは ``*.absent``）。

配布元は 1 接続あたり約 0.5 MB/s に絞っている（2026-10-28 に手元とクラスタで測った）。10° の zip を丸ごと
取ると 62 個で数十時間かかるので、zip の中の要るタイルだけを範囲指定の要求で読み、タイルごとに並列に取る。
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
import traceback
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))

import world_eval as ev

from topsem.world.b9 import (
    FABDEM_BASE,
    RIVER_PAD_DEG,
    fabdem_elevation_features,
    fabdem_tile_name,
    fabdem_zip_name,
    river_heights,
    tiles_for,
    warp_fabdem,
)
from topsem.world.terrain import grid_lonlat_bounds

ROOT = ev.ROOT
CACHE = ROOT / "fabdem"
TILES = CACHE / "tiles"
WORKERS = 32  # 並列の接続の数（1 接続 ≈ 0.5 MB/s）
FETCH_MIN = 1 << 16  # 範囲要求の最小の大きさ（zip の頭書きを読むときの往復を減らす）
CHUNK = 1 << 23  # タイルを読む 1 回の範囲要求の大きさ
RETRIES = 6
OUT_NAME = "features_fabdem.parquet"
COVER_RUN = 0.95  # 実行の条件: 解析するセルのうち FABDEM の標高が得られた割合（B-9 の C7）
COVER_LOW = 0.9  # AOI ごとの記録で「低い」とする割合


def _gn() -> pd.DataFrame:
    return pd.read_parquet(ROOT / "geonames" / "aoi_subset.parquet")  # world_eval と同じ


def included_dirs() -> list[Path]:
    dirs, _, _ = ev.select_aois(_gn())
    return dirs


def _padded(b: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    w, s, e, n = b
    return w - RIVER_PAD_DEG, s - RIVER_PAD_DEG, e + RIVER_PAD_DEG, n + RIVER_PAD_DEG


def needed_tiles(dirs: list[Path]) -> list[tuple[int, int]]:
    """AOI の範囲 ± RIVER_PAD_DEG にかかる 1° タイル（河道の頂点の標高も読むので広げた範囲）。"""
    tiles = set()
    for d in dirs:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        tiles |= set(tiles_for(_padded(grid_lonlat_bounds(meta))))
    return sorted(tiles)


class _HttpRange:
    """範囲指定の要求で読む、読み出し専用・シーク可能なファイル（zipfile に渡す）。"""

    def __init__(self, url: str) -> None:
        self.url, self.pos = url, 0
        self.buf, self.buf_start = b"", 0
        req = urllib.request.Request(url, method="HEAD", headers=_HEADERS)  # noqa: S310
        with _retry(lambda: urllib.request.urlopen(req, timeout=120)) as resp:  # noqa: S310
            self.size = int(resp.headers["Content-Length"])

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, off: int, whence: int = 0) -> int:
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def read(self, n: int = -1) -> bytes:
        n = self.size - self.pos if n is None or n < 0 else min(n, self.size - self.pos)
        if n <= 0:
            return b""
        off = self.pos - self.buf_start
        if not (off >= 0 and off + n <= len(self.buf)):
            end = min(self.pos + max(n, FETCH_MIN), self.size) - 1
            self.buf, self.buf_start, off = self._get(self.pos, end), self.pos, 0
        out = self.buf[off:off + n]
        self.pos += len(out)
        return out

    def _get(self, a: int, b: int) -> bytes:
        def once() -> bytes:
            req = urllib.request.Request(self.url, headers=_HEADERS | {"Range": f"bytes={a}-{b}"})  # noqa: S310
            with urllib.request.urlopen(req, timeout=300) as resp:  # noqa: S310
                data = resp.read()
            if len(data) != b - a + 1:
                msg = f"short read {len(data)} of {b - a + 1}"
                raise OSError(msg)
            return data
        return _retry(once)

    def close(self) -> None:
        self.buf = b""


_HEADERS = {"User-Agent": "toponym-semantics research"}


def _retry(fn):  # noqa: ANN001, ANN202
    for k in range(RETRIES):
        try:
            return fn()
        except urllib.error.HTTPError as e:
            if e.code == 404 or k == RETRIES - 1:  # noqa: PLR2004
                raise
        except (OSError, TimeoutError):
            if k == RETRIES - 1:
                raise
        time.sleep(2 ** k)
    msg = "unreachable"
    raise RuntimeError(msg)


_ZIP_MEMBERS: dict[str, set[str] | None] = {}
_ZIP_LOCK = threading.Lock()


def _members(zname: str) -> set[str] | None:
    """zip の中のタイルの名前（zip が無ければ None、海域だけの区画）。zip ごとに一度だけ読む。"""
    with _ZIP_LOCK:
        if zname not in _ZIP_MEMBERS:
            try:
                with zipfile.ZipFile(_HttpRange(f"{FABDEM_BASE}/{zname}")) as zf:
                    _ZIP_MEMBERS[zname] = {Path(x).name for x in zf.namelist()}
            except urllib.error.HTTPError as e:
                if e.code != 404:  # noqa: PLR2004
                    raise
                _ZIP_MEMBERS[zname] = None
        return _ZIP_MEMBERS[zname]


def fetch_tile(lat: int, lon: int) -> tuple[str, int]:
    """1° タイルを zip から取り出して TILES に置く。戻り値: 状態と大きさ（バイト）。"""
    z, t = fabdem_zip_name(lat, lon), fabdem_tile_name(lat, lon)
    dst, absent = TILES / t, TILES / f"{t}.absent"
    if dst.is_file():
        return "cached", dst.stat().st_size
    if absent.is_file():
        return "absent", 0
    names = _members(z)
    if names is None or t not in names:
        absent.write_text(z, encoding="utf-8")
        return "absent", 0
    tmp = TILES / f"{t}.part"
    with zipfile.ZipFile(_HttpRange(f"{FABDEM_BASE}/{z}")) as zf:
        member = next(x for x in zf.namelist() if Path(x).name == t)
        with zf.open(member) as src, tmp.open("wb") as f:  # 最後まで読むと CRC を確かめる
            while chunk := src.read(CHUNK):
                f.write(chunk)
        if tmp.stat().st_size != zf.getinfo(member).file_size:
            msg = f"{t}: size mismatch"
            raise OSError(msg)
    tmp.rename(dst)
    return "downloaded", dst.stat().st_size


def prepare(workers: int = WORKERS, limit: int | None = None) -> None:
    TILES.mkdir(parents=True, exist_ok=True)
    tiles = needed_tiles(included_dirs())[:limit]
    print(f"{len(tiles)} tiles needed ({len({fabdem_zip_name(*x) for x in tiles})} zips)", flush=True)
    status: dict[str, str] = {}
    t0, done, nbytes = time.time(), 0, 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_tile, la, lo): fabdem_tile_name(la, lo) for la, lo in tiles}
        for fu in as_completed(futs):
            name = futs[fu]
            try:
                st, size = fu.result()
            except Exception as e:  # noqa: BLE001
                st, size = f"error: {type(e).__name__}: {e}", 0
            status[name], done, nbytes = st, done + 1, nbytes + size
            el = time.time() - t0
            print(f"  [{done}/{len(tiles)}] {name}: {st} {size / 1e6:.1f} MB "
                  f"(total {nbytes / 1e9:.2f} GB, {nbytes / 1e6 / max(el, 1):.1f} MB/s, {el:.0f}s)", flush=True)
    counts = {k: sum(v.startswith(k) for v in status.values())
              for k in ("downloaded", "cached", "absent", "error")}
    print(json.dumps(counts), flush=True)
    (CACHE / "prepare.json").write_text(json.dumps({"counts": counts, "tiles": status}, indent=1),
                                       encoding="utf-8")


def river_vertices(bounds: tuple[float, float, float, float]) -> pd.DataFrame:
    """AOI の範囲 ± RIVER_PAD_DEG にかかる HydroRIVERS の河道の頂点と集水面積（B-8 の B7 と同じ読み方）。"""
    import pyogrio.raw
    import shapely

    w, s, e, n = _padded(bounds)
    parts = []
    for z in sorted((ROOT / "hydrorivers").glob("HydroRIVERS_v10_*_shp.zip")):
        with zipfile.ZipFile(z) as zf:
            shp = next(x for x in zf.namelist() if x.endswith(".shp"))
        _, _, wkb, fields = pyogrio.raw.read(f"/vsizip/{z}/{shp}", columns=["UPLAND_SKM"], bbox=(w, s, e, n))
        if len(wkb) == 0:
            continue
        geoms = shapely.from_wkb(wkb)
        xy, ix = shapely.get_coordinates(geoms, return_index=True)
        up = np.asarray(fields[0], dtype=np.float64)[ix]
        parts.append(np.column_stack([xy[:, 0], xy[:, 1], up]))
    if not parts:
        return pd.DataFrame(columns=["lon", "lat", "upland"])
    return pd.DataFrame(np.concatenate(parts), columns=["lon", "lat", "upland"])


def one(d: Path, *, force: bool = False) -> dict:
    out = d / OUT_NAME
    if not force and out.is_file():
        return {"aoi": d.name, "status": "cached"}
    t0 = time.time()
    try:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        cells = pq.read_table(d / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        g = meta["grid"]
        flat = (cells["iy"].to_numpy() * g["nx"] + cells["ix"].to_numpy()).astype(np.int64)
        dem, n_ok, n_miss = warp_fabdem(meta, CACHE)
        feat = fabdem_elevation_features(dem, meta)
        vert = river_vertices(grid_lonlat_bounds(meta))
        feat |= river_heights(feat["elev_mean"], meta, vert, CACHE)
        cols = {k: np.asarray(v)[flat] for k, v in feat.items()}
        pq.write_table(pa.table(cols), out, compression="zstd")
        info = {"aoi": d.name, "status": "ok", "tiles_read": n_ok, "tiles_missing": n_miss,
                "n_cells": len(cells), "valid_share": float(np.mean(np.isfinite(cols["elev_mean"]))),
                "h100_share": float(np.mean(np.isfinite(cols["h_above_100"]))),
                "n_river_vertices": len(vert), "sec": round(time.time() - t0, 1),
                "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        (d / "features_fabdem_meta.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    except Exception as e:  # noqa: BLE001
        return {"aoi": d.name, "status": "error", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-1500:]}
    else:
        return info


def summary() -> dict:
    """実行の条件（B-9 の C7）: 解析するセル（確証と同じ上限つきの表）のうち FABDEM の標高が得られた割合。"""
    dirs, _, infos = ev.select_aois(_gn())
    df = ev.build_table(dirs, seed=20261025, infos=infos, extra_cols=("ix", "iy"))
    ok, n = 0, 0
    per = {}
    for aoi, x in df.groupby("aoi"):
        f = ROOT / "labels" / aoi / OUT_NAME
        if not f.is_file():
            per[aoi] = None
            n += len(x)
            continue
        cells = pq.read_table(ROOT / "labels" / aoi / "cells.parquet", columns=["ix", "iy"]).to_pandas()
        fe = pq.read_table(f, columns=["elev_mean", "h_above_100"]).to_pandas()
        fe["ix"], fe["iy"] = cells["ix"].to_numpy(), cells["iy"].to_numpy()
        m = x[["ix", "iy"]].merge(fe, on=["ix", "iy"], how="left")
        v = np.isfinite(m["elev_mean"].to_numpy(np.float64))
        per[aoi] = float(v.mean())
        ok += int(v.sum())
        n += len(x)
    res = {"rule": "amendment B-9 item C7: run only if FABDEM covers at least 95 % of the analysed cells",
           "n_cells": n, "n_valid": ok, "share_valid": ok / max(n, 1),
           "condition_met": ok / max(n, 1) >= COVER_RUN,
           "aois_missing": sorted(a for a, v in per.items() if v is None),
           "aois_below_90pct": sorted(a for a, v in per.items() if v is not None and v < COVER_LOW)}
    (ROOT / "fabdem" / "coverage.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--limit", type=int, default=None, help="--prepare で取るタイルの数の上限（試験用）")
    ap.add_argument("--aois", default="", help="採用 AOI（名前の順）の範囲 a:b")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if args.prepare:
        prepare(args.workers, args.limit)
    if args.aois:
        a, b = (int(x) for x in args.aois.split(":"))
        dirs = included_dirs()[a:b]
        for d in dirs:
            r = one(d, force=args.force)
            print(json.dumps({k: v for k, v in r.items() if k != "trace"}), flush=True)
            if r["status"] == "error":
                print(r["trace"], flush=True)
    if args.summary:
        print(json.dumps(summary(), indent=1), flush=True)


if __name__ == "__main__":
    main()
