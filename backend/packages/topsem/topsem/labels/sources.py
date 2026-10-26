"""観測浸水ラベルのソース別ローダー。

各関数は (geometry, attrs, meta) を返す。meta には出所 URL・取得 SHA・
版日・利用条件・既知の弱点を書き、`data/labels/<event>/derived/<river>/
metadata.json` にそのまま保存する（CLAUDE.md「データの来歴」）。
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union


@dataclass(frozen=True)
class FloodExtents:
    """1 系統の浸水範囲（輪郭多角形の和集合）。"""

    river: str
    union_geom: BaseGeometry
    n_polygons: int
    dn_hist: dict[str, int]
    meta: dict[str, Any]


def load_gsi_shinsui_dir(
    extracted_dir: Path,
    version_date: str = "20191018",
    source_url: str = "https://www1.gsi.go.jp/geowww/201910/shinsui/shinsui_rinkaku.zip",
    retrieved_utc: str = "",
) -> list[FloodExtents]:
    """国土地理院「浸水推定段彩図 輪郭線」展開済みディレクトリを読む。

    弱点（免責事項.txt 由来、protocol §3 の注記と対応）:
    - 情報+標高・航空写真+標高による**推定**混じり（純観測ではない）
    - 見落とし・見かけの浸水を地理院自身が宣言
    → 主ラベルに使わず robustness-check 枠に限定する。
    """
    out: list[FloodExtents] = []
    for gj_path in sorted(extracted_dir.rglob("*.geojson")):
        data = json.loads(gj_path.read_text(encoding="utf-8"))
        feats = data.get("features", [])
        polys = [shape(f["geometry"]) for f in feats if f.get("geometry")]
        if not polys:
            continue
        # 2020 系ファイルは crs に EPSG:6668（JGD2011 測地系）を宣言することがある。
        # 数値は lon/lat 順のまま WGS84 と測地異常の差のみ（<0.1 m、250 m 格子で無視可）。
        crs_name = (data.get("crs") or {}).get("properties", {}).get("name") or "CRS84"
        crs_note = (
            "CRS84 (lon/lat)" if "CRS84" in crs_name
            else f"{crs_name} (lon/lat 値・WGS84 と格子精度上同一)"
        )
        base = gj_path.parent.name
        head, _, tail = base.rpartition("_")
        river = head if (tail and tail.isdigit()) else base
        # ファイル名に YYYYMMDD を持つ 2020 系はそれを版日として優先（zip 内で異なる日付混在対策）
        stamp = re.search(r"(?<!\d)(20\d{6})(?!\d)", gj_path.stem)
        file_version = stamp.group(1) if stamp else version_date
        dn_hist: dict[str, int] = {}
        for f in feats:
            key = str((f.get("properties") or {}).get("DN"))
            dn_hist[key] = dn_hist.get(key, 0) + 1
        out.append(
            FloodExtents(
                river=river,
                union_geom=unary_union(polys),
                n_polygons=len(polys),
                dn_hist=dn_hist,
                meta={
                    "source": "gsi_shinsui_rinkaku",
                    "source_url": source_url,
                    "version_date": file_version,
                    "retrieved_utc": retrieved_utc,
                    "crs": crs_note,
                    "file": str(gj_path.relative_to(extracted_dir)),
                    "caveat": "情報/航空写真+標高による推定混じり（国土地理院免責）",
                    "license": "国土地理院コンテンツ利用規約 2.0",
                    "role": "robustness_check_only",
                },
            )
        )
    return out
