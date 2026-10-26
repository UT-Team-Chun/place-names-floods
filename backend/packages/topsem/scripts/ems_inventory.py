"""Copernicus EMS Rapid Mapping の洪水 activation と delineation 製品の台帳を作る。

メタデータのみを取得する（製品 zip は落とさない）。出力:
``data/global/ems/inventory.json`` — activation ごとに国・発生時刻・名称と、
活動ページから抽出した DEL（delineation）系ベクトル zip の URL 一覧。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/ems_inventory.py
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.config import DATA_DIR

API = "https://mapping.emergency.copernicus.eu/activations/api/activations/"
PAGE = "https://mapping.emergency.copernicus.eu/activations/{code}/"
# 旧系の命名は 2 世代ある:
#   2019– : EMSR517_AOI04_DEL_PRODUCT_r1_RTP01_v1 / ..._DEL_MONIT01_...
#   –2018 : EMSR271_01VALTINO_01DELINEATION_MAP_v1 / EMSR167_01LIER_DELINEATION_OVERVIEW_v2 /
#           EMSR107_01PTUJ_DELINEATION_DETAIL01_v2 / ..._OVERVIEW-MONIT01_v1
# REFERENCE（発災前の参照図）・GRADING（被害判読）は対象外。
ZIP_RE = re.compile(
    r"https://cems-mapping-website\.s3\.eu-west-1\.amazonaws\.com/static/activations/"
    r"(EMSR\d+)/(EMSR\d+_[A-Za-z0-9_-]*?(?:_DEL_|DELINEATION)[A-Za-z0-9_-]*?_v\d+)_vector\.zip"
)
NEW_API = "https://rapidmapping.emergency.copernicus.eu/backend/dashboard-api/public-activations/?code={code}"
UA = {"User-Agent": "toponym-semantics research inventory (ryoutaokauchi@g.ecc.u-tokyo.ac.jp)"}
# 新系（2022 年頃以降）の dashboard-api は Origin/Referer つき JSON 要求のみ匿名応答する
# （rm-viewer の fetch と同じ条件）。旧系 activation は 401 → 活動ページの S3 zip に fallback
WEB = {
    "Accept": "application/json",
    "Origin": "https://mapping.emergency.copernicus.eu",
}


def _get(url: str, retries: int = 4, headers: dict | None = None) -> bytes:
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers={**UA, **(headers or {})})  # noqa: S310
            with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310
                return r.read()
        except Exception:
            if k == retries - 1:
                raise
            time.sleep(2 * (k + 1))
    return b""


def list_activations() -> list[dict]:
    out, url = [], f"{API}?limit=200"
    while url:
        d = json.loads(_get(url))
        out.extend(d["results"])
        url = d.get("next")
        if url:
            url = url.replace("http://", "https://")
    return out


def products_new_api(code: str) -> tuple[list[dict], dict] | None:
    """新系 activation の DEL 製品と subCategory/eventTime。旧系なら None。"""
    try:
        raw = _get(NEW_API.format(code=code), retries=2,
                   headers={**WEB, "Referer": PAGE.format(code=code)})
    except Exception:  # noqa: BLE001
        return None
    d = json.loads(raw)
    if not d.get("results"):
        return None
    r = d["results"][0]
    prods = []
    for a in r.get("aois") or []:
        for p in a.get("products") or []:
            if p.get("type") != "DEL" or not p.get("downloadPath"):
                continue
            v = p.get("version") or {}
            prods.append({
                "aoi": f"{int(a.get('number', 0)):02d}",
                "kind": f"MONIT{p['monitoringNumber']:02d}" if p.get("monitoring") else "PRODUCT",
                "version": int(v.get("number") or 0),
                "status": v.get("statusCode", ""),
                "url": p["downloadPath"],
                "image_times": [i.get("acquisitionTime") for i in p.get("images") or []],
            })
    meta = {"sub_category": r.get("subCategory", ""), "event_time": r.get("eventTime", ""),
            "continent": r.get("continent", ""), "system": "new"}
    return prods, meta


def products_for(code: str) -> list[dict]:
    html = _get(PAGE.format(code=code)).decode("utf-8", "replace")
    seen: dict[str, dict] = {}
    for m in ZIP_RE.finditer(html):
        stem = m.group(2)
        a_new = re.search(r"_AOI(\d+)_", stem)
        a_old = re.match(r"EMSR\d+_(\d{2})", stem)
        aoi = a_new.group(1) if a_new else (a_old.group(1) if a_old else "00")
        kind = "MONIT" if "MONIT" in stem else "PRODUCT"
        # 旧世代は同一 AOI に OVERVIEW（広域）と DETAILnn（拡大）がある → scale を記録し
        # ラベル化で AOI ごと 1 製品に絞る（二重計上の防止）
        scale = next((t for t in ("OVERVIEW", "DETAIL", "MAP") if t in stem), "PRODUCT")
        ver = int(re.search(r"_v(\d+)$", stem).group(1))
        seen[stem] = {
            "aoi": aoi,
            "kind": kind,
            "scale": scale,
            "version": ver,
            "url": m.group(0),
        }
    return sorted(seen.values(), key=lambda p: (p["aoi"], p["kind"], p["version"]))


def main() -> None:
    dst = DATA_DIR / "global" / "ems" / "inventory.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    acts = [a for a in list_activations() if (a.get("category") or {}).get("slug") == "flood"]
    print(f"flood activations: {len(acts)}")
    inv = []
    for i, a in enumerate(acts):
        meta = {"system": "old", "sub_category": "", "event_time": ""}
        err = ""
        try:
            got = products_new_api(a["code"])
            if got is not None:
                prods, meta = got
            else:
                prods = products_for(a["code"])
        except Exception as e:  # noqa: BLE001
            prods, err = [], str(e)
        inv.append(
            {
                "code": a["code"],
                "name": a.get("name", ""),
                "countries": [c.get("short_name", "") for c in a.get("countries") or []],
                "activation_time": a.get("activationTime", ""),
                "centroid": a.get("centroid", ""),
                "n_aois": a.get("n_aois"),
                **meta,
                "del_products": prods,
                "error": err,
            }
        )
        if i % 25 == 0:
            print(f"  {i:4d}/{len(acts)} {a['code']} del={len(prods)}", flush=True)
        time.sleep(0.4)
    payload = {
        "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "api": API,
        "n_flood_activations": len(inv),
        "activations": inv,
    }
    dst.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    n_del = sum(1 for a in inv for p in a["del_products"] if p["kind"] == "PRODUCT")
    print(f"→ {dst}  DEL_PRODUCT zips: {n_del}")


if __name__ == "__main__":
    main()
