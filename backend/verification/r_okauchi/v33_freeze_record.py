# %% v33: 事前登録 v2 の凍結記録（§10）の値を成果物から機械的に集める（手で写さない）
"""凍結の直前に走らせ、出力の表を preregistration_v2.md §10 に貼る。値は全てファイルから計算する:

- 凍結物の sha256: N-lex 辞書・EMS 台帳・GeoNames zip（NAS と同一の手元コピー）・海岸線 zip・
  N-read の入力（names.parquet）と出力（nread.parquet、あれば）
- 発生表（incidence.json）の要約: 採用 AOI・activation・場所グループ・国・セル数・Kish 実効数
- N-read: モデル・snapshot・プロンプト sha256（コードから計算）
- 検出力: power_summary.json（帰無・SESOI）の要約

使い方（backend/ から）:
    uv run python verification/r_okauchi/v33_freeze_record.py --power-null power_v2_null \
        --power-sesoi power_v2_sesoi [--nread-run data/global/nread/nread_meta.json]
出力: stdout の Markdown 表と docs/verification/v33_freeze_record.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))

from topsem.world.lexicon import Lexicon
from topsem.world.nread import MODEL, REVISION, prompt_sha256


def sha(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def incidence_summary(inc: dict) -> dict:
    rows = [r for r in inc["aois"] if r["included"]]
    w: Counter = Counter()
    for r in rows:
        w[r["group"]] += r["quota"]
    s = inc["summary"]
    return {"n_aoi_labelled": s["n_aoi_labelled"], "n_aoi_included": s["n_aoi_included"],
            "n_activations": s["n_activations_included"], "n_groups": s["n_location_groups"],
            "n_countries": len({c for r in rows for c in r["countries"]}),
            "group_cap": s["group_cap"], "n_cells": s["n_cells_after_cap"],
            "kish_groups": round(sum(w.values()) ** 2 / sum(v * v for v in w.values()), 1),
            "years": [s["years"][0], s["years"][-1]], "reasons": s["reasons"]}


def power_summary(tag: str) -> dict | None:
    p = ROOT / "data" / "runs" / tag / "power_summary.json"
    if not p.is_file():
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    return {k: d.get(k) for k in ("n_reps", "rate_positive", "rate_equivalent_TOST",
                                  "decision_counts", "mean_R", "median_ci90_halfwidth")}


# %%
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--power-null", default="power_v3_null")
    ap.add_argument("--power-sesoi", default="power_v3_sesoi")
    ap.add_argument("--nread-run", default="data/global/nread/nread_meta.json")
    args = ap.parse_args()
    g = ROOT / "data" / "global"
    lex = Lexicon(ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")
    inc = json.loads((g / "incidence.json").read_text(encoding="utf-8"))
    nmeta_p = ROOT / args.nread_run
    rec = {
        "lexicon": {"version": lex.version, "sha256": lex.sha256},
        "ems_inventory_sha256": sha(ROOT / "docs" / "data" / "ems_inventory_20261025.json"),
        "geonames_allCountries_zip_sha256": sha(g / "geonames" / "allCountries.zip"),
        "coastline_zip_sha256": sha(g / "coast" / "ne_10m_coastline.zip"),
        "incidence": incidence_summary(inc),
        "incidence_json_sha256": sha(g / "incidence.json"),
        "nread": {"model": MODEL, "revision": REVISION, "prompt_sha256": prompt_sha256(),
                  "names_parquet_sha256": sha(g / "nread" / "names.parquet"),
                  "output_parquet_sha256": sha(g / "nread" / "nread.parquet"),
                  "meta": json.loads(nmeta_p.read_text(encoding="utf-8"))
                  if nmeta_p.is_file() else None},
        "power": {"null": {"tag": args.power_null, **(power_summary(args.power_null) or {})},
                  "sesoi": {"tag": args.power_sesoi, **(power_summary(args.power_sesoi) or {})}},
    }
    (ROOT / "docs" / "verification" / "v33_freeze_record.json").write_text(
        json.dumps(rec, indent=1, ensure_ascii=False), encoding="utf-8")
    i = rec["incidence"]
    print("| 項目 | 値 |\n|---|---|")
    print(f"| lexicon_hydro_v1.yaml（版 {lex.version}）sha256 | `{lex.sha256}` |")
    print(f"| ems_inventory_20261025.json sha256 | `{rec['ems_inventory_sha256']}` |")
    print(f"| GeoNames allCountries.zip sha256 | `{rec['geonames_allCountries_zip_sha256']}` |")
    print(f"| Natural Earth ne_10m_coastline.zip sha256 | `{rec['coastline_zip_sha256']}` |")
    print(f"| 発生表 | 採用 {i['n_aoi_included']}/{i['n_aoi_labelled']} AOI・{i['n_activations']} "
          f"activation・{i['n_groups']} 場所グループ・{i['n_countries']} か国・"
          f"{i['n_cells']:,} セル（グループ上限 {i['group_cap']:,}、Kish 実効数 "
          f"{i['kish_groups']}）・{i['years'][0]}–{i['years'][1]} |")
    n = rec["nread"]
    print(f"| N-read | {n['model']} snapshot `{n['revision'][:7]}`・プロンプト sha256 "
          f"`{n['prompt_sha256'][:16]}…`・入力 names.parquet `{(n['names_parquet_sha256'] or '')[:16]}…`"
          f"・出力 `{(n['output_parquet_sha256'] or '（推論待ち）')[:16]}` |")
    for k, v in rec["power"].items():
        print(f"| 検出力（{k}） | run `{v['tag']}`: {json.dumps({x: v.get(x) for x in v if x != 'tag'}, ensure_ascii=False)} |")


if __name__ == "__main__":
    main()
