# %% v62（事後・amendment B-8 追記 7）: 厳密な読みの精度を言語別に（査読者 3 の細かい指摘 10）
"""確証判定は変えない。最初の符号化者の符号（v44）で、厳密な読み（B-6 の 5: 泉・井戸を除く水・湿地・河川）の精度を、
符号のある名が 30 以上の言語ごとに（Wilson の 95% 区間）。広い読みの言語別の精度と並べる。記述のみ。

出力: stdout と docs/verification/v62_strict_by_language.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

from topsem.world.b6 import STRICT_CATEGORIES, strict_lexicon
from topsem.world.lexicon import Lexicon

OUT = ROOT / "docs" / "verification" / "v62_strict_by_language.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    v51 = _load("v51_readings")
    v44 = v51.v44
    d, _sizes = v51.coded()
    c = d[d["status"] == "coded"].copy()
    ls = strict_lexicon(Lexicon(v51.LEX))
    c["strict_flag"] = [any(ls.match(r.name, r.cc)[k] for k in STRICT_CATEGORIES) for r in c.itertuples()]
    c["spring_only"] = [("water" in r.cats_set) and v51.spring_only(r.quote) for r in c.itertuples()]
    c["strict_truth"] = [bool((r.cats_set - ({"water"} if r.spring_only else set())) & set(STRICT_CATEGORIES))
                         for r in c.itertuples()]
    c["flag"] = c["lex_wetside"].astype(bool)
    out: dict = {"rule": "amendment B-8 addendum 7 (descriptive)", "codes": "first coder (v44)",
                 "min_group": v44.MIN_GROUP, "by_language": {}}
    for lang, x in c.groupby("lang"):
        if len(x) < v44.MIN_GROUP:
            continue
        fl, sf = x["flag"], x["strict_flag"]
        out["by_language"][lang] = {
            "n_coded": len(x),
            "broad": v44._rate(int((fl & x["truth"]).sum()), int(fl.sum())),
            "strict": v44._rate(int((sf & x["strict_truth"]).sum()), int(sf.sum())),
        }
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    for lang, r in out["by_language"].items():
        print(lang, r["n_coded"], "broad", r["broad"], "strict", r["strict"])


if __name__ == "__main__":
    main()
