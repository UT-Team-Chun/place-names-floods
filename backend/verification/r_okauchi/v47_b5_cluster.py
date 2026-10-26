# %% v47（事後・amendment B-3〜B-5）: クラスタで回した事後の実行の集計
"""確証判定は変えない。B-5 の 1（置換の分布）・2（二段目の係数）・3（fold の割り当ての繰り返し）・5（オラクル対照）・
7（ラベルの変種）、B-3 の C2（GFM で学習し直し）、B-4（辞書を直した再実行）の値を、run の成果物から集める。

入力（data/runs/、world_aggregate の集計と held-out 予測。NAS から fetch_nfs_run.sh で引き上げる）:
- world_v2_b5perm/predictions/primary/fold*_perms.npz（置換 200 回の held-out 予測）
- world_v2_b5coef/predictions/primary/fold*.npz・fold*.json（M1co・M3o の予測、二段目の係数）
- world_v2_lexfix・world_v2_fs1b・world_v2_fs2b・world_v2_fs3 … world_v2_fs12 の metrics_primary.json
  （シード 0 は world_v2_lexfix。B-5 追記 3）、比較のため修正前の world_v2・world_v2_fs1・world_v2_fs2
- world_v2_lv_notrace・world_v2_lv_nonwater・world_v2_gfm の metrics_primary.json
- world_v2_lexfix の metrics_{primary,country,temporal,llmcut}.json

出力: stdout と docs/verification/v47_b5_cluster.json（無い入力は "missing" と記録する）
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_aggregate as wa

from topsem.config import RUNS_DIR

OUT = ROOT / "docs" / "verification" / "v47_b5_cluster.json"
SEEDS_FIXED = {0: "world_v2_lexfix", 1: "world_v2_fs1b", 2: "world_v2_fs2b",
               **{s: f"world_v2_fs{s}" for s in range(3, 13)}}
SEEDS_FROZEN = {0: "world_v2", 1: "world_v2_fs1", 2: "world_v2_fs2"}
Z90 = 1.645


def metrics(run: str, split: str = "primary") -> dict | None:
    p = RUNS_DIR / run / f"metrics_{split}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else None


def brief(m: dict | None) -> dict | str:
    if not m:
        return "missing"
    p = m["primary"]
    return {"R": p["R"], "ci95_R": p["ci95_R"], "ci90_R": p["ci90_R"], "delta_bs": p["delta_bs_M1c_M3n"],
            "ci95_delta": p["ci95_delta"], "yardstick": p["delta_bs_M0_M1"],
            "ci95_yardstick": p["ci95_delta_M0_M1"], "decision": p["decision"],
            "n_cells": m.get("n_cells"), "n_units": m.get("n_units")}


# %% 3. fold の割り当ての繰り返し
def fold_seeds() -> dict:
    out: dict = {}
    for name, seeds in (("fixed_code", SEEDS_FIXED), ("frozen_code", SEEDS_FROZEN)):
        rows = {s: brief(metrics(r)) for s, r in seeds.items()}
        got = {s: v for s, v in rows.items() if isinstance(v, dict)}
        blk: dict = {"runs": {s: seeds[s] for s in seeds}, "n_available": len(got),
                     "missing": [seeds[s] for s, v in rows.items() if not isinstance(v, dict)],
                     "R": {s: v["R"] for s, v in got.items()},
                     "decision": {s: v["decision"] for s, v in got.items()}}
        if 0 in got and len(got) >= 2:
            rs = np.array([v["R"] for v in got.values()])
            se = (got[0]["ci90_R"][1] - got[0]["ci90_R"][0]) / (2 * Z90)
            sd = float(np.std(rs, ddof=1))
            half = Z90 * float(np.sqrt(se ** 2 + sd ** 2))
            blk |= {"R_min": float(rs.min()), "R_max": float(rs.max()), "R_mean": float(rs.mean()),
                    "R_sd": sd, "boot_se_seed0": se,
                    "combined_ci90": [got[0]["R"] - half, got[0]["R"] + half],
                    "n_equivalent_null": int(sum(v["decision"] == "equivalent-null" for v in got.values()))}
        out[name] = blk
    return out


# %% 1. 置換の分布
def permutations(yardstick: float | None) -> dict | str:
    d = RUNS_DIR / "world_v2_b5perm" / "predictions" / "primary"
    files = sorted(d.glob("fold*_perms.npz")) if d.is_dir() else []
    if len(files) < 10 or yardstick is None:
        return f"missing ({len(files)} of 10 folds)"
    y, m1c, m3n, perms = [], [], [], []
    for f in files:
        with np.load(f) as z:
            y.append(z["y"].astype(np.float64))
            m1c.append(z["p_M1c"].astype(np.float64))
            m3n.append(z["p_M3n"].astype(np.float64))
            perms.append(z["p_perms"].astype(np.float64))
    y, m1c, m3n = np.concatenate(y), np.concatenate(m1c), np.concatenate(m3n)
    pp = np.concatenate(perms, axis=1)
    bs_c = float(np.mean((m1c - y) ** 2))
    d_real = bs_c - float(np.mean((m3n - y) ** 2))
    d_perm = bs_c - np.mean((pp - y[None, :]) ** 2, axis=1)
    r_perm = d_perm / yardstick
    return {"n_perms": int(pp.shape[0]), "n_cells": len(y), "yardstick": yardstick,
            "delta_bs_real": d_real, "R_real": d_real / yardstick,
            "R_perm_mean": float(r_perm.mean()), "R_perm_q025": float(np.quantile(r_perm, 0.025)),
            "R_perm_q975": float(np.quantile(r_perm, 0.975)),
            "share_perm_ge_real": float(np.mean(d_perm >= d_real)),
            "percentile_of_real": float(100 * np.mean(d_perm < d_real))}


# %% 5. オラクル対照と 2. 二段目の係数
def oracle_and_coefs() -> dict | str:
    run = RUNS_DIR / "world_v2_b5coef"
    d = run / "predictions" / "primary"
    if not d.is_dir() or len(list(d.glob("fold[0-9].json"))) < 10:
        return "missing"
    data = wa.load_split(run, "primary")
    preds = {k[2:]: v for k, v in data.items() if k.startswith("p_")}
    s, n, arms = wa.unit_sums(data["y"].astype(np.float64), preds, data["unit"])
    hat = s.sum(axis=0) / n.sum()
    boot = wa.bootstrap_bs(s, n)
    delta = wa.make_delta(hat, boot, {a: i for i, a in enumerate(arms)})
    d01, b01 = delta("M0", "M1")
    res: dict = {"n_cells": int(n.sum()), "n_units": len(n), "yardstick": d01,
                 "ci95_yardstick": wa.ci(b01, 0.95)}
    for key, base, alt in (("oracle", "M1co", "M3o"), ("names", "M1c", "M3n")):
        dn, bn = delta(base, alt)
        r_b = bn / d01
        ci95, ci90 = wa.ci(r_b, 0.95), wa.ci(r_b, 0.90)
        res[key] = {"comparison": f"{base}->{alt}", "delta_bs": dn, "ci95_delta": wa.ci(bn, 0.95),
                    "R": dn / d01, "ci95_R": ci95, "ci90_R": ci90,
                    "decision": wa.decide(dn / d01, ci95, ci90, den_ci95=wa.ci(b01, 0.95))["label"],
                    "reaches_sesoi": bool(dn / d01 >= wa.SESOI)}
    coefs: dict[str, dict] = {}
    for f in sorted(d.glob("fold[0-9].json")):
        sc = json.loads(f.read_text(encoding="utf-8")).get("stage2_coefs", {})
        for stack, v in sc.items():
            c = coefs.setdefault(stack, {"names": v["names"], "by_fold": []})
            c["by_fold"].append(v["coef"])
    res["coefs"] = {}
    for stack, c in coefs.items():
        a = np.array(c["by_fold"], dtype=np.float64)
        res["coefs"][stack] = {nm: {"mean": float(a[:, i].mean()), "min": float(a[:, i].min()),
                                    "max": float(a[:, i].max())} for i, nm in enumerate(c["names"])}
    return res


def c2_detail(m: dict | None) -> dict | str:
    if not m:
        return "missing"
    pl = m.get("placebo") or {}
    return {"brier": {a: m["brier"][a] for a in ("M0", "M0s", "M1", "M1s", "M1c", "M3n", "M3np", "M3s", "M3d")
                      if a in m["brier"]},
            "delta_bs_placebo": pl.get("delta_bs_M1c_M3np"), "ci90_delta": m["primary"].get("ci90_delta"),
            "wet_rate": m.get("wet_rate")}


def main() -> None:
    lexfix = metrics("world_v2_lexfix")
    res: dict = {"rule": "amendments B-3 (C2), B-4, B-5 items 1, 2, 3, 5, 7",
                 "lexfix": {s: brief(metrics("world_v2_lexfix", s))
                            for s in ("primary", "country", "temporal", "llmcut")},
                 "frozen": {s: brief(metrics("world_v2", s)) for s in ("primary", "country", "temporal", "llmcut")},
                 "fold_seeds": fold_seeds(),
                 "permutations": permutations(lexfix["primary"]["delta_bs_M0_M1"] if lexfix else None),
                 "oracle_and_coefs": oracle_and_coefs(),
                 # 照合: 置換・係数の実行の M1c → M3n は、同じコードの確証の再実行（lexfix）と同じになるはず
                 "check_R_lexfix_primary": lexfix["primary"]["R"] if lexfix else None,
                 "label_variants": {v: brief(metrics(f"world_v2_lv_{v}")) for v in ("notrace", "nonwater")},
                 "gfm_c2": brief(metrics("world_v2_gfm")),
                 # C2 は R が未定義になりうる（物差しの区間が 0 を含む）ので、腕ごとの Brier と placebo の増分も残す
                 "gfm_c2_detail": c2_detail(metrics("world_v2_gfm"))}
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    for k in ("lexfix", "label_variants"):
        for s, v in res[k].items():
            print(k, s, v if isinstance(v, str) else f"R {v['R']:+.4f} 90% {v['ci90_R']} {v['decision']}")
    print("gfm_c2", res["gfm_c2"] if isinstance(res["gfm_c2"], str) else
          f"R {res['gfm_c2']['R']:+.4f} {res['gfm_c2']['ci90_R']} {res['gfm_c2']['decision']}")
    for name, blk in res["fold_seeds"].items():
        print("fold seeds", name, {k: blk.get(k) for k in ("n_available", "R_min", "R_max", "R_sd",
                                                           "combined_ci90", "n_equivalent_null")})
    print("permutations", res["permutations"], "| lexfix R", res["check_R_lexfix_primary"])
    oc = res["oracle_and_coefs"]
    print("oracle", oc if isinstance(oc, str) else {k: oc[k] for k in ("oracle", "names")})


if __name__ == "__main__":
    main()
