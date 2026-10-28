r"""事後 B-8 の C2: 主比較の二段目を、場所グループの bootstrap の中で当て直す区間。

``world_eval --dump-stage2`` が保存した各 fold の二段目の設計（基礎の logit・地域構成・名前の特徴・場所グループ）を
読み、回 first .. first + n − 1 について、場所グループの多項の重み（種 20261028 + 回）で、各 fold の基準 M1c と
名前 M3n の二段目を学習側の重みで当て直し、評価側の重みで二乗誤差を合計する。基礎の予測（M0・M1）は固定。
同じ重みで、モデルを固定した（確証と同じ）予測の二乗誤差も合計する。

    uv run python packages/topsem/scripts/world_boot_refit.py --src world_v2_b8stage2 \\
        --tag world_v2_b8bootref --first 0 --n 10

出力: ``data/runs/<tag>/predictions/<split>/boot_<first>.npz``（回ごとの和）と同名の json（来歴）。
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from topsem.config import RUNS_DIR
from topsem.runs import start_run, write_json
from topsem.world.b8 import group_weights, stage2_weighted

SEED_BASE = 20261028
N_FOLDS = 10


def load_folds(src: str, split: str) -> list[dict]:
    d = RUNS_DIR / src / "predictions" / split
    folds = []
    for k in range(N_FOLDS):
        z = np.load(d / f"fold{k}_stage2.npz", allow_pickle=False)
        folds.append({key: z[key] for key in z.files})
    return folds


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default="world_v2_b8stage2")
    ap.add_argument("--tag", default="world_v2_b8bootref")
    ap.add_argument("--split", default="primary")
    ap.add_argument("--first", type=int, required=True)
    ap.add_argument("--n", type=int, required=True)
    args = ap.parse_args()
    t0 = time.time()
    run = start_run(args.tag)
    folds = load_folds(args.src, args.split)
    units = np.unique(np.concatenate([f["unit_te"] for f in folds]).astype(str))
    code = {u: i for i, u in enumerate(units)}
    g = len(units)
    prep = []
    check = []
    for f in folds:
        c_tr = np.array([code[u] for u in f["unit_tr"].astype(str)])
        c_te = np.array([code[u] for u in f["unit_te"].astype(str)])
        xn_tr = np.column_stack([f["xc_tr"], f["xn_tr"]])
        xn_te = np.column_stack([f["xc_te"], f["xn_te"]])
        pc, init_c = stage2_weighted(f["z_tr"], f["z_te"], f["xc_tr"], f["xc_te"], f["y_tr"])
        pn, init_n = stage2_weighted(f["z_tr"], f["z_te"], xn_tr, xn_te, f["y_tr"])
        # 重みの無い当て直しは確証の予測（保存した p_M1c・p_M3n）を再現するはず
        check.append([float(np.abs(pc - f["p_M1c_te"]).max()), float(np.abs(pn - f["p_M3n_te"]).max())])
        prep.append({"c_tr": c_tr, "c_te": c_te, "xn_tr": xn_tr, "xn_te": xn_te, "init_c": init_c,
                     "init_n": init_n, "y_te": f["y_te"].astype(np.float64)})
    print(f"loaded {len(folds)} folds, {g} groups; max |refit − saved| {np.max(check):.2e} "
          f"({time.time() - t0:.0f}s)", flush=True)
    names = ("refit_c", "refit_n", "fixed_c", "fixed_n", "m0", "m1", "w")
    acc = {k: np.zeros(args.n) for k in names}
    for i, b in enumerate(range(args.first, args.first + args.n)):
        wg = group_weights(g, SEED_BASE + b)
        for f, p in zip(folds, prep, strict=True):
            w_tr, w_te = wg[p["c_tr"]], wg[p["c_te"]]
            if w_te.sum() == 0:
                continue
            pc, _ = stage2_weighted(f["z_tr"], f["z_te"], f["xc_tr"], f["xc_te"], f["y_tr"], w_tr=w_tr,
                                    init=p["init_c"])
            pn, _ = stage2_weighted(f["z_tr"], f["z_te"], p["xn_tr"], p["xn_te"], f["y_tr"], w_tr=w_tr,
                                    init=p["init_n"])
            y = p["y_te"]
            acc["refit_c"][i] += float((w_te * (pc - y) ** 2).sum())
            acc["refit_n"][i] += float((w_te * (pn - y) ** 2).sum())
            acc["fixed_c"][i] += float((w_te * (f["p_M1c_te"] - y) ** 2).sum())
            acc["fixed_n"][i] += float((w_te * (f["p_M3n_te"] - y) ** 2).sum())
            acc["m0"][i] += float((w_te * (f["p_M0_te"] - y) ** 2).sum())
            acc["m1"][i] += float((w_te * (f["p_M1_te"] - y) ** 2).sum())
            acc["w"][i] += float(w_te.sum())
        print(f"replicate {b} done ({time.time() - t0:.0f}s)", flush=True)
    dst = run.path("predictions") / args.split / f"boot_{args.first:04d}.npz"
    dst.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dst, replicates=np.arange(args.first, args.first + args.n), **acc)
    write_json(dst.with_suffix(".json"),
               {**run.provenance(), "src": args.src, "split": args.split, "first": args.first, "n": args.n,
                "n_groups": g, "seed_base": SEED_BASE, "check_max_abs_refit_vs_saved": check,
                "rule": "amendment B-8 item C2", "sec": round(time.time() - t0)})
    print(f"done ({time.time() - t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
