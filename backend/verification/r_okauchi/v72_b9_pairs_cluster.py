# %% v72（事後・amendment B-13 の追記 7）: B9 の (b)(c)(g) の pairs の再抽出をクラスタで
"""v70b の B9 と同じ部分集合を作り、v68.pairs_or を同じ鍵で呼ぶ（保存は data/verification/v68_pairs_<鍵>.npy）。

確証判定は変えない。表・部分集合・種・回数は v70b の b9_near と同じ（部分集合の式はそこから写した）。手元の v70b の B9 は
この保存を読む（同じ鍵の保存があれば再計算しない）。

使い方（backend/ から）: python verification/r_okauchi/v72_b9_pairs_cluster.py --part b --workers 12
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import sys
import time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--part", required=True, choices=["b", "c", "g"])
ap.add_argument("--workers", type=int, default=12)
args = ap.parse_args()
os.environ["V68_WORKERS"] = str(args.workers)  # v68 の作業プロセスの数（読み込む前に決める）
sys.argv = [sys.argv[0]]  # v70b・v68 の読み込みに引数を渡さない

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v70b = _load("v70b_b13_part2")
v68, v69 = v70b.v68, v70b.v69


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    src = v68.VER / "v68_uncapped.parquet"
    print(f"uncapped table {src} sha256 {_sha256(src)}", flush=True)
    u = v70b._uncapped()
    t0 = time.time()
    if args.part == "b":  # v70b の b9_near と同じ: EMSR052・EMSR130 を除く
        act = np.array([a.split("_")[0] for a in u["aoi"].to_numpy(dtype=str)])
        keep2 = ~np.isin(act, ("EMSR052", "EMSR130"))
        o = v68.pairs_or(u.loc[keep2].reset_index(drop=True), "v70_b9b_no052_130")
    elif args.part == "c":  # 世代 0・1 を除く
        gen = {a: int(v70b._json(v69.LABELS / a / "meta.json").get("schema_generation", -1))
               for a in u["aoi"].unique()}
        ug = u["aoi"].map(gen).to_numpy()
        o = v68.pairs_or(u.loc[~np.isin(ug, v70b.GEN_EXCLUDE)].reset_index(drop=True), "v70_b9c_gen2plus")
    else:  # 地物の符号 PPLX・PPLH・PPLQ・PPLW を除く
        fc = v70b._gn()["fcode"].reindex(u["P_gid1"].to_numpy(np.int64)).astype(str).to_numpy()
        o = v68.pairs_or(u.loc[~np.isin(fc, v70b.FCODE_EXCLUDE)].reset_index(drop=True), "v70_b9g_fcodes")
    print(f"part {args.part} ({time.time() - t0:.0f} s): OR {o['odds_ratio']:.4f} pairs 95% {o['pairs_ci95']} "
          f"reps ok {o['pairs_reps_ok']}", flush=True)


if __name__ == "__main__":
    main()
