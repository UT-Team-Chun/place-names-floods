# %% v40（事後・amendment B-2）: 勾配ブースティングの早期終了 — 主分割の 1 fold で M0・M1 は何回反復したか
"""学習器は scikit-learn の HistGradientBoostingClassifier（max_iter=200, learning_rate=0.1、他は既定）。
既定の early_stopping='auto' は訓練セルが 10,000 を超えると有効になり、訓練セルの無作為な 10% で損失が
10 回改善しなければ止まる。検証の 10% は場所でまとめていないので、位置だけの M0 が場所を覚える方向に
働き得る（内部審査の指摘 10）。確証解析と同じ表・同じ fold 0 の訓練側全体で、M0 と M1 を 5 シードで
fit し、反復回数（n_iter_）と、早期終了の有無を記録する。予測や判定には使わない。

クラスタで実行（確証解析と同じ amd64 のコンテナ）:
    world_submit.py --stage script --extra "verification/r_okauchi/v40_hgb_iterations.py"
出力: data/verification/v40_hgb_iterations.json（NAS）
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# コンテナ（/workspace）には CLAUDE.md が無いので、見つからなければ backend の 1 つ上をルートとみなす
ROOT = next(
    (p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file()),
    Path(__file__).resolve().parents[3],
)
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))

import world_eval as ev
from sklearn.ensemble import HistGradientBoostingClassifier

from topsem.config import DATA_DIR
from topsem.world.dataset import continent_map

OUT = DATA_DIR / "verification" / "v40_hgb_iterations.json"
FOLD = 0
SEEDS = (42, 43, 44, 45, 46)


def main() -> None:
    t0 = time.time()
    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    sens = ev.Sensitivity()
    dirs, _, infos = ev.select_aois(gn, sens)
    df = ev.build_table(dirs, seed=20261025, infos=infos, sens=sens)
    df, _, _ = ev.add_channels(df, gn, None)
    cont = continent_map(ev.ROOT / "geonames" / "countryInfo.txt")
    df["fold"] = ev.folds_for(
        "primary", df, infos, cont, buffer_km=sens.buffer_km
    ).to_numpy()
    tr = df[(df["fold"] != FOLD) & (df["fold"] >= 0)]
    res: dict = {"fold": FOLD, "n_train": len(tr), "arms": {}}
    for arm in ("M0", "M1"):
        x = tr[ev.ARMS[arm]].to_numpy(dtype=np.float32)
        y = tr["y"].to_numpy()
        rows = []
        for s in SEEDS:
            m = HistGradientBoostingClassifier(
                random_state=s, max_iter=200, learning_rate=0.1
            )
            m.fit(x, y)
            rows.append(
                {
                    "seed": s,
                    "n_iter": int(m.n_iter_),
                    "early_stopping": bool(m.do_early_stopping_),
                    "validation_fraction": m.validation_fraction,
                    "n_iter_no_change": m.n_iter_no_change,
                }
            )
        res["arms"][arm] = {
            "features": ev.ARMS[arm],
            "fits": rows,
            "n_iter": [r["n_iter"] for r in rows],
        }
        print(
            arm, [r["n_iter"] for r in rows], f"({time.time() - t0:.0f}s)", flush=True
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
