"""本体解析の足場: 複数イベント × 複数分割系の評価を run 契約で出す。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/run_main_eval.py --tag main_v1 \
        --splits watershed-loro-1,event-loo-1,temporal-1

各分割系について fold ごとに M0/M1/M2/M3 を訓練し、プールホールドアウト予測で
二方向クラスタ bootstrap CI を出す。predictions/*.npz に cell 単位の生予測と
event/watershed/role を保存するので、ラベル role 別の感度分析は後から
（run_pilot の npz を読む検証ファイル方式で）可能。現データは全 unit が
robustness 級（GSI 推定段彩図）— 主ラベル（実績図）は Phase 1 の残課題。
"""

from __future__ import annotations

import argparse
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from topsem.analysis.cell_features import (
    discover_label_units,
    load_population,
    population_tag,
    river_cells,
)
from topsem.analysis.models import fit_arms, fit_arms_ensemble
from topsem.corpus.snapshot import load_snapshot
from topsem.evaluation.bootstrap import (
    two_way_anchored_test,
    two_way_cluster_bootstrap,
    two_way_wild_test,
    wild_cluster_test,
)
from topsem.evaluation.metrics import brier_score, log_loss, pr_auc, top_area_capture
from topsem.features.route_a import ElementIndex
from topsem.runs import RunContext, start_run, write_json
from topsem.splits.geo import exclude_near_holdout
from topsem.splits.manifest import FoldDef, SplitSpec, load_split


def build_dataset(units: list[dict], index: ElementIndex) -> pd.DataFrame:
    """全ラベル unit × 母集団 × DEM を窓格子まで接合した学習 table。"""
    pop = load_population()
    feat_cache: dict = {}
    emb_cache: dict = {}
    return pd.concat(
        [river_cells(u, pop, index, feat_cache, emb_cache) for u in units],
        ignore_index=True,
    )


def fold_row(fold: FoldDef, tr: pd.DataFrame, te: pd.DataFrame, arms: dict[str, np.ndarray]) -> dict:
    """1 fold の指標行。ΔBS は base − alt（正 = alt の利得）。"""
    y = te["y"].to_numpy()
    row = {
        "fold": fold.name,
        "n_cells": len(te),
        "wet_rate": float(y.mean()),
        "events": sorted(te["event"].unique().tolist()),
        "watersheds": sorted(te["watershed"].unique().tolist()),
        "n_train_events": len(tr["event"].unique()),
        "n_train_watersheds": len(tr["watershed"].unique()),
    }
    if row["n_train_events"] <= 1:
        # event-loo 系の現状: 訓練が単一イベント = 異イベント汎化の実質 n/a
        row["note"] = "trained on single event — cross-event generalization not tested"
    # 監査 B-2 開示: held-out cell の水系が訓練側にも出る比率（位置汎化でなく
    # イベント汎化であることをメトリクス側でも可視化）
    row["same_watershed_train_frac"] = float(
        te["watershed"].isin(set(tr["watershed"].tolist())).mean()
    )
    for arm, pa in arms.items():
        row[f"bs_{arm}"] = brier_score(y, pa)
        row[f"ll_{arm}"] = log_loss(y, pa)
        if y.sum() > 0:
            row[f"pr_auc_{arm}"] = pr_auc(y, pa)
            row[f"cap10_{arm}"] = top_area_capture(y, pa)
    row["delta_bs_M3_M1"] = row["bs_M1"] - row["bs_M3"]
    row["delta_bs_M3_M0"] = row["bs_M0"] - row["bs_M3"]
    return row


def post_pooled(
    ci: dict,
    pooled: dict[str, list],
    *,
    yy: np.ndarray,
    pe: np.ndarray,
    pw: np.ndarray,
) -> tuple[dict, dict]:
    """プール予測に対する常設後処理: M1r 差（A-4）・role 層別・水系 leave-one-out wild。

    pooled は腕予測配列に加えて "role" キー（label role 配列のリスト）を含む。
    """
    p1, p3 = np.concatenate(pooled["M1"]), np.concatenate(pooled["M3"])
    if "M1r" in pooled:  # A-4 発効時のみ: 河道距離だけの基礎腕との差
        ci["M3-M1r"] = two_way_cluster_bootstrap(
            yy, np.concatenate(pooled["M1r"]), np.concatenate(pooled["M3"]), pe, pw
        )
        # r3-F4: 地名の単独増分（地形抜きの世界）と地形の増分の分解
        ci["M1rN-M1r"] = two_way_cluster_bootstrap(
            yy, np.concatenate(pooled["M1r"]), np.concatenate(pooled["M1rN"]), pe, pw
        )
        ci["M3-M1rN"] = two_way_cluster_bootstrap(
            yy, np.concatenate(pooled["M1rN"]), np.concatenate(pooled["M3"]), pe, pw
        )
        ci["M1rN_wild"] = wild_cluster_test(
            yy, np.concatenate(pooled["M1r"]), np.concatenate(pooled["M1rN"]), pe, pw
        )
        ci["M1rN_wild2"] = two_way_wild_test(
            yy, np.concatenate(pooled["M1r"]), np.concatenate(pooled["M1rN"]), pe, pw
        )
        ci["M1rN_anchored"] = two_way_anchored_test(
            yy, np.concatenate(pooled["M1r"]), np.concatenate(pooled["M1rN"]), pe, pw
        )
    # label role 層別（snapshot_only と robustness を混ぜて語らないための主成果物行）
    pr = np.concatenate(pooled["role"])
    strata = {}
    for rname in sorted(set(pr.tolist())):
        m = pr == rname
        if not m.any():
            continue
        bs1, bs3 = brier_score(yy[m], p1[m]), brier_score(yy[m], p3[m])
        strata[rname] = {"n_cells": int(m.sum()), "wet_rate": float(yy[m].mean()),
                         "bs_M1": bs1, "bs_M3": bs3, "delta_bs_M3_M1": bs1 - bs3}
    # prereg §5 常設感度: 水系を一つ除いた wild p（wet 率外れ値水系への依存度）。
    # 再訓練は不要で、プールホールドアウト予測のセルを落とすだけ。
    leaveone = {}
    for wname in sorted(set(pw.tolist())):
        m = pw != wname
        if not m.any():
            continue
        blk = wild_cluster_test(yy[m], p1[m], p3[m], pe[m], pw[m])
        blk["n_cells"] = int(m.sum())
        leaveone[f"no_{wname}"] = blk
    return strata, leaveone


def run_split(sv: str, data: pd.DataFrame, run: RunContext, seeds: list[int], spec: SplitSpec) -> dict:
    """1 分割系を全 fold 回し、fold 表 + プール bootstrap を返す。"""
    leak = spec.validate(data["event"].tolist(), data["watershed"].tolist())
    if not leak["ok"]:
        msg = f"分割 manifest とデータが不一致: {leak['problems']}"
        raise RuntimeError(msg)
    (run.path("predictions") / sv).mkdir(parents=True, exist_ok=True)
    fold_rows = []
    pooled: dict[str, list] = {}
    pooled_y, pooled_e, pooled_w = [], [], []
    for fold in spec.folds:
        evs, wss = data["event"].tolist(), data["watershed"].tolist()
        mask = np.array(spec.fold_mask(fold, evs, wss))
        # temporal では訓練側も spec.train_mask が未来イベントを除外（レビュー A-1）
        tr = data[np.array(spec.train_mask(fold, evs, wss))]
        te = data[mask]
        n_geo_excl = 0
        if spec.geo_exclusion_m > 0 and spec.kind in ("watershed_holdout", "strict"):
            # A-7 (r4' F3): 隣接水系窓と物理地面を共有する訓練セルを幾何除去
            tr, n_geo_excl = exclude_near_holdout(tr, te, spec.geo_exclusion_m)
        if len(te) == 0 or te["y"].sum() == 0:
            continue
        arms = fit_arms(tr, te, seeds[0]) if len(seeds) == 1 else fit_arms_ensemble(tr, te, seeds)
        row = fold_row(fold, tr, te, arms)
        row["n_geo_excluded_train"] = n_geo_excl  # A-7: 幾何除去された訓練セル数
        fold_rows.append(row)
        for arm in arms:
            pooled.setdefault(arm, []).append(arms[arm])
        pooled_y.append(te["y"].to_numpy())
        pooled_e.append(te["event"].to_numpy())
        pooled_w.append(te["watershed"].to_numpy())
        pooled.setdefault("role", []).append(te["role"].to_numpy())
        np.savez_compressed(
            run.path("predictions") / sv / f"{fold.name}.npz",
            **{f"p_{k}": v for k, v in arms.items()},
            y=te["y"].to_numpy(),
            event=te["event"].to_numpy(),
            watershed=te["watershed"].to_numpy(),
            role=te["role"].to_numpy(),
            cell_id=te.index.to_numpy(),
        )
        print(f"{sv}/{fold.name:30s} ΔBS(M3-M1)={row['delta_bs_M3_M1']:+.6f} n={row['n_cells']}")
    if not pooled_y:
        msg = f"{sv}: ホールドアウトに wet cell が無い（fold 全滅）"
        raise RuntimeError(msg)
    yy = np.concatenate(pooled_y)
    pe = np.concatenate(pooled_e)
    pw = np.concatenate(pooled_w)
    # prereg §5 の wet 率外れ値水系を機械生成物として保存（r3-M6b: 表生成側が
    # 中央値再計算にフォールバックせず正本の flag を読めるように）
    wr = {str(w): float(yy[pw == w].mean()) for w in set(pw.tolist())}
    med = float(np.median(list(wr.values()))) if wr else 0.0
    outliers = sorted(w for w, v in wr.items() if v > 3 * med)
    ci = {
        "M3-M1": two_way_cluster_bootstrap(
            yy, np.concatenate(pooled["M1"]), np.concatenate(pooled["M3"]), pe, pw
        ),
        "M3-M0": two_way_cluster_bootstrap(
            yy, np.concatenate(pooled["M0"]), np.concatenate(pooled["M3"]), pe, pw
        ),
        # 小クラスタ系の正しい不確実表示（レビュー A-3・CGM4 p 値）
        # A-8: 主 = two_way_wild（加法重み・完全 studentized）。旧 product 重みは
        # 歴史比較用に残置（v6 数値の再現経路）。
        "M3-M1_wild": wild_cluster_test(
            yy, np.concatenate(pooled["M1"]), np.concatenate(pooled["M3"]), pe, pw
        ),
        "M3-M1_wild2": two_way_wild_test(
            yy, np.concatenate(pooled["M1"]), np.concatenate(pooled["M3"]), pe, pw
        ),
        # A-8 主判定: 既登録衝撃族 anchored（臨界値は data 非依存の閉形式）
        "M3-M1_anchored": two_way_anchored_test(
            yy, np.concatenate(pooled["M1"]), np.concatenate(pooled["M3"]), pe, pw
        ),
    }
    strata, leaveone = post_pooled(ci, pooled, yy=yy, pe=pe, pw=pw)
    # realized 節数（r3-F4: 「10–27 節」の幅ではなく実数を保存）
    realized = len({(e, w) for e, w in zip(pe.tolist(), pw.tolist(), strict=True)})
    return {
        "folds": fold_rows,
        "leakage_checks": leak,
        "bootstrap": ci,
        "strata_by_role": strata,
        "sensitivity_leaveone_wild": leaveone,
        "wet_rate_outliers": outliers,
        "n_sections_realized": realized,
        "n_cells_pooled": len(yy),
        "n_events_pooled": len(set(pe.tolist())),
        "n_watersheds_pooled": len(set(pw.tolist())),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--seed", type=int, default=42, help="単一シード（--seeds 指定時は無視）")
    ap.add_argument("--seeds", default=None, help="カンマ区切り = 予測平均アンサンブル腕")
    ap.add_argument("--splits", default="watershed-loro-1,event-loo-1,temporal-1")
    args = ap.parse_args()

    run = start_run(args.tag or datetime.now(UTC).strftime("main_%Y%m%dT%H%M%SZ"))
    snap = load_snapshot()
    index = ElementIndex.from_elements(snap.elements)
    units = discover_label_units()
    data = build_dataset(units, index)
    river_on = "dist_river" in data.columns

    split_versions = [s.strip() for s in args.splits.split(",") if s.strip()]
    seed_list = (
        [int(s) for s in args.seeds.split(",")] if args.seeds else [args.seed]
    )
    per_split = {sv: run_split(sv, data, run, seed_list, load_split(sv)) for sv in split_versions}
    summary = {
        **run.provenance(),
        "corpus_snapshot": snap.manifest(),
        # prereg v1.0 凍結（2026-09-25）。flag で A-4 発効状態が切り替わるため
        # 腕記述はデータ列から動的に作る（凍結記述との乖離を防ぐ）
        "protocol_version": "v1.0" + ("+A4" if river_on else ""),
        "splits_versions": split_versions,
        "arms": {
            "M0": "ix,iy（公開位置情報のみ）",
            "M1": "M0+地形(SRTM30) elev, 相対低地度"
            + (" +河道距離 log1p（A-4）" if river_on else "") + "【主対照】",
            "M2": "M1+辞書 route_a 3 特征（機構判別・副次）",
            "M3": "M1+n-gram hashed embedding 64 次元【主比較 M3 vs M1】",
            **({"M1r": "M0+河道距離 log1p のみ（A-4 消融去・水文統制）",
                "M1rN": "M0+河道距離+地名埋め込み（地形なし；r3-F4 の分解腕）"}
               if river_on else {}),
        },
        "seeds": seed_list,
        # r3-M2: 住所・地名層の版来歴（現行=単一最新版、retrospective 評価である
        # ことの機械可読な記録。point_in_time.md §追記参照）
        "population_source_tag": population_tag(),
        "ensemble": "mean-of-seeds predict_proba" if len(seed_list) > 1 else "single seed",
        "label_roles": sorted({str(u["role"]) for u in units}),
        "n_cells_total": len(data),
        "results": per_split,
    }
    write_json(run.path("metrics.json"), per_split)
    write_json(run.path("summary.json"), summary)
    print(f"→ {run.run_dir}")


if __name__ == "__main__":
    main()
