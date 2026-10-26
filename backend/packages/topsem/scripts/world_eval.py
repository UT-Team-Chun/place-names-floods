"""v2 全球検定の評価（preregistration_v2.md §4–§6）と合成ラベルの検出力シミュレーション（§8）。

1 プロセス = 1 分割系 × 1 fold（または全 fold）。予測は
``RUNS_DIR/<tag>/predictions/<split>/fold<k>.npz`` に書き、``world_aggregate.py`` がプールして
ΔBS・比 R・クラスタ bootstrap・TOST を出す。

合成モード（``--synthetic null|sesoi``）は**実ラベルを読まない**（AOI の実浸水率だけを使って
合成ラベルの切片を合わせる）。凍結前に許される唯一のモデル fit。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/world_eval.py --tag world_v2_primary --split primary --fold 3
    uv run python packages/topsem/scripts/world_eval.py --tag power_v2_null --synthetic null \
        --reps 0:20 --seeds 42
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression

from topsem.config import DATA_DIR, REPO_ROOT
from topsem.runs import start_run, write_json
from topsem.world.channels import (
    INNER_FOLDS,
    NREAD_FEATURES,
    gid_lookup,
    inner_folds,
    lexicon_block,
    nearest_names,
    nlearn_scores,
    nread_block,
    nread_lookup,
    placebo_lookups,
    placebo_perms,
)
from topsem.world.dataset import (
    GROUP_BUFFER_KM,
    GROUP_CAP,
    CoastIndex,
    aoi_info,
    assign_folds,
    continent_map,
    group_quotas,
    include,
    load_aoi_table,
    location_groups,
)
from topsem.world.lexicon import CATEGORIES, Lexicon

ROOT = DATA_DIR / "global"
LEXICON = REPO_ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"
TEMPORAL_CUT = "2024-01-01"
LLM_CUT = "2025-02-01"  # Gemma-4 の事前学習データ打切り（2025-01、モデルカード）の翌月初
CUTS = {"temporal": TEMPORAL_CUT, "llmcut": LLM_CUT}
PLACEBO_SEED = 20261025

POS = ["lon", "lat"]
HYDRO = ["ldr_all", "ldr_major"]
# 土地被覆は洪水の前年版の ESA CCI / C3S（§2.3、point-in-time 準拠）。WorldCover 2021 の lc_* は
# 2012–2020 年の洪水には事後情報なので使わない
LC = ["cci_built", "cci_crop", "cci_tree", "cci_shrub", "cci_grass", "cci_bare", "cci_wetland",
      "cci_water"]
SETTLE = ["lP_d1", "P_n1km", "P_n3km"]
TER_1 = ["elev_1km"]
TER_2 = [*TER_1, "elev_mean", "depress"]
TER_3 = [*TER_2, "elev_min", "elev_std", "slope_mean"]
TER_4 = [*TER_3, "hand_mean", "hand_min"]
LEX = [*(f"lex_{c}_n1" for c in CATEGORIES), *(f"lex_{c}_k" for c in CATEGORIES),
       "lex_wetside_n1", "lex_covered_n1"]
NAMES = [*LEX, "nlearn"]
NAMES_P = [*(f"p_{c}" for c in LEX), "nlearn_p"]
NREAD_P = [f"p_{c}" for c in NREAD_FEATURES]
T0 = [*POS, *HYDRO, *LC, *SETTLE]
ARMS: dict[str, list[str]] = {
    "M0": POS,
    "M0N": [*POS, *NAMES],
    "T0": T0, "T0N": [*T0, *NAMES],
    "T1": [*T0, *TER_1], "T1N": [*T0, *TER_1, *NAMES],
    "T2": [*T0, *TER_2], "T2N": [*T0, *TER_2, *NAMES],
    "T3": [*T0, *TER_3], "T3N": [*T0, *TER_3, *NAMES],
    "M1": [*T0, *TER_4],
    "M3": [*T0, *TER_4, *NAMES],
    "M3p": [*T0, *TER_4, *NAMES_P],
    "M3x": [*T0, *TER_4, *LEX],
    "M3l": [*T0, *TER_4, "nlearn"],
    "M3r": [*T0, *TER_4, *NREAD_FEATURES],
    "M3rp": [*T0, *TER_4, *NREAD_P],
}
NREAD_ARMS = ("M3r", "M3rp")  # --nread を与えたときだけ計算（H4）
# 主比較（§4・§6）: 入れ子の増分モデル。基礎モデル（上の HistGB 腕）の予測 — 学習側は内側
# 交差適合、評価側は学習側全体で fit — のロジットを二段目のロジスティック回帰に入れ、名前の有無を
# 比べる。別々に学習した柔軟モデル同士の比較は、モデルの揺らぎが名前の効果を覆う（pilot2:
# ΔBS(M1→M3) の揺らぎは地形の増分の最大 29%）ので副にする。
# さらに pilot3 で、名前を足した改善の大半は「その地域にどんな名前が並ぶか」（言語・地域の構成が
# 地域ごとの浸水率を拾う）で、置換 placebo でも同じだけ出ると分かった。主比較は名前の地域構成
# （AOI 平均、c_*）を基準側にも入れ、名前が場所ごとに持つ情報（主張の中身）を測る。
COMP = [f"c_{f}" for f in NAMES]
COMP_NR = [f"c_{f}" for f in NREAD_FEATURES]
STACK: dict[str, tuple[str, list[str]]] = {
    # 主: 場所内の増分（基準 = 基礎予測 + 名前の地域構成）
    "M1c": ("M1", COMP), "M3n": ("M1", [*COMP, *NAMES]), "M3np": ("M1", [*COMP, *NAMES_P]),
    "M3nx": ("M1", [*COMP, *LEX]), "M3nl": ("M1", [*COMP, "nlearn"]),
    "M0c": ("M0", COMP), "M0n": ("M0", [*COMP, *NAMES]),
    "T0c": ("T0", COMP), "T0n": ("T0", [*COMP, *NAMES]),
    "T1c": ("T1", COMP), "T1n": ("T1", [*COMP, *NAMES]),
    "T2c": ("T2", COMP), "T2n": ("T2", [*COMP, *NAMES]),
    "T3c": ("T3", COMP), "T3n": ("T3", [*COMP, *NAMES]),
    "M1rc": ("M1", COMP_NR), "M3rn": ("M1", [*COMP_NR, *NREAD_FEATURES]),
    "M3rnp": ("M1", [*COMP_NR, *NREAD_P]),
    # 副: 地域構成込みの全体の増分（基準 = 基礎予測の再較正だけ）
    "M1s": ("M1", []), "M3s": ("M1", NAMES), "M3ps": ("M1", NAMES_P),
    "M3rs": ("M1", NREAD_FEATURES), "M3rps": ("M1", NREAD_P),
}
NREAD_STACK = ("M1rc", "M3rn", "M3rnp", "M3rs", "M3rps")
# §8: 主判定（M1c → M3n、placebo M3np）と副の全体増分（M1s → M3s、M3ps）、分母（M0 → M1）
SIM_ARMS = ["M0", "M1"]
SIM_STACK = ["M1c", "M3n", "M3np", "M1s", "M3s", "M3ps"]
LOGIT_EPS = 1e-5


def _p_names_inside(gn: pd.DataFrame, b: tuple[float, float, float, float]) -> int:
    w, s, e, n = b
    m = (gn["fclass"] == "P") & gn["lon"].between(w, e) & gn["lat"].between(s, n)
    return int(m.sum())


def select_aois(gn: pd.DataFrame, sens: "Sensitivity | None" = None,  # noqa: UP037
                ) -> tuple[list[Path], dict, dict]:
    """採用 AOI と、除外理由・AOI 要約（ラベル単独）。sens は §7 の感度解析の設定。"""
    sens = sens or Sensitivity()
    dirs = sorted(p for p in (ROOT / "labels").glob("EMSR*")
                  if all((p / f).is_file() for f in ("cells.parquet", "features.parquet",
                                                       "landcover.parquet")))
    infos, reasons, keep = {}, Counter(), []
    coast = CoastIndex(ROOT / "coast" / "ne_10m_coastline.zip")
    for d in dirs:
        info = aoi_info(d, coast, all_types=sens.all_types)
        ok, why = include(info, n_names_inside=_p_names_inside(gn, info.bounds))
        if ok and sens.exclude_no_footprint and (info.footprint_missing
                                                 or info.footprint_implicit):
            ok, why = False, "sens:no-footprint"
        if ok and sens.min_event_year and int(info.event_time[:4]) < sens.min_event_year:
            ok, why = False, f"sens:before-{sens.min_event_year}"
        reasons[why] += 1
        infos[d.name] = info
        if ok:
            keep.append(d)
    return keep, dict(reasons), infos


@dataclass(frozen=True)
class Sensitivity:
    """§7 の感度解析の設定（既定 = 主解析）。"""

    all_types: bool = False  # §7-1 全浸水種別
    group_cap: int = GROUP_CAP  # §7-2 場所グループ上限（0 = 上限なし）
    aoi_cap: int = 0  # §7-2 AOI あたり上限（> 0 ならグループ上限の代わり）
    name_class: str = "P"  # §7-3 名前チャネルの地名クラス（P / PLTS）
    exclude_no_footprint: bool = False  # §7-4
    buffer_km: float = GROUP_BUFFER_KM  # §7-7
    min_event_year: int = 0  # §7-8


def quotas_for(dirs: list[Path], infos: dict, sens: Sensitivity) -> dict[str, int]:
    """AOI ごとの抽出セル数（§2.4-4。感度 §7-2 では上限なし / AOI あたり上限）。"""
    n_eval = {d.name: infos[d.name].n_eval for d in dirs}
    if sens.aoi_cap > 0:
        return {a: min(n, sens.aoi_cap) for a, n in n_eval.items()}
    if sens.group_cap <= 0:
        return n_eval
    grp = location_groups({d.name: infos[d.name].bounds for d in dirs}, sens.buffer_km)
    return group_quotas(n_eval, grp, sens.group_cap)


def build_table(dirs: list[Path], seed: int, infos: dict,
                sens: Sensitivity | None = None) -> pd.DataFrame:
    """採用 AOI の評価セル表（§2.4-4 の上限・§7 の感度設定つき）。"""
    sens = sens or Sensitivity()
    name_cols = [f"{c}{j}" for c in ("P_gid", "P_d") for j in (1, 2, 3)]
    plts_cols = [c.replace("P_", "PLTS_") for c in name_cols] if sens.name_class == "PLTS" else []
    cols = ["y", "aoi", "lon", "lat", "elev_1km", "elev_mean", "depress", "elev_min",
            "elev_std", "slope_mean", "hand_mean", "hand_min", "dist_river_all",
            "dist_river_major", *LC, *name_cols, *plts_cols,
            "P_n1km", "P_n3km", "f_riverine", "f_other", "f_flood_any"]
    quota = quotas_for(dirs, infos, sens)
    parts = [load_aoi_table(d, cap=quota[d.name], seed=seed, all_types=sens.all_types,
                            unspecified_as_riverine=infos[d.name].unspecified_as_riverine)[cols]
             for d in dirs]
    df = pd.concat(parts, ignore_index=True)
    df["ldr_all"] = np.log1p(df["dist_river_all"]).astype(np.float32)
    df["ldr_major"] = np.log1p(df["dist_river_major"]).astype(np.float32)
    df["lP_d1"] = np.log1p(np.minimum(df["P_d1"], 1e6)).astype(np.float32)
    if plts_cols:  # §7-3: 名前チャネルだけ P+L+T+S に（集落幾何 lP_d1・P_n* は P のまま）
        for c, pc in zip(name_cols, plts_cols, strict=True):
            df[c] = df[pc]
    return df


def add_channels(df: pd.DataFrame, gn: pd.DataFrame, nread: pd.DataFrame | None = None,
                 ) -> tuple[pd.DataFrame, dict, dict]:
    lex = Lexicon(LEXICON)
    look = gid_lookup(gn)
    per_aoi = placebo_lookups(df, look, PLACEBO_SEED)
    real = lexicon_block(df, look, lex)
    plac = lexicon_block(df, look, lex, per_aoi=per_aoi, prefix="p_")
    blocks = [df, real, plac]
    if nread is not None:  # H4: N-lex と同じ置換で placebo を作る
        nr = nread_lookup(nread)
        perms = placebo_perms(df, PLACEBO_SEED)
        blocks += [nread_block(df, nr), nread_block(df, nr, perms=perms, prefix="p_")]
    df = pd.concat(blocks, axis=1)
    df["name_real"] = nearest_names(df, look)
    df["name_plac"] = nearest_names(df, look, per_aoi)
    df["cc"] = [look.get(int(g), ("", ""))[1] if g >= 0 else "" for g in df["P_gid1"]]
    return df, {"sha256": lex.sha256, "version": lex.version}, per_aoi


def time_folds(aois: list[str], infos: dict, cut: str,
               unit: dict[str, str]) -> dict[str, int]:
    """時間外し系（§5）の fold: 評価 0 = cut 以後、訓練 1 = cut 以前、除外 −1。

    訓練は cut 以前のうち評価 AOI と同じ場所グループに属さないものに限る。同じ場所の
    以前の洪水は両側から外す（場所の暗記を時間外しでも防ぐ）。
    """
    post = {a for a in aois if infos[a].event_time[:10] >= cut}
    post_units = {unit[a] for a in post}
    return {a: 0 if a in post else (-1 if unit[a] in post_units else 1) for a in aois}


def country_units(aois: list[str], grp: dict[str, int], aoi_cc: dict[str, str],
                  n_cells: dict[str, int]) -> dict[str, str]:
    """国外しの単位: 場所グループごとにセル数で多数派の国（国境をまたぐ場所を分けない）。"""
    votes: dict[int, Counter] = {}
    for a in aois:
        votes.setdefault(grp[a], Counter())[aoi_cc[a]] += n_cells[a]
    return {a: votes[grp[a]].most_common(1)[0][0] for a in aois}


def folds_for(split: str, df: pd.DataFrame, infos: dict, cont: dict[str, str],
              buffer_km: float = GROUP_BUFFER_KM) -> pd.Series:
    """AOI → fold（§5）。fold −1 はどちらの側にも使わない（時間外し系の同一場所）。

    buffer_km は場所グループの緩衝幅（§7-6 の感度では 10 km）。引数で明示的に渡す: 以前は
    DataFrame.attrs で渡しており、pandas 3 の concat が attrs を捨てるため感度 s6 の fold と
    bootstrap の単位が 5 km のままになっていた（amendment B-1）。
    """
    aois = sorted(df["aoi"].unique())
    n_cells = df.groupby("aoi").size().to_dict()
    aoi_cc = df.groupby("aoi")["cc"].agg(lambda s: Counter(s[s != ""]).most_common(1)[0][0]
                                         if (s != "").any() else "??").to_dict()
    grp = location_groups({a: infos[a].bounds for a in aois}, buffer_km)
    if split == "country":
        unit = country_units(aois, grp, aoi_cc, n_cells)
    else:
        unit = {a: f"g{grp[a]:05d}" for a in aois}
    if split in CUTS:
        df["unit"] = df["aoi"].map(unit)
        return df["aoi"].map(time_folds(aois, infos, CUTS[split], unit))
    size = Counter()
    ucont: dict[str, Counter] = {}
    for a in aois:
        size[unit[a]] += n_cells[a]
        ucont.setdefault(unit[a], Counter())[cont.get(aoi_cc[a], "??")] += n_cells[a]
    sort = {u: (ucont[u].most_common(1)[0][0], -size[u], u) for u in size}
    f = assign_folds(unit, sort)
    df["unit"] = df["aoi"].map(unit)
    return df["aoi"].map(f)


def fit_arm(tr: pd.DataFrame, te: pd.DataFrame, feats: list[str], seeds: list[int]) -> np.ndarray:
    xtr = tr[feats].to_numpy(dtype=np.float32)
    xte = te[feats].to_numpy(dtype=np.float32)
    acc = np.zeros(len(te), dtype=np.float64)
    for s in seeds:
        m = HistGradientBoostingClassifier(random_state=s, max_iter=200, learning_rate=0.1)
        m.fit(xtr, tr["y"].to_numpy())
        acc += m.predict_proba(xte)[:, 1]
    return (acc / len(seeds)).astype(np.float32)


def _logit(p: np.ndarray) -> np.ndarray:
    q = np.clip(p.astype(np.float64), LOGIT_EPS, 1 - LOGIT_EPS)
    return np.log(q / (1 - q))


def crossfit_base(tr: pd.DataFrame, te: pd.DataFrame, feats: list[str], seeds: list[int],
                  fold_tr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """基礎モデルの、学習側の内側交差適合予測（場所グループ単位）と評価側予測（学習側全体で fit）。"""
    oof = np.zeros(len(tr), dtype=np.float32)
    for k in range(INNER_FOLDS):
        m = fold_tr == k
        if m.any():
            oof[m] = fit_arm(tr[~m], tr[m], feats, seeds)
    return oof, fit_arm(tr, te, feats, seeds)


def stage2(z_tr: np.ndarray, z_te: np.ndarray, x_tr: np.ndarray, x_te: np.ndarray,
           y_tr: np.ndarray) -> np.ndarray:
    """二段目のロジスティック回帰（L2、C = 1）。

    入力は [基礎予測のロジット, 名前特徴] を学習側の平均・標準偏差で標準化したもの。
    名前特徴が空なら基礎予測の再較正だけになる（入れ子の比較基準）。
    """
    a = np.column_stack([z_tr, x_tr]).astype(np.float64)
    b = np.column_stack([z_te, x_te]).astype(np.float64)
    mu = a.mean(axis=0)
    sd = a.std(axis=0)
    sd[sd == 0] = 1.0
    m = LogisticRegression(C=1.0, max_iter=1000)
    m.fit((a - mu) / sd, y_tr)
    return m.predict_proba((b - mu) / sd)[:, 1].astype(np.float32)


def add_composition(tr: pd.DataFrame, te: pd.DataFrame) -> None:
    """名前の地域構成: 名前特徴（N-learn 得点を含む）の AOI 平均 c_*（ラベルは使わない）。"""
    feats = [f for f in (*NAMES, *NREAD_FEATURES) if f in tr.columns]
    for frame in (tr, te):
        means = frame.groupby("aoi")[feats].transform("mean")
        for f in feats:
            frame[f"c_{f}"] = means[f].to_numpy(dtype=np.float32)


def run_fold(df: pd.DataFrame, fold: int, arms: list[str], seeds: list[int],
             stacks: tuple[str, ...] | list[str] = ()) -> dict:
    tr = df[(df["fold"] != fold) & (df["fold"] >= 0)].copy()
    te = df[df["fold"] == fold].copy()
    if split_is_temporal(df):
        tr = df[df["fold"] == 1].copy()
        te = df[df["fold"] == 0].copy()
    groups_tr = tr["unit"].to_numpy() if "unit" in tr else tr["aoi"].to_numpy()
    for src, col in (("name_real", "nlearn"), ("name_plac", "nlearn_p")):
        oof, sc = nlearn_scores(tr[src].to_numpy(), tr["y"].to_numpy(), groups_tr,
                                te[src].to_numpy())
        tr[col] = oof
        te[col] = sc
    add_composition(tr, te)
    fold_tr = inner_folds(groups_tr)
    preds: dict[str, np.ndarray] = {}
    z: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for b in sorted({STACK[s][0] for s in stacks}):  # 基礎の全体 fit は柔軟腕の予測も兼ねる
        oof, p_te = crossfit_base(tr, te, ARMS[b], seeds, fold_tr)
        preds[b] = p_te
        z[b] = (_logit(oof), _logit(p_te))
    for a in arms:
        if a not in preds:
            preds[a] = fit_arm(tr, te, ARMS[a], seeds)
    y_tr = tr["y"].to_numpy()
    for s in stacks:
        b, extra = STACK[s]
        x_tr = tr[extra].to_numpy(dtype=np.float64) if extra else np.empty((len(tr), 0))
        x_te = te[extra].to_numpy(dtype=np.float64) if extra else np.empty((len(te), 0))
        preds[s] = stage2(z[b][0], z[b][1], x_tr, x_te, y_tr)
    return {"te": te, "preds": preds, "n_train": len(tr)}


def strata_columns(te: pd.DataFrame) -> dict[str, np.ndarray]:
    """H5（記述のみ）と診断図に使う評価セルの層別変数（ラベルを含まない共変量だけ）。"""
    lex_any = te[[f"lex_{c}_n1" for c in CATEGORIES]].to_numpy().max(axis=1)
    return {"s_lex_any": lex_any.astype(np.int8),
            "s_wetside": te["lex_wetside_n1"].to_numpy(dtype=np.int8),
            "s_slope": te["slope_mean"].to_numpy(dtype=np.float32),
            "s_built": te["cci_built"].to_numpy(dtype=np.float32),
            "s_hand": te["hand_min"].to_numpy(dtype=np.float32),
            "lon": te["lon"].to_numpy(dtype=np.float32),
            "lat": te["lat"].to_numpy(dtype=np.float32)}


PROTOCOL = REPO_ROOT / "docs" / "protocol" / "preregistration_v2.md"
FROZEN_LINE = re.compile(r"^状態: \*\*FROZEN\*\* \d{4}-\d{2}-\d{2}T[\d:]+Z / commit [0-9a-f]{7,40}",
                         re.MULTILINE)


def protocol_frozen(path: Path = PROTOCOL) -> bool:
    """事前登録 v2 の状態行が FROZEN か（実ラベルでの評価は凍結後にだけ許す、§9）。"""
    return path.is_file() and bool(FROZEN_LINE.search(path.read_text(encoding="utf-8")))


def _sha256(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def split_is_temporal(df: pd.DataFrame) -> bool:
    return bool(df.attrs.get("split") in CUTS)


def synth_labels(df: pd.DataFrame, delta: float, rep: int) -> np.ndarray:
    """§8: 実特徴 + AOI 実浸水率で合成ラベル（名前効果 delta は N-lex wetside に載せる）。"""
    rng = np.random.default_rng([20261025, rep])

    def z(x: pd.Series) -> np.ndarray:
        v = x.to_numpy(dtype=float)
        v = np.where(np.isfinite(v), v, np.nanmedian(v))
        return (v - v.mean()) / (v.std() + 1e-9)

    lin = (-1.2 * z(np.log1p(df["hand_min"].clip(lower=0)))
           - 0.4 * z(df["slope_mean"]) - 0.4 * z(df["ldr_all"]) + 0.2 * z(df["cci_crop"]))
    field = np.zeros(len(df))
    for idx in df.groupby("aoi", sort=False).indices.values():
        sub = df.iloc[idx]
        bx = (sub["lon"].to_numpy() * 50).astype(int)
        by = (sub["lat"].to_numpy() * 50).astype(int)
        keys = pd.Series(list(zip(bx, by, strict=True)))
        codes, uniq = pd.factorize(keys)
        field[idx] = rng.normal(0, 0.8, len(uniq))[codes]
    nm = df["lex_wetside_n1"].to_numpy(dtype=float)
    lin = lin + field + delta * (nm - nm.mean())
    y = np.zeros(len(df), dtype=np.int8)
    for idx in df.groupby("aoi", sort=False).indices.values():
        target = float(df["y"].to_numpy()[idx].mean())  # 実浸水率（AOI 単位の集計値のみ使用）
        lo, hi = -15.0, 15.0
        for _ in range(40):
            mid = (lo + hi) / 2
            if (1 / (1 + np.exp(-(lin[idx] + mid)))).mean() < target:
                lo = mid
            else:
                hi = mid
        p = 1 / (1 + np.exp(-(lin[idx] + (lo + hi) / 2)))
        y[idx] = (rng.random(len(idx)) < p).astype(np.int8)
    return y


def main() -> None:  # noqa: PLR0915
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--split", default="primary",
                    choices=["primary", "country", "temporal", "llmcut"])
    ap.add_argument("--fold", type=int, default=-1, help="-1 = 全 fold")
    ap.add_argument("--seeds", default="42,43,44,45,46")
    ap.add_argument("--arms", default=",".join(ARMS), help="柔軟モデル（HistGB）の腕（副）")
    ap.add_argument("--stacks", default=",".join(STACK), help="入れ子の増分モデルの腕（主）")
    ap.add_argument("--group-cap", type=int, default=GROUP_CAP,
                    help="場所グループあたりの評価セル上限（§2.4-4、0 = 上限なし §7-2）")
    ap.add_argument("--aoi-cap", type=int, default=0, help="§7-2: AOI あたり上限（> 0 で有効）")
    ap.add_argument("--flood-types", default="riverine", choices=["riverine", "all"],
                    help="§7-1: 全浸水種別のラベル")
    ap.add_argument("--name-class", default="P", choices=["P", "PLTS"], help="§7-3")
    ap.add_argument("--exclude-no-footprint", action="store_true", help="§7-4")
    ap.add_argument("--group-buffer-km", type=float, default=GROUP_BUFFER_KM, help="§7-7")
    ap.add_argument("--min-event-year", type=int, default=0, help="§7-8: これより前の事象を除外")
    ap.add_argument("--synthetic", default="none", choices=["none", "null", "sesoi"])
    ap.add_argument("--delta", type=float, default=0.0, help="sesoi 時の名前効果（logit）")
    ap.add_argument("--reps", default="0:1")
    ap.add_argument("--nread", default="", help="N-read 推論結果 parquet（与えると M3r/M3rp も回す）")
    args = ap.parse_args()
    if args.synthetic == "none" and not protocol_frozen():
        msg = ("preregistration_v2.md is not FROZEN: real-label evaluation is allowed only after "
               "the freeze (synthetic labels only before it)")
        raise SystemExit(msg)

    t0 = time.time()
    run = start_run(args.tag)
    gn = pd.read_parquet(ROOT / "geonames" / "aoi_subset.parquet")
    sens = Sensitivity(all_types=args.flood_types == "all", group_cap=args.group_cap,
                       aoi_cap=args.aoi_cap, name_class=args.name_class,
                       exclude_no_footprint=args.exclude_no_footprint,
                       buffer_km=args.group_buffer_km, min_event_year=args.min_event_year)
    dirs, reasons, infos = select_aois(gn, sens)
    print(f"AOIs included {len(dirs)} / reasons {reasons}", flush=True)
    df = build_table(dirs, seed=20261025, infos=infos, sens=sens)
    nread = pd.read_parquet(args.nread) if args.nread else None
    df, lexinfo, _ = add_channels(df, gn, nread)
    cont = continent_map(ROOT / "geonames" / "countryInfo.txt")
    df.attrs["split"] = args.split
    df["fold"] = folds_for(args.split, df, infos, cont, buffer_km=sens.buffer_km).to_numpy()
    seeds = [int(s) for s in args.seeds.split(",")]
    print(f"table {len(df):,} cells / {df['aoi'].nunique()} AOIs ({time.time() - t0:.0f}s)",
          flush=True)

    folds = [0] if args.split in CUTS else (
        list(range(int(df["fold"].max()) + 1)) if args.fold < 0 else [args.fold])
    common = {**run.provenance(), "protocol": "preregistration_v2 (DRAFT until frozen)",
              "lexicon": lexinfo, "group_cap": args.group_cap, "seeds": seeds,
              "sensitivity": dataclasses.asdict(sens),
              "aoi_inclusion": reasons, "n_aoi": len(dirs), "split": args.split,
              "nread": {"path": args.nread, "sha256": _sha256(args.nread)} if args.nread else None,
              "n_aoi_by_fold": {int(k): int(v) for k, v in
                                df.groupby("fold")["aoi"].nunique().items()}}

    if args.synthetic != "none":
        a, b = (int(v) for v in args.reps.split(":"))
        out = []
        for rep in range(a, b):
            df["y"] = synth_labels(df, args.delta if args.synthetic == "sesoi" else 0.0, rep)
            ys, ps, us = [], {k: [] for k in (*SIM_ARMS, *SIM_STACK)}, []
            for k in folds:
                if not (df["fold"] == k).any():
                    continue
                r = run_fold(df, k, SIM_ARMS, seeds, SIM_STACK)
                ys.append(r["te"]["y"].to_numpy())
                us.append(r["te"]["unit"].to_numpy())
                for arm, acc in ps.items():
                    acc.append(r["preds"][arm])
            res = {"rep": rep, "y": np.concatenate(ys), "unit": np.concatenate(us).astype(str),
                   **{arm: np.concatenate(v) for arm, v in ps.items()}}
            path = run.path("sim") / f"rep{rep:04d}.npz"
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, **res)
            out.append(rep)
            print(f"rep {rep} done ({time.time() - t0:.0f}s)", flush=True)
        write_json(run.path(f"sim_meta_{a}_{b}.json"),
                   {**common, "synthetic": args.synthetic, "delta": args.delta, "reps": out})
        return

    arms = [x for x in args.arms.split(",") if x and (nread is not None or x not in NREAD_ARMS)]
    stacks = [x for x in args.stacks.split(",")
              if x and (nread is not None or x not in NREAD_STACK)]
    for k in folds:
        if not (df["fold"] == k).any():
            continue
        r = run_fold(df, k, arms, seeds, stacks)
        te = r["te"]
        dst = run.path("predictions") / args.split / f"fold{k}.npz"
        dst.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            dst, y=te["y"].to_numpy(), aoi=te["aoi"].to_numpy(dtype=str),
            unit=te["unit"].to_numpy(dtype=str) if "unit" in te else te["aoi"].to_numpy(dtype=str),
            cc=te["cc"].to_numpy(dtype=str),
            **strata_columns(te),
            **{f"p_{a}": v for a, v in r["preds"].items()},
        )
        write_json(run.path("predictions") / args.split / f"fold{k}.json",
                   {**common, "fold": k, "n_test": len(te), "n_train": r["n_train"],
                    "arms": {a: ARMS[a] for a in arms},
                    "stacks": {s: {"base": STACK[s][0], "extra": STACK[s][1]} for s in stacks}})
        print(f"fold {k}: test {len(te):,} ({time.time() - t0:.0f}s)", flush=True)
    # 配分 0 件の AOI（巨大な場所グループ内の小さな AOI）は表に行が無いので、行のある AOI だけ書く
    fold_of = df.groupby("aoi")["fold"].first()
    n_of = df.groupby("aoi").size()
    write_json(run.path(f"aoi_table_{args.split}.json"),
               {a: {"fold": int(fold_of[a]), "n_cells": int(n_of[a])} for a in fold_of.index})
    print(json.dumps({"done": args.tag, "sec": round(time.time() - t0)}))


if __name__ == "__main__":
    main()
