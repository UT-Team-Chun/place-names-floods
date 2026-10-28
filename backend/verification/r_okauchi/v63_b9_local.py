# %% v63（事後・amendment B-9 の B3〜B8）: 既存の成果物から手元で計算するもの
"""確証判定は変えない。B-9 で計算の前に固定した量を、直した実行（world_v2_lexfix）の held-out の予測から計算する。

- B3: M1c と M3n の Brier の Murphy の分解（信頼性・分解能・不確実性。B-8 の B8 と同じ 10 の分位の区間）と、
  その差の場所グループの cluster bootstrap（4,000 回、種 20261029）の 95% 区間
- B4: AOI の中だけの評価: 各 AOI で予測の logit を一つずらして平均の予測を観測の浸水率に合わせた M0・M1・M1c・M3n
  の R_within = ΔBS(M1c → M3n) / ΔBS(M0 → M1) と cluster bootstrap の 90%・95% 区間
- B5: 浸水の地図の技能（各 AOI で浸水したセルの数だけ予測の高い順に地図にする）: CSI・的中率・空振り率の AOI の
  中央値。HAND だけの比較は HAND の最小、同じなら平均、それも同じなら乱数の順。CSI の差（M3n − M1c）の AOI の
  中央値と cluster bootstrap の 95% 区間
- B6: 集落の近くの推定（B-8 の B6 の主の式、水系側 × 500 m、EMS）を、水系側の名の近くに浸水セルのある activation を
  一つずつ除いて
- B7: 6 区分の係数（500 m 以内、B-8 の B6 の区分のモデル）の Benjamini–Hochberg の q
- B8: 解析した AOI のうち、凍結した目録に詳細図の製品があるもの・概観図と両方あるものの数

出力: stdout と docs/verification/v63_b9_local.json
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "backend").is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem"))
sys.path.insert(0, str(ROOT / "backend" / "packages" / "topsem" / "scripts"))
sys.path.insert(0, str(HERE))

OUT = ROOT / "docs" / "verification" / "v63_b9_local.json"
VER = ROOT / "docs" / "verification"
INVENTORY = ROOT / "docs" / "data" / "ems_inventory_20261025.json"
SEED = 20261029
B = 4000
NEAR_M = 500.0
N_BINS = 10


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v55 = _load("v55_b7_local")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _boot_counts(units: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """場所グループの番号と、B 回の再標本での各グループの回数 [B, G]（多項の重み）。"""
    lab, inv = np.unique(units, return_inverse=True)
    g = len(lab)
    draws = rng.integers(0, g, size=(B, g))
    cnt = np.zeros((B, g), dtype=np.float64)
    for b in range(B):
        cnt[b] = np.bincount(draws[b], minlength=g)
    return inv, cnt


# %% B3: Murphy の分解
def _bins(p: np.ndarray) -> np.ndarray:
    edges = np.quantile(p, np.linspace(0, 1, N_BINS + 1))
    return np.clip(np.searchsorted(edges, p, side="right") - 1, 0, N_BINS - 1)


def _decomp(n: np.ndarray, sp: np.ndarray, sy: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """区間ごとの数・予測の和・浸水の和（[..., 区間]）から信頼性・分解能・不確実性。"""
    nt = n.sum(axis=-1)
    ybar = sy.sum(axis=-1) / nt
    pb = sp / np.maximum(n, 1e-12)
    yb = sy / np.maximum(n, 1e-12)
    rel = (n * (pb - yb) ** 2).sum(axis=-1) / nt
    res = (n * (yb - ybar[..., None]) ** 2).sum(axis=-1) / nt
    return rel, res, ybar * (1 - ybar)


def b3_murphy(t: pd.DataFrame, inv: np.ndarray, cnt: np.ndarray) -> dict:
    y = t["y"].to_numpy(np.float64)
    g = cnt.shape[1]
    out: dict = {"bins": N_BINS, "replicates": B, "seed": SEED}
    per = {}
    for m in ("M1c", "M3n"):
        p = t[f"p_lf_{m}"].to_numpy(np.float64)
        b = _bins(p)
        key = inv * N_BINS + b
        n = np.bincount(key, minlength=g * N_BINS).reshape(g, N_BINS)
        sp = np.bincount(key, weights=p, minlength=g * N_BINS).reshape(g, N_BINS)
        sy = np.bincount(key, weights=y, minlength=g * N_BINS).reshape(g, N_BINS)
        rel0, res0, unc0 = _decomp(n.sum(0), sp.sum(0), sy.sum(0))
        rel_b, res_b, _ = _decomp(cnt @ n, cnt @ sp, cnt @ sy)
        per[m] = (rel_b, res_b)
        out[m] = {"brier": float(np.mean((p - y) ** 2)), "reliability": float(rel0), "resolution": float(res0),
                  "uncertainty": float(unc0)}
    d_rel = per["M3n"][0] - per["M1c"][0]
    d_res = per["M3n"][1] - per["M1c"][1]
    out["diff_M3n_minus_M1c"] = {
        "reliability": out["M3n"]["reliability"] - out["M1c"]["reliability"],
        "reliability_ci95": np.quantile(d_rel, [0.025, 0.975]).tolist(),
        "resolution": out["M3n"]["resolution"] - out["M1c"]["resolution"],
        "resolution_ci95": np.quantile(d_res, [0.025, 0.975]).tolist(),
        "brier": out["M3n"]["brier"] - out["M1c"]["brier"]}
    return out


# %% B4: AOI の中だけの評価
def _shift_to_rate(z: np.ndarray, aoi: np.ndarray, rate: np.ndarray, n_a: np.ndarray) -> np.ndarray:
    """AOI ごとに logit を c だけずらして、平均の予測を観測の浸水率に合わせる（Newton、全 AOI を同時に）。"""
    na = len(n_a)
    c = np.zeros(na)
    for _ in range(100):
        p = 1 / (1 + np.exp(-(z + c[aoi])))
        f = np.bincount(aoi, weights=p, minlength=na) / n_a - rate
        d = np.bincount(aoi, weights=p * (1 - p), minlength=na) / n_a
        step = f / np.maximum(d, 1e-12)
        c -= np.clip(step, -5, 5)
        if np.abs(f).max() < 1e-10:
            break
    return 1 / (1 + np.exp(-(z + c[aoi])))


def b4_within(t: pd.DataFrame, inv: np.ndarray, cnt: np.ndarray) -> dict:
    y = t["y"].to_numpy(np.float64)
    aoi = pd.factorize(t["aoi"])[0]
    n_a = np.bincount(aoi).astype(np.float64)
    rate = np.bincount(aoi, weights=y) / n_a
    inner = (rate > 0) & (rate < 1)
    g = cnt.shape[1]
    sse = {}
    for m in ("M0", "M1", "M1c", "M3n"):
        p = np.clip(t[f"p_lf_{m}"].to_numpy(np.float64), 1e-7, 1 - 1e-7)
        q = _shift_to_rate(np.log(p / (1 - p)), aoi, rate, n_a)
        q = np.where(inner[aoi], q, rate[aoi])  # 全て浸水・全て乾いた AOI は誤差 0
        sse[m] = np.bincount(inv, weights=(q - y) ** 2, minlength=g)
    n_u = np.bincount(inv, minlength=g).astype(np.float64)
    nn = cnt @ n_u

    def r_of(w: np.ndarray) -> np.ndarray:
        den = (w @ sse["M0"] - w @ sse["M1"])
        return (w @ sse["M1c"] - w @ sse["M3n"]) / den

    one = np.ones(g)
    r0 = float(r_of(one[None, :])[0])
    rb = r_of(cnt)
    yard = float((sse["M0"].sum() - sse["M1"].sum()) / n_u.sum())
    d = float((sse["M1c"].sum() - sse["M3n"].sum()) / n_u.sum())
    db = (cnt @ sse["M1c"] - cnt @ sse["M3n"]) / nn
    return {"R_within": r0, "ci90_R": np.quantile(rb, [0.05, 0.95]).tolist(),
            "ci95_R": np.quantile(rb, [0.025, 0.975]).tolist(), "yardstick_within": yard,
            "delta_bs_within": d, "ci95_delta_bs": np.quantile(db, [0.025, 0.975]).tolist(),
            "n_aoi_inner": int(inner.sum()), "replicates": B}


# %% B5: 浸水の地図の技能
def _csi_rows(y: np.ndarray, order: np.ndarray, k: int) -> tuple[float, float, float]:
    tp = int(y[order[:k]].sum())
    return tp / (2 * k - tp), tp / k, (k - tp) / k


def b5_skill(t: pd.DataFrame, rng: np.random.Generator) -> dict:
    rows: dict[str, list] = {"M1": [], "M1c": [], "M3n": [], "HAND": []}
    aoi_unit, d_csi = [], []
    for _, x in t.groupby("aoi", sort=True):
        y = x["y"].to_numpy()
        k = int(y.sum())
        if k in (0, len(y)):
            continue
        res = {}
        for m in ("M1", "M1c", "M3n"):
            res[m] = _csi_rows(y, np.argsort(-x[f"p_lf_{m}"].to_numpy(np.float64), kind="stable"), k)
        u = rng.random(len(y))
        order = np.lexsort((u, x["hand_mean"].to_numpy(np.float64), x["hand_min"].to_numpy(np.float64)))
        res["HAND"] = _csi_rows(y, order, k)
        for m, r in res.items():
            rows[m].append(r)
        d_csi.append(res["M3n"][0] - res["M1c"][0])
        aoi_unit.append(x["unit"].iloc[0])
    out: dict = {}
    for m, r in rows.items():
        a = np.array(r)
        out[m] = {s: {"median": float(np.median(a[:, i])), "q25_q75": np.quantile(a[:, i], [0.25, 0.75]).tolist()}
                  for i, s in enumerate(("csi", "hit_rate", "false_alarm_ratio"))}
        out[m]["n_aoi"] = len(a)
    d = np.array(d_csi)
    lab, inv = np.unique(np.array(aoi_unit), return_inverse=True)
    g = len(lab)
    meds = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, g, size=g)
        w = np.bincount(pick, minlength=g)[inv]
        meds[b] = np.median(np.repeat(d, w)) if w.sum() else np.nan
    out["delta_csi_M3n_minus_M1c"] = {"median": float(np.median(d)), "mean": float(d.mean()),
                                      "ci95_median": np.nanquantile(meds, [0.025, 0.975]).tolist(),
                                      "share_aoi_positive": float(np.mean(d > 0)),
                                      "share_aoi_negative": float(np.mean(d < 0)), "replicates": B}
    out["rule"] = "HAND ties broken by mean HAND, then at random (seed 20261029)"
    return out


# %% B6: activation を一つずつ除く
def b6_leave_activation(t: pd.DataFrame) -> dict:
    near = (t["P_d1"].to_numpy() <= NEAR_M) & (t["lex_wetside_n1"].to_numpy() > 0) & (t["y"].to_numpy() == 1)
    acts = sorted(t.loc[near, "activation"].unique())
    rows = []
    for a in acts:
        o = v55.odds(t, cls="wetside", mask=(t["activation"] != a).to_numpy())
        rows.append({"activation": a, "n_flooded_near": int((near & (t["activation"] == a).to_numpy()).sum()),
                     "odds_ratio": o["odds_ratio"], "ci95": o["ci95"]})
    ors = np.array([r["odds_ratio"] for r in rows])
    lo, hi = rows[int(np.argmin(ors))], rows[int(np.argmax(ors))]
    return {"n_activations": len(rows), "min": lo, "max": hi, "rows": rows}


# %% B7: 区分の BH
def b7_bh() -> dict:
    ca = _json(VER / "v57_b8_local.json")["B6"]["class_adjusted"]
    names = list(ca)
    p = np.array([ca[c]["p"] for c in names])
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        prev = min(prev, p[i] * m / rank)
        q[i] = prev
    return {c: {"p": float(p[i]), "q_bh": float(q[i]), "p_holm": float(ca[c]["p_holm"])} for i, c in enumerate(names)}


# %% B8: 詳細図の製品
def b8_detail(aois: list[str]) -> dict:
    inv = _json(INVENTORY)
    kinds: dict[str, set] = {}
    for a in inv["activations"]:
        for p in a.get("del_products", []):
            if p.get("kind") != "PRODUCT":
                continue
            key = f"{a['code']}_AOI{p['aoi']}"
            kinds.setdefault(key, set()).add("detail" if re.search(r"DETAIL", p["url"]) else "overview")
    have = {a: kinds.get(a, set()) for a in aois}
    return {"n_aoi": len(aois), "with_detail": sum("detail" in v for v in have.values()),
            "with_detail_and_overview": sum({"detail", "overview"} <= v for v in have.values()),
            "detail_only": sum(v == {"detail"} for v in have.values()),
            "inventory_aois_with_detail": sum("detail" in v for v in kinds.values()),
            "inventory_aois_with_both": sum({"detail", "overview"} <= v for v in kinds.values())}


# %% A7: 方法の細部の数（記述）
def a7_counts() -> dict:
    """50 m の副格子の AOI の数、河川以外の浸水で除いたセルの数（上限の前の評価できるセル）、二人目の符号化者のモデル。"""
    import world_eval as ev

    from topsem.world.dataset import evaluable, riverine_labels

    gn = pd.read_parquet(ev.ROOT / "geonames" / "aoi_subset.parquet")
    dirs, _, infos = ev.select_aois(gn)
    n50, ex_other, n_ev = 0, 0, 0
    for d in dirs:
        meta = _json(d / "meta.json")
        n50 += int(float(meta["grid"]["fine_res_m"]) >= 50)
        c = pd.read_parquet(d / "cells.parquet")
        _, keep = riverine_labels(c, unspecified_as_riverine=infos[d.name].unspecified_as_riverine)
        ev_ = evaluable(c)
        ex_other += int((ev_ & ~keep).sum())
        n_ev += int(ev_.sum())
    out = {"n_aoi": len(dirs), "n_aoi_50m_subgrid": n50, "n_cells_excluded_other_flood": ex_other,
           "n_cells_evaluable_before_cap": n_ev}
    # 二人目の符号化者（v57b）の Codex のセッションの記録（手元の ~/.codex、2026-10-27 の v57b の依頼文を含むもの）
    models, versions, n = set(), set(), 0
    for f in sorted((Path.home() / ".codex" / "sessions" / "2026" / "10" / "27").glob("rollout-*.jsonl")):
        txt = f.read_text(encoding="utf-8", errors="ignore")
        if "coding what Wikipedia says about the ORIGIN of place names" not in txt:
            continue
        n += 1
        models |= set(re.findall(r'"model":"([^"]+)"', txt))
        versions |= set(re.findall(r'"cli_version":"([^"]+)"', txt))
    out["coder2"] = {"sessions": n, "models": sorted(models), "cli_versions": sorted(versions),
                     "date": "2026-10-27"}
    return out


def main() -> None:
    if "--a7" in sys.argv:  # 方法の細部の数だけを足す
        res = _json(OUT)
        res["A7"] = a7_counts()
        print("A7", res["A7"], flush=True)
        OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
        return
    t = v55.table()
    rng = np.random.default_rng(SEED)
    inv, cnt = _boot_counts(t["unit"].to_numpy(), rng)
    res: dict = {"rule": "amendment B-9 items B3-B8 (local, from the corrected run world_v2_lexfix)",
                 "seed": SEED, "replicates": B}
    res["B3"] = b3_murphy(t, inv, cnt)
    print("B3", json.dumps(res["B3"], default=float)[:600], flush=True)
    res["B4"] = b4_within(t, inv, cnt)
    print("B4", res["B4"], flush=True)
    res["B5"] = b5_skill(t, np.random.default_rng(SEED))
    print("B5", json.dumps({k: v for k, v in res["B5"].items()}, default=float)[:900], flush=True)
    res["B6"] = b6_leave_activation(t)
    print("B6", {k: v for k, v in res["B6"].items() if k != "rows"}, flush=True)
    res["B7"] = b7_bh()
    print("B7", res["B7"], flush=True)
    res["B8"] = b8_detail(sorted(t["aoi"].unique()))
    print("B8", res["B8"], flush=True)
    OUT.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")


if __name__ == "__main__":
    main()
