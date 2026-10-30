"""名前チャネル（preregistration_v2.md §3）: N-lex・N-learn・placebo。

- N-lex: 凍結辞書を最近傍 K=3 の P 名に照合（場所を知らない＝暗記不可能）
- N-learn: 文字 1–4 gram（ハッシュ 2^20）の L2 ロジスティック得点。**外 fold の訓練側だけ**で、
  場所グループ単位の内側 5 分割交差適合により訓練セルの得点を作り、評価セルは訓練側全体で
  学習したモデルで採点する。各名前の重みは合計 1（多くのセルを覆う名前が支配しない）。
  v30 の過学習（64 次元を勾配ブースティングへ直入れ）を避けるための 1 次元化
- placebo: AOI 内の P 名集合の中で名前を無作為置換（集落の位置・数は不変、意味だけ破壊）
- N-read（H4）: 言語モデルの語義判定（``topsem.world.nread``）を最近傍 K=3 名に割り当てる。
  placebo は N-lex と**同じ置換**を使う
"""

from __future__ import annotations

import zlib

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.utils import murmurhash3_32

from topsem.world.lexicon import Lexicon, cell_lexicon_features
from topsem.world.nread import NREAD_CATEGORIES

NREAD_HYDRO = ("water", "wetland", "riverine", "crossing", "lowland", "highland")

K = 3
GID_COLS = [f"P_gid{j}" for j in range(1, K + 1)]
DIST_COLS = [f"P_d{j}" for j in range(1, K + 1)]
_VEC = HashingVectorizer(analyzer="char_wb", ngram_range=(1, 4), n_features=2**20,
                         alternate_sign=False, norm="l2", lowercase=True)
INNER_FOLDS = 5
C_REG = 1.0
MIN_NAMES = 2  # 名前モデルを学習する最小の異なり名前数


def gid_lookup(gn: pd.DataFrame) -> dict[int, tuple[str, str]]:
    return dict(zip(gn["gid"].astype(np.int64).tolist(),
                    zip(gn["asciiname"].astype(str).tolist(), gn["cc"].astype(str).tolist(),
                        strict=True), strict=True))


def placebo_perms(df: pd.DataFrame, seed: int) -> dict[str, dict[int, int]]:
    """AOI ごとに、その AOI で参照される P 名の gid を置換する表（gid → 名前を借りる gid）。"""
    out: dict[str, dict[int, int]] = {}
    for aoi, sub in df.groupby("aoi", sort=True):
        gids = np.unique(sub[GID_COLS].to_numpy().ravel())
        gids = gids[gids >= 0]
        # 組込み hash() はプロセスごとに乱数化されるので使わない（再現性）
        rng = np.random.default_rng([seed, zlib.crc32(str(aoi).encode("utf-8"))])
        perm = gids[rng.permutation(len(gids))]
        out[str(aoi)] = {int(g): int(p) for g, p in zip(gids, perm, strict=True)}
    return out


def placebo_lookups(df: pd.DataFrame, lookup: dict[int, tuple[str, str]], seed: int,
                    ) -> dict[str, dict[int, tuple[str, str]]]:
    """AOI ごとに、その AOI で参照される P 名の間で名前を置換した lookup。"""
    return {aoi: {g: lookup.get(p, ("", "")) for g, p in m.items()}
            for aoi, m in placebo_perms(df, seed).items()}


def site_strata(gids: np.ndarray, value: dict[int, float]) -> np.ndarray:
    """事後 B-16 の C3: AOI の gid を、値（集落の地点の HAND）の AOI 内の三分位で層 0・1・2 に。

    分位点は 1/3・2/3、同じ値は低い層へ。値の無い gid は層 −1。
    """
    v = np.array([value.get(int(g), np.nan) for g in gids], dtype=np.float64)
    lab = np.full(len(gids), -1, dtype=np.int64)
    ok = np.isfinite(v)
    if ok.any():
        q1, q2 = np.quantile(v[ok], [1 / 3, 2 / 3])
        lab[ok] = np.where(v[ok] <= q1, 0, np.where(v[ok] <= q2, 1, 2))
    return lab


def placebo_perms_strata(df: pd.DataFrame, seed: int,
                         value: dict[int, float]) -> dict[str, dict[int, int]]:
    """事後 B-16 の C3: placebo_perms と同じ置換を、AOI の中の層（site_strata）の中に限ったもの。

    乱数は placebo_perms と同じ作り方に層の番号（−1〜2 に 1 を足したもの）を足す。
    """
    out: dict[str, dict[int, int]] = {}
    for aoi, sub in df.groupby("aoi", sort=True):
        gids = np.unique(sub[GID_COLS].to_numpy().ravel())
        gids = gids[gids >= 0]
        lab = site_strata(gids, value)
        m: dict[int, int] = {}
        for s in np.unique(lab):
            g = gids[lab == s]
            rng = np.random.default_rng([seed, zlib.crc32(str(aoi).encode("utf-8")), int(s) + 1])
            m |= {int(a): int(b) for a, b in zip(g, g[rng.permutation(len(g))], strict=True)}
        out[str(aoi)] = m
    return out


def placebo_lookups_strata(df: pd.DataFrame, lookup: dict[int, tuple[str, str]], seed: int,
                           value: dict[int, float]) -> dict[str, dict[int, tuple[str, str]]]:
    """placebo_lookups の層別の版（B-16 の C3）。"""
    return {aoi: {g: lookup.get(p, ("", "")) for g, p in m.items()}
            for aoi, m in placebo_perms_strata(df, seed, value).items()}


def lexicon_block(df: pd.DataFrame, lookup: dict[int, tuple[str, str]], lex: Lexicon,
                  *, per_aoi: dict[str, dict[int, tuple[str, str]]] | None = None,
                  prefix: str = "") -> pd.DataFrame:
    """N-lex 特徴（per_aoi を渡すと placebo 版）。"""
    if per_aoi is None:
        f = cell_lexicon_features(df[GID_COLS].to_numpy(), df[DIST_COLS].to_numpy(), lookup, lex)
        return pd.DataFrame({prefix + k: v for k, v in f.items()}, index=df.index)
    parts = []
    for aoi, sub in df.groupby("aoi", sort=False):
        f = cell_lexicon_features(sub[GID_COLS].to_numpy(), sub[DIST_COLS].to_numpy(),
                                  per_aoi[str(aoi)], lex)
        parts.append(pd.DataFrame({prefix + k: v for k, v in f.items()}, index=sub.index))
    return pd.concat(parts).loc[df.index]


def nearest_names(df: pd.DataFrame, lookup: dict[int, tuple[str, str]],
                  per_aoi: dict[str, dict[int, tuple[str, str]]] | None = None) -> np.ndarray:
    """各セルの最近傍 P 名の asciiname（placebo 時は AOI 内置換後）。"""
    g = df["P_gid1"].to_numpy()
    if per_aoi is None:
        return np.array([lookup.get(int(x), ("", ""))[0] if x >= 0 else "" for x in g],
                        dtype=object)
    out = np.empty(len(df), dtype=object)
    for aoi, idx in df.groupby("aoi", sort=False).indices.items():
        m = per_aoi[str(aoi)]
        out[idx] = [m.get(int(x), ("", ""))[0] if x >= 0 else "" for x in g[idx]]
    return out


def _fit_name_model(names: np.ndarray, y: np.ndarray) -> LogisticRegression | float:
    """名前単位に集約し、各名前の重みを合計 1 にして L2 ロジスティック回帰。"""
    frame = pd.DataFrame({"n": names, "y": y})
    agg = frame.groupby("n", sort=True)["y"].agg(["sum", "count"])
    agg = agg[agg.index != ""]
    if len(agg) < MIN_NAMES or agg["sum"].sum() == 0 or (agg["sum"] == agg["count"]).all():
        return float(np.mean(y)) if len(y) else 0.0
    p = (agg["sum"] / agg["count"]).to_numpy()
    txt = agg.index.to_numpy(dtype=object)
    x = _VEC.transform(np.concatenate([txt, txt]))
    yy = np.concatenate([np.ones(len(txt)), np.zeros(len(txt))])
    w = np.concatenate([p, 1.0 - p])
    keep = w > 0
    model = LogisticRegression(C=C_REG, max_iter=2000, solver="lbfgs")
    model.fit(x[keep], yy[keep], sample_weight=w[keep])
    return model


def _score(model: LogisticRegression | float, names: np.ndarray) -> np.ndarray:
    if isinstance(model, float):
        p = min(max(model, 1e-6), 1 - 1e-6)
        return np.full(len(names), np.log(p / (1 - p)), dtype=np.float32)
    # 同じ名前が多数のセルに出るので、一意な名前だけ変換して配る（結果は同一）
    uniq, inv = np.unique(np.asarray(names, dtype=str), return_inverse=True)
    return model.decision_function(_VEC.transform(uniq)).astype(np.float32)[inv]


def inner_folds(groups: np.ndarray, k: int = INNER_FOLDS) -> np.ndarray:
    """学習側の内側 k 分割（場所グループ単位、ソート順の巡回割当で決定的）。"""
    ug = np.unique(groups)
    inner = {g: i % k for i, g in enumerate(ug)}
    return np.array([inner[g] for g in groups])


def nlearn_scores(names_tr: np.ndarray, y_tr: np.ndarray, groups_tr: np.ndarray,
                  names_te: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """交差適合した訓練側得点と、訓練側全体で学習した評価側得点。"""
    fold_tr = inner_folds(groups_tr)
    oof = np.zeros(len(names_tr), dtype=np.float32)
    for k in range(INNER_FOLDS):
        m = fold_tr == k
        if not m.any():
            continue
        model = _fit_name_model(names_tr[~m], y_tr[~m])
        oof[m] = _score(model, names_tr[m])
    full = _fit_name_model(names_tr, y_tr)
    return oof, _score(full, names_te)


def nread_lookup(tab: pd.DataFrame) -> dict[int, tuple[str, float]]:
    """N-read の推論結果（gid, category, confidence, parsed_ok）→ {gid: (category, confidence)}。

    応答を解釈できなかった名前は登録しない（特徴上は「判定なし」= 被覆 0）。
    """
    ok = tab[tab["parsed_ok"].astype(bool)]
    return dict(zip(ok["gid"].astype(np.int64).tolist(),
                    zip(ok["category"].astype(str).tolist(),
                        ok["confidence"].astype(float).tolist(), strict=True), strict=True))


def nread_block(df: pd.DataFrame, nr: dict[int, tuple[str, float]], *,
                perms: dict[str, dict[int, int]] | None = None, prefix: str = "",
                ) -> pd.DataFrame:
    """N-read 特徴（§3）: 最近傍 1 名の 9 類 one-hot・確信度・被覆と、距離重み K 名の水文 6 類。

    重みは N-lex と同じ 1/(1 + d/1 km)。perms（placebo_perms）を渡すと placebo 版。
    """
    gids = df[GID_COLS].to_numpy(dtype=np.int64)
    if perms is not None:
        gids = gids.copy()
        for aoi, idx in df.groupby("aoi", sort=False).indices.items():
            m = perms[str(aoi)]
            sub = gids[idx]
            gids[idx] = np.array([[m.get(g, g) if g >= 0 else g for g in row]
                                  for row in sub.tolist()], dtype=np.int64).reshape(sub.shape)
    dists = df[DIST_COLS].to_numpy(dtype=np.float64)
    n, k = gids.shape
    cat_idx = {c: i for i, c in enumerate(NREAD_CATEGORIES)}
    uniq, inv = np.unique(np.asarray(gids, dtype=np.int64).ravel(), return_inverse=True)
    oh_u = np.zeros((len(uniq), len(NREAD_CATEGORIES)), dtype=np.float32)
    conf_u = np.zeros(len(uniq), dtype=np.float32)
    cov_u = np.zeros(len(uniq), dtype=np.float32)
    for i, g in enumerate(uniq.tolist()):
        hit = nr.get(g) if g >= 0 else None
        if hit is None:
            continue
        oh_u[i, cat_idx[hit[0]]] = 1.0
        conf_u[i] = hit[1]
        cov_u[i] = 1.0
    onehot = oh_u[inv].reshape(n, k, len(NREAD_CATEGORIES))
    conf = conf_u[inv].reshape(n, k)
    cov = cov_u[inv].reshape(n, k)
    w = 1.0 / (1.0 + np.where(np.isfinite(dists), dists, 1e9) / 1000.0)
    out: dict[str, np.ndarray] = {}
    for c in NREAD_CATEGORIES:
        out[f"nr_{c}_n1"] = onehot[:, 0, cat_idx[c]]
    for c in NREAD_HYDRO:
        out[f"nr_{c}_k"] = (onehot[:, :, cat_idx[c]] * w).sum(axis=1).astype(np.float32)
    out["nr_conf_n1"] = conf[:, 0]
    out["nr_covered_n1"] = cov[:, 0]
    return pd.DataFrame({prefix + key: v for key, v in out.items()}, index=df.index)


NREAD_FEATURES = [*(f"nr_{c}_n1" for c in NREAD_CATEGORIES),
                  *(f"nr_{c}_k" for c in NREAD_HYDRO), "nr_conf_n1", "nr_covered_n1"]


def top_ngrams(names: np.ndarray, y: np.ndarray, n: int = 30) -> list[tuple[str, float]]:
    """解釈用: 全データで学習した名前モデルの係数上位（候補 n-gram をハッシュで逆引き）。"""
    model = _fit_name_model(names, y)
    if isinstance(model, float):
        return []
    analyzer = _VEC.build_analyzer()
    cands: dict[int, str] = {}
    for s in pd.unique(names):
        for g in analyzer(str(s)):
            idx = abs(murmurhash3_32(g, seed=0)) % _VEC.n_features
            cands.setdefault(idx, g)
    coef = model.coef_.ravel()
    order = np.argsort(-coef)
    return [(cands.get(int(j), f"#{j}"), float(coef[j])) for j in order[:n]]
