"""事後 B-8（初見の査読 第 4 回への対応）の部品。

計算と規則は amendment B-8 で計算の前に固定したもの。

- B3: 決着しなかった判定に要る場所グループの数
- B5(b): 場所グループを一つずつ除いた推定からの jackknife の標準誤差（CV3）
- B6(d)・C5: 掃除した辞書（2 文字の接辞と英語の 6 つの接尾辞を除く）
- B6(f): 誤分類の補正（稀な結果の近似）
- B12: 言語に依らない水名の読み
- C2: 二段目を bootstrap の重みで当て直す
"""

from __future__ import annotations

import copy
import re
from collections import Counter
from collections.abc import Iterable

import numpy as np
from sklearn.linear_model import LogisticRegression

SHORT_AFFIX_MAX = 2  # B6(d): この長さ以下の語頭・語末の項目を除く（語は残す）
# B6(d): 移植名・人名と衝突しうる英語の接尾辞
EN_TRANSFER_SUFFIXES = ("burn", "bourne", "borne", "pool", "well", "wells")
HYDRONYM_RADIUS_M = 5000.0  # B12
HYDRONYM_MIN_CORE = 4
HYDRONYM_MIN_REST = 2
HYDRONYM_TOP_GENERIC = 200  # B12 (i)(ii): 最も多いトークンを総称として除く数
HYDRONYM_CODES = ("STM", "LK")  # B12: GeoNames の H 類の符号の頭（川・小川・湖）
_TOKEN = re.compile(r"[a-z]+")


# %% B3: 要る場所グループの数
def groups_needed(est: float, lo: float, hi: float, n_groups: int, *,
                  margin: float | None = None) -> float:
    """要る場所グループの数（区間の半幅が数の平方根に反比例し、点推定が変わらないとして）。

    margin を与えると「区間が ±margin に入る」、与えなければ「区間の下限が 0 を超える」
    に要る数。届かなければ無限大。lo・hi は判定に使う水準の区間（同等性は 90%、正は 95%）。
    """
    half = (hi - lo) / 2.0
    room = (margin - abs(est)) if margin is not None else est
    if half <= 0:
        return float(n_groups)
    if room <= 0:
        return float("inf")
    return float(n_groups) * (half / room) ** 2


# %% B5(b): jackknife（CV3）
def jackknife_se(loo: Iterable[float]) -> float:
    """一つずつ除いた推定 β_(g) からの jackknife の標準誤差。

    sqrt((G − 1)/G Σ(β_(g) − β̄)²)。
    """
    b = np.asarray(list(loo), dtype=np.float64)
    g = len(b)
    if g < 2:  # noqa: PLR2004
        return float("nan")
    return float(np.sqrt((g - 1) / g * ((b - b.mean()) ** 2).sum()))


# %% B6(f): 誤分類の補正
def correct_misclassification(or_obs: np.ndarray | float, p1: np.ndarray | float,
                              p0: np.ndarray | float) -> np.ndarray:
    """非差別的な誤分類の補正（稀な結果の近似）。

    OR = (OR*(1 − p0) − (1 − p1)) / (p1 − OR* p0)。p1: 水系側と読んだ名のうち真に水に
    関わる割合、p0: 他の名のうち真に水に関わる割合。分子か分母が 0 以下なら NaN
    （観測の値が誤分類で説明しきれない）。
    """
    o = np.asarray(or_obs, dtype=np.float64)
    a = np.asarray(p1, dtype=np.float64)
    b = np.asarray(p0, dtype=np.float64)
    num = o * (1 - b) - (1 - a)
    den = a - o * b
    with np.errstate(divide="ignore", invalid="ignore"):
        out = num / den
    return np.where((num > 0) & (den > 0), out, np.nan)


# %% B6(d)・C5: 掃除した辞書
def _drop(form: str, pos: str, lang: str) -> bool:
    if pos in ("prefix", "suffix") and len(form) <= SHORT_AFFIX_MAX:
        return True
    return lang == "en" and pos == "suffix" and form in EN_TRANSFER_SUFFIXES


def clean_lexicon_doc(doc: dict) -> tuple[dict, list[tuple[str, str, str]]]:
    """B6(d) の形を語頭・語末の位置だけ除いた辞書の文書と、除いた (言語, 形, 位置)。

    位置の一覧を持つ項目は位置ごとに分ける（語の位置の形は全て残す）。元の文書は変えない。
    """
    out = copy.deepcopy(doc)
    dropped: list[tuple[str, str, str]] = []
    for lang, spec in (out.get("languages") or {}).items():
        new = []
        for e in spec.get("entries") or []:
            poss = e.get("pos", "word")
            poss = [poss] if isinstance(poss, str) else list(poss)
            forms = e["form"] if isinstance(e["form"], list) else [e["form"]]
            for pos in poss:
                keep = []
                for f in forms:
                    if _drop(str(f).lower(), pos, lang):
                        dropped.append((lang, str(f).lower(), pos))
                    else:
                        keep.append(f)
                if keep:
                    new.append({"form": keep, "pos": pos, "cat": e["cat"]})
        spec["entries"] = new
    return out, dropped


# %% B12: 水名の読み
def name_tokens(name: str) -> list[str]:
    """名前を小文字の英字のトークンに分ける（GeoNames の asciiname を想定）。"""
    return _TOKEN.findall(str(name).lower())


def frequent_tokens(names: Iterable[str], top: int = HYDRONYM_TOP_GENERIC) -> set[str]:
    """名前の集まりで最も多いトークン（名前ごとに 1 回数える）の上位 top 個。"""
    c: Counter[str] = Counter()
    for n in names:
        c.update(set(name_tokens(n)))
    return {t for t, _ in c.most_common(top)}


def hydronym_cores(name: str, generic: set[str], frequent_local: set[str]) -> set[str]:
    """水の地物の名の核。

    4 文字以上で、総称（generic）にも、その国の集落の名で多いトークンにも入らないトークン。
    """
    return {t for t in name_tokens(name)
            if len(t) >= HYDRONYM_MIN_CORE and t not in generic and t not in frequent_local}


def contains_hydronym(settlement_name: str, cores: Iterable[str]) -> bool:
    """集落の名のトークンが核と一致するか、核で始まるか終わり残りが 2 文字以上なら真。"""
    cs = set(cores)
    if not cs:
        return False
    for t in name_tokens(settlement_name):
        if t in cs:
            return True
        for c in cs:
            if len(t) - len(c) >= HYDRONYM_MIN_REST and (t.startswith(c) or t.endswith(c)):
                return True
    return False


# %% C2: 二段目を bootstrap の重みで当て直す
def group_weights(n_groups: int, seed: int) -> np.ndarray:
    """場所グループの多項の重み（合計 n_groups、各確率 1/n_groups）。"""
    rng = np.random.default_rng(seed)
    return rng.multinomial(n_groups, np.full(n_groups, 1.0 / n_groups)).astype(np.float64)


def stage2_weighted(  # noqa: PLR0913
    z_tr: np.ndarray, z_te: np.ndarray, x_tr: np.ndarray, x_te: np.ndarray, y_tr: np.ndarray, *,
    w_tr: np.ndarray | None = None, init: tuple[np.ndarray, float] | None = None,
) -> tuple[np.ndarray, tuple[np.ndarray, float]]:
    """確証の二段目（L2、C = 1、学習側で標準化）を、学習側の重みつきで当てる。

    重みが無ければ確証の stage2 と同じ計算。重みは複製した行と同じ扱い（標準化も重みつきの
    平均・標準偏差）。init（標準化した尺度の係数・切片）を与えると、そこから lbfgs を始める
    （結果は収束先で、初期値には依らない）。
    """
    a = np.column_stack([z_tr, x_tr]).astype(np.float64)
    b = np.column_stack([z_te, x_te]).astype(np.float64)
    w = np.ones(len(a)) if w_tr is None else np.asarray(w_tr, dtype=np.float64)
    keep = w > 0
    a, y, w = a[keep], np.asarray(y_tr)[keep], w[keep]
    mu = np.average(a, axis=0, weights=w)
    sd = np.sqrt(np.average((a - mu) ** 2, axis=0, weights=w))
    sd[sd == 0] = 1.0
    m = LogisticRegression(C=1.0, max_iter=1000, warm_start=init is not None)
    if init is not None:
        m.coef_ = np.asarray(init[0], dtype=np.float64).reshape(1, -1).copy()
        m.intercept_ = np.array([float(init[1])])
    m.fit((a - mu) / sd, y, sample_weight=None if w_tr is None else w)
    return m.predict_proba((b - mu) / sd)[:, 1], (m.coef_[0].copy(), float(m.intercept_[0]))
