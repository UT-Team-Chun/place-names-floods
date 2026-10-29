"""amendment B-13（初見の査読 第 9 回、事後）の実装。

- B1: 比の余白に対する Fieller 型の同等性の検定
  （L± = ΔBS ∓ 0.05·Δ_base の片側の限界）
- B15: 相対の修飾語の対（同じ基の名で、下と上の修飾語だけが違う集落の対）
- C1(b): 再起の記録を、その事象より前の事象だけから
- C2: 誤りを直した辞書（査読者 3 が名指した項目と、分けるよう求めた三つの類）
- C3: 基礎のモデルの反復の回数を、場所グループを分けた内側の検証で選ぶ
"""

from __future__ import annotations

import re
from datetime import UTC, datetime

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

from topsem.world.b8 import clean_lexicon_doc

SESOI_R = 0.05

# %% C2: 誤りを直した辞書
# (b) 同形異義。位置が None なら全ての位置で外す
HOMOGRAPHS: tuple[tuple[str, str, tuple[str, ...] | None], ...] = (
    ("de", "ach", ("suffix",)), ("de", "ried", None), ("fr", "mont", ("prefix",)),
    ("no", "ore", None), ("no", "ora", None), ("fa", "kacha", None), ("my", "taung", None),
    ("pt", "monte", None), ("es", "monte", None), ("en", "bank", None), ("en", "banks", None),
)
# (c) 洪水の意味を持たないと名指された川沿いの項目
RIVERINE_NAMED: dict[str, tuple[str, ...]] = {
    "sv": ("ang", "angen", "angar", "vall", "sand", "sanden"),
    "no": ("eng", "enga", "engene", "voll", "vollen", "oy", "oya", "sand", "sanden"),
    "da": ("eng", "enge"),
    "fr": ("ile", "iles", "isle", "ilot"),
}
# (d) 相対の修飾語（下の側, 上の側）
QUALIFIERS: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {
    "en": (("low", "lower", "nether"), ("high", "upper")),
    "de": (("unter", "nieder"), ("ober",)),
    "nl": (("laag", "neder", "beneden"), ("hoog", "boven")),
    "fr": (("bas", "basse", "sous"), ("haut", "haute", "hauts")),
    "it": (("basso", "bassa", "bassi", "sotto", "inferiore"),
           ("alto", "alta", "alti", "sopra", "superiore")),
    "es": (("bajo", "baja", "bajos", "bajas", "abajo"),
           ("alto", "alta", "altos", "altas", "arriba")),
    "pt": (("baixa", "baixo"), ("alto", "cima")),
    "el": (("kato",), ("ano", "epano")),
    "sh": (("donji", "donja", "donje"), ("gornji", "gornja", "gornje")),
    "sl": (("spodnji", "spodnja", "spodnje"),
           ("gorenji", "gorenja", "gorenje", "zgornji", "zgornja", "zgornje", "gornji")),
    "bg": (("dolna", "dolni", "dolno", "dolen"), ("gorna", "gorni", "gorno", "goren")),
    "sq": (("poshte", "poshtme"), ("siper", "siperme")),
    "ro": (("jos",), ("sus",)),
    "hu": (("also",), ("felso",)),
    "cs": (("dolni", "dolny", "dolne", "dolna"), ("horni", "horny", "horne", "horna")),
    "pl": (("dolny", "dolna", "dolne"), ("gorny", "gorna", "gorne")),
    "uk": (("nyzhnii", "nyzhnia", "nyzhnie", "nizhniy", "nizhnyaya", "nizhneye", "nizhne",
            "nizhni"),
           ("verkhnii", "verkhnia", "verkhnie", "verkhniy", "verkhnyaya", "verkhneye",
            "verkhne")),
    "lt": (("zemutinis", "zemutine", "zemieji"), ("aukstutinis", "aukstoji")),
    "et": (("alam",), ("ulem",)),
    "fi": (("ala",), ("yla",)),
    "sv": (("nedre",), ("ovre",)),
    "no": (("nedre",), ("ovre",)),
    "da": (("nedre",), ("ovre",)),
    "tr": (("asagi",), ("yukari",)),
    "ar": (("tahta", "tahtani", "tahtaniya", "sufla", "sufli"),
           ("fawqani", "foqani", "fouqani", "fawqa", "oulia", "ulya", "ulia")),
    "fa": (("nichla", "nichli", "nichle", "payan", "payin", "paeen", "pain", "payeen"),
           ("bala",)),
    "hi": (("nichla",), ("upar",)),
    "id": (("bawah", "hilir", "ilir"), ("atas", "hulu", "ulu")),
    "tl": (("ibaba", "baba"), ("itaas",)),
    "sw": (("chini",), ("juu",)),
    "mg": (("ambany",), ("ambony",)),
    "ka": (("kvemo",), ("zemo",)),
}
# (e) 泉・井戸
SPRINGS: dict[str, tuple[str, ...]] = {
    "en": ("well", "wells", "spring", "springs", "soak"),
    "de": ("born", "brunn", "brunnen", "quell", "quelle"),
    "fr": ("fontaine", "fontaines", "source", "sources", "puits"),
    "it": ("fonte", "fontana", "fontane", "sorgente", "pozzo", "pozzi"),
    "es": ("fuente", "fuentes", "pozo", "pozos", "manantial", "ojo"),
    "pt": ("fonte", "fontes", "poco", "nascente", "cacimba"),
    "el": ("pigi", "pigadi", "pigadia", "vrysi", "vrisi", "vryses", "vrises", "krini",
           "kefalovryso", "vrys", "vris"),
    "sh": ("vrelo", "vrela", "izvor", "izvori", "bunar", "studenac", "vrutak"),
    "sl": ("studenec", "izvir", "studenci"),
    "bg": ("izvor", "izvori", "kladenets", "kladenci", "bunar", "cheshma"),
    "sq": ("burim", "burimi", "pus", "puse"),
    "ro": ("izvor", "izvorul", "izvoare", "fantana", "fantanele"),
    "hu": ("kut", "forras"),
    "cs": ("studnice", "studna", "studnicka", "pramen", "prameny"),
    "pl": ("zrodlo", "zrodla", "studnia", "studzianka"),
    "uk": ("krynytsia", "krinitsa", "kolodez", "kolodiazi", "dzherelo", "istok"),
    "lt": ("saltinis", "saltiniai"),
    "lv": ("avots", "avoti"),
    "et": ("allikas",),
    "fi": ("lahde",),
    "sv": ("kalla", "kallan"),
    "no": ("kilde",),
    "da": ("kilde",),
    "tr": ("pinar", "kuyu", "bulak", "kaynak"),
    "ar": ("ain", "ayn", "aioun", "ayoun", "bir", "biar", "abar", "hassi"),
    "fa": ("chashma", "chashmeh", "chah", "chahi"),
    "vi": ("gieng",),
    "id": ("sumur",),
    "tl": ("bukal",),
    "sw": ("chemchemi", "kisima"),
    "bantu": ("chitsime",),
    "af": ("fontein", "put"),
    "mg": ("loharano",),
    "qu": ("puquio", "pukyu"),
    "nah": ("ameyal",),
    "ka": ("tskaro", "tskarosi"),
    "kk": ("bulak", "bulaq"),
    "so": ("ceel",),
    "am": ("minch", "minchi"),
}
# (f) 早瀬・滝
RAPIDS: dict[str, tuple[str, ...]] = {
    "en": ("falls",), "fi": ("koski",), "sv": ("fors", "forsen"), "no": ("foss", "fossen"),
}
# (g) 名指された欠け: (言語, 形, 位置, 区分)
ADDITIONS: tuple[tuple[str, tuple[str, ...], tuple[str, ...], str], ...] = (
    ("fr", ("sur",), ("word",), "riverine"),
    ("cs", ("nad",), ("word",), "riverine"),
    ("pl", ("nad",), ("word",), "riverine"),
    ("de", ("hamm", "damm", "deich"), ("word", "suffix"), "riverine"),
    ("es", ("rambla",), ("word",), "water"),
    ("hu", ("fok",), ("word", "suffix"), "water"),
    ("km", ("trapeang",), ("word",), "water"),
)
# (h) 国から外すブロック
BLOCK_REMOVALS: dict[str, tuple[str, ...]] = {"BE": ("de",), "IT": ("de",)}


def _removals() -> dict[tuple[str, str], tuple[str, tuple[str, ...] | None]]:
    """(言語, 形) → (規則の名, 外す位置。None は全て)。"""
    out: dict[tuple[str, str], tuple[str, tuple[str, ...] | None]] = {}
    for lang, form, pos in HOMOGRAPHS:
        out[(lang, form)] = ("b_homograph", pos)
    tables = (("c_riverine_named", RIVERINE_NAMED), ("e_spring_well", SPRINGS),
              ("f_rapids", RAPIDS))
    for rule, table in tables:
        for lang, forms in table.items():
            for f in forms:
                out[(lang, f)] = (rule, None)
    for lang, (lo, hi) in QUALIFIERS.items():
        for f in (*lo, *hi):
            out[(lang, f)] = ("d_qualifier", None)
    return out


def _drop_forms(out: dict, rem: dict) -> tuple[list[dict], set[tuple[str, str]]]:
    """(b)〜(f): 各項目から外す形を除く（位置を限った形はその位置だけ）。"""
    hit: set[tuple[str, str]] = set()
    dropped: list[dict] = []
    for lang, spec in (out.get("languages") or {}).items():
        new = []
        for e in spec.get("entries") or []:
            keep = []
            for f in e["form"]:
                k = (lang, str(f).lower())
                r = rem.get(k)
                if r is not None and (r[1] is None or e["pos"] in r[1]):
                    hit.add(k)
                    dropped.append({"lang": lang, "form": k[1], "pos": e["pos"],
                                    "cat": e["cat"], "rule": r[0]})
                else:
                    keep.append(f)
            if keep:
                new.append({"form": keep, "pos": e["pos"], "cat": e["cat"]})
        spec["entries"] = new
    return dropped, hit


def _add_forms(out: dict) -> list[dict]:
    """(g): 名指された欠けを足す（言語のブロックが無ければ飛ばす）。"""
    added = []
    for lang, forms, poss, cat in ADDITIONS:
        spec = out["languages"].get(lang)
        if spec is None:
            continue
        for pos in poss:
            spec["entries"].append({"form": list(forms), "pos": pos, "cat": cat})
            added.extend({"lang": lang, "form": f, "pos": pos, "cat": cat, "rule": "g_addition"}
                         for f in forms)
    return added


def _remove_blocks(out: dict) -> dict:
    """(h): 国のブロックを外す。"""
    changed = {}
    for cc, langs in BLOCK_REMOVALS.items():
        if cc not in out.get("countries", {}):
            continue
        before = list(out["countries"][cc])
        out["countries"][cc] = [x for x in before if x not in langs]
        changed[cc] = {"before": before, "after": out["countries"][cc]}
    return changed


def correct_lexicon_doc(doc: dict, *, strict: bool = True) -> tuple[dict, dict]:
    """C2: 凍結した辞書の文書から、誤りを直した辞書の文書と記録を作る（元は変えない）。

    (a) v57 の掃除の後に (b)〜(f) の形を外し、(g) を足し、(h) で国のブロックを外す。
    strict なら、外すと決めた (言語, 形) が辞書に無いときに止める（書き間違いを見逃さない）。
    """
    out, dropped_a = clean_lexicon_doc(doc)
    rem = _removals()
    dropped, hit = _drop_forms(out, rem)
    missing = sorted(k for k in rem if k not in hit)
    if strict and missing:
        msg = f"C2 removals not found in the lexicon: {missing}"
        raise ValueError(msg)
    first = [{"lang": a, "form": b, "pos": c, "rule": "a_v57_clean"} for a, b, c in dropped_a]
    return out, {"dropped": first + dropped, "added": _add_forms(out),
                 "countries": _remove_blocks(out)}


def qualifier_side(lang_forms: dict[str, tuple[tuple[str, ...], tuple[str, ...]]],
                   langs: tuple[str, ...]) -> dict[str, str]:
    """国のブロックの修飾語の形 → 'lower' / 'upper'（同じ形が両側にあれば使わない）。"""
    side: dict[str, set[str]] = {}
    for lang in langs:
        lo, hi = lang_forms.get(lang, ((), ()))
        for f in lo:
            side.setdefault(f, set()).add("lower")
        for f in hi:
            side.setdefault(f, set()).add("upper")
    return {f: next(iter(s)) for f, s in side.items() if len(s) == 1}


# %% B15: 相対の修飾語の対
_SPLIT = re.compile(r"[^a-z]+")
PREFIX_QUALIFIERS = {"de": ("unter", "nieder", "ober"),
                     "nl": ("laag", "neder", "beneden", "hoog", "boven"),
                     "hu": ("also", "felso"), "en": ("nether",)}
MIN_BASE = 3
MIN_REST = 2


def split_qualifier(name: str, side_of: dict[str, str],
                    prefix_forms: tuple[str, ...]) -> tuple[str, str] | None:
    """名から修飾語を一つだけ外し（側, 基の名）を返す。

    修飾語が無い・二つ以上・基が 3 文字未満なら None。語の位置: トークンが修飾語に一致。
    頭の位置（prefix_forms の形だけ）: トークンがその形で始まり残りが 2 文字以上
    （例: Oberhausen → ober + hausen）。
    """
    toks = [t for t in _SPLIT.split(str(name).lower()) if t]
    found: list[tuple[int, str, str]] = []
    for i, t in enumerate(toks):
        if t in side_of:
            found.append((i, side_of[t], ""))
            continue
        for f in prefix_forms:
            if f in side_of and t.startswith(f) and len(t) - len(f) >= MIN_REST:
                found.append((i, side_of[f], t[len(f):]))
                break
    if len(found) != 1:
        return None
    i, side, rest = found[0]
    base = " ".join([*toks[:i], *([rest] if rest else []), *toks[i + 1:]])
    if len(base.replace(" ", "")) < MIN_BASE:
        return None
    return side, base


def match_pairs(lower_xy: np.ndarray, upper_xy: np.ndarray,
                max_m: float) -> list[tuple[int, int, float]]:
    """下と上の点の対（互いに max_m 以内）を近い順に貪欲にとる（各点は一回だけ）。"""
    if not len(lower_xy) or not len(upper_xy):
        return []
    dx = lower_xy[:, None, 0] - upper_xy[None, :, 0]
    dy = lower_xy[:, None, 1] - upper_xy[None, :, 1]
    d = np.hypot(dx, dy)
    cand = sorted((float(d[i, j]), i, j) for i, j in zip(*np.nonzero(d <= max_m), strict=True))
    used_l: set[int] = set()
    used_u: set[int] = set()
    out = []
    for dist, i, j in cand:
        if i in used_l or j in used_u:
            continue
        used_l.add(i)
        used_u.add(j)
        out.append((int(i), int(j), dist))
    return out


# %% C1(b): 事象の前後
def _utc(t: str) -> datetime:
    d = datetime.fromisoformat(t)
    return d.replace(tzinfo=UTC) if d.tzinfo is None else d.astimezone(UTC)


def event_before(t: str, t0: str, min_days: float) -> bool:
    """事象の日付 t が t0 より min_days 日を超えて前か（読めなければ False）。"""
    try:
        return (_utc(t0) - _utc(t)).total_seconds() / 86400.0 > min_days
    except ValueError:
        return False


# %% C3: 場所グループを分けた早期終了
ES_VAL_SHARE = 0.1
ES_MAX_ITER = 200
ES_LEARNING_RATE = 0.1


def grouped_n_iter(x: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
                   seed: int, max_iter: int = ES_MAX_ITER) -> tuple[int, dict]:
    """訓練側の場所グループの 1 割を内側の検証にし、検証の log loss が最小の反復の回数。

    検証のグループは seed で選ぶ。検証に浸水と乾きの両方が無ければ max_iter を返す。
    """
    rng = np.random.default_rng(seed)
    g = np.asarray(groups).astype(str)
    u = np.unique(g)
    k = max(1, round(len(u) * ES_VAL_SHARE))
    val = np.isin(g, rng.choice(u, size=k, replace=False))
    yv = np.asarray(y)[val]
    if val.all() or not val.any() or len(np.unique(yv)) < 2:  # noqa: PLR2004
        return max_iter, {"n_val_groups": int(k), "reason": "degenerate validation"}
    m = HistGradientBoostingClassifier(random_state=seed, max_iter=max_iter,
                                       learning_rate=ES_LEARNING_RATE, early_stopping=False)
    m.fit(x[~val], np.asarray(y)[~val])
    eps = 1e-12
    losses = []
    for p in m.staged_predict_proba(x[val]):
        q = np.clip(p[:, 1], eps, 1 - eps)
        losses.append(float(-np.mean(yv * np.log(q) + (1 - yv) * np.log(1 - q))))
    best = int(np.argmin(losses)) + 1
    return best, {"n_val_groups": int(k), "n_val_cells": int(val.sum()), "best": best,
                  "loss_best": losses[best - 1], "loss_last": losses[-1]}


def fit_grouped(x_tr: np.ndarray, y_tr: np.ndarray, groups: np.ndarray, x_te: np.ndarray, *,
                seed: int) -> tuple[np.ndarray, dict]:
    """C3: 反復の回数を grouped_n_iter で選び、訓練側の全体で当て直した予測。"""
    n, info = grouped_n_iter(x_tr, y_tr, groups, seed=seed)
    m = HistGradientBoostingClassifier(random_state=seed, max_iter=n,
                                       learning_rate=ES_LEARNING_RATE, early_stopping=False)
    m.fit(x_tr, y_tr)
    return m.predict_proba(x_te)[:, 1], info


# %% B1: Fieller 型の同等性の検定
def l_pm(delta_bs: np.ndarray | float, delta_base: np.ndarray | float,
         margin: float = SESOI_R) -> tuple[np.ndarray, np.ndarray]:
    """L+ = ΔBS − margin·Δ_base と L− = ΔBS + margin·Δ_base。

    ΔBS = BS(M1c) − BS(M3n)、Δ_base = BS(M0) − BS(M1)。
    """
    e = np.asarray(delta_bs, dtype=np.float64)
    d = np.asarray(delta_base, dtype=np.float64)
    return e - margin * d, e + margin * d


def tost_from_replicates(l_plus: np.ndarray, l_minus: np.ndarray,
                         alpha: float = 0.05) -> dict:
    """再標本の L± から片側の限界と同等の判定（L+ の上限 < 0 かつ L− の下限 > 0）。"""
    lp = np.asarray(l_plus, dtype=np.float64)
    lm = np.asarray(l_minus, dtype=np.float64)
    ok = np.isfinite(lp) & np.isfinite(lm)
    up = float(np.quantile(lp[ok], 1 - alpha))
    lo = float(np.quantile(lm[ok], alpha))
    return {"l_plus_upper": up, "l_minus_lower": lo, "equivalent": bool(up < 0 < lo),
            "reps_ok": int(ok.sum())}


Z90 = 1.6448536269514722


def fieller_interval(e: float, d: float, var: tuple[float, float, float],
                     z: float = Z90) -> tuple[float, float] | None:
    """R = e/d の Fieller の区間（正規近似、既定は 90%）。var = (Var e, Var d, Cov(e, d))。

    d が 0 と区別できなければ None（区間は有界でない）。
    """
    var_e, var_d, cov_ed = var
    a = d * d - z * z * var_d
    if a <= 0:
        return None
    b = e * d - z * z * cov_ed
    c = e * e - z * z * var_e
    disc = b * b - a * c
    if disc < 0:
        return None
    r = float(np.sqrt(disc))
    return (b - r) / a, (b + r) / a


def nested_weights(rng: np.random.Generator, country_of_act: np.ndarray) -> np.ndarray:
    """国を引き、引いた国ごとにその国の activation を引き直す二段の再抽出の重み。

    引くのは activation を一つ以上持つ国だけ（その数だけ等確率で。B-13 の追記 6）。
    符号の表の国には、他の国に入れた activation の AOI にセルがあるだけの国があり、
    それを引くと引き直す activation が無い。
    """
    countries = np.unique(country_of_act)
    k = len(countries)
    wc = rng.multinomial(k, np.full(k, 1.0 / k))
    w = np.zeros(len(country_of_act), dtype=np.float64)
    for c, times in zip(countries, wc, strict=True):
        acts = np.nonzero(country_of_act == c)[0]
        for _ in range(int(times)):
            w[acts] += rng.multinomial(len(acts), np.full(len(acts), 1.0 / len(acts)))
    return w


def nadeau_bengio_var(values: np.ndarray, test_train_ratio: float = 1 / 9) -> float:
    """Nadeau–Bengio の補正した分散 (1/J + n_test/n_train) S²（fold の割り当てを変えた値）。"""
    v = np.asarray(values, dtype=np.float64)
    return float((1 / len(v) + test_train_ratio) * np.var(v, ddof=1))
