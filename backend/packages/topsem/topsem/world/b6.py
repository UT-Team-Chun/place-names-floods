"""事後 B-6（初見の査読 第 2 回への対応）の部品。

計算と規則は amendment B-6 で計算の前に固定したもの。

- 2: 実データのラベルに名前の効果を注入する模擬（plasmode）の注入と、読み違いを入れた潜在の印
- 3: 浸水の情報を集落ごとに持たせた対照（事象のラベルの割合、leave-one-out）
- 5: 厳密な読み（「水」の区分から泉・井戸・水場の形を除き、水域・湿地・河川だけ）
- 6(a): 地名辞典の地物コード（PPLX・PPLH・PPLQ・PPLW）を除き、近くの同名の重複を 1 つに
"""

from __future__ import annotations

import copy
import zlib
from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from topsem.world.lexicon import Lexicon

# B-6 付表: 泉・井戸・水場とみなす形（言語の符号 → 形）。語の意味だけから決め、名前に当てずに固定
_SPRING = {
    "en": "well wells spring springs soak waterhole",
    "de": "born brunn brunnen quell quelle",
    "fr": "fontaine fontaines source sources puits",
    "it": "fontana fontane fonte pozzo pozzi sorgente",
    "es": "fuente fuentes manantial pozo pozos ojo",
    "pt": "fonte fontes nascente poco cacimba",
    "el": "krini pigi pigadi pigadia vris vrises vrisi vrys vryses vrysi kefalovryso",
    "sh": "bunar izvor izvori studenac vrelo vrela vrutak",
    "sl": "izvir studenec studenci",
    "bg": "bunar cheshma izvor izvori kladenets kladenci",
    "sq": "burim burimi pus puse",
    "ro": "fantana fantanele izvor izvorul izvoare",
    "hu": "forras kut",
    "cs": "pramen prameny studna studnice studnicka",
    "pl": "studnia studzianka zrodla zrodlo",
    "uk": "dzherelo istok kolodez kolodiazi krinitsa krynytsia",
    "lt": "saltinis saltiniai",
    "lv": "avots avoti",
    "et": "allikas",
    "fi": "lahde",
    "sv": "kalla kallan",
    "no": "kilde",
    "da": "kilde",
    "tr": "bulak kaynak kuyu pinar",
    "ar": "ain ayn aioun ayoun bir biar abar hassi",
    "fa": "chah chahi chashma chashmeh",
    "vi": "gieng",
    "id": "sumur",
    "tl": "bukal",
    "sw": "chemchemi kisima",
    "bantu": "chitsime",
    "af": "fontein put",
    "mg": "loharano",
    "qu": "pukyu puquio",
    "nah": "ameyal",
    "ka": "tskaro tskarosi",
    "kk": "bulak bulaq",
    "so": "ceel",
    "am": "minch minchi",
}
SPRING_FORMS: dict[str, frozenset[str]] = {k: frozenset(v.split()) for k, v in _SPRING.items()}
STRICT_CATEGORIES = ("water", "wet", "riverine")  # 厳密な読み（水は泉・井戸を除いた水域）
GN_EXCLUDE_FCODES = ("PPLX", "PPLH", "PPLQ", "PPLW")
DEDUP_RADIUS_M = 100.0
EARTH_R = 6_371_008.8


# %% 5: 厳密な読み
def strict_lexicon(lex: Lexicon) -> Lexicon:
    """「水」の区分から泉・井戸・水場の形（SPRING_FORMS）を除いた辞書の写し（他の区分はそのまま）。"""
    out = copy.copy(lex)
    out.by_lang = {
        lang: tuple(e for e in ents
                    if not (e.cat == "water" and e.form in SPRING_FORMS.get(lang, frozenset())))
        for lang, ents in lex.by_lang.items()}
    out.sha256 = f"{lex.sha256}+B-6-strict"
    return out


# %% 2: plasmode
def _splitmix64(x: np.ndarray) -> np.ndarray:
    z = x.astype(np.uint64) + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


def hashed_uniform(keys: np.ndarray, salt: int) -> np.ndarray:
    """鍵（uint64）と塩から決まる [0, 1) の一様乱数（fold をまたいで同じセル・集落に同じ値）。"""
    with np.errstate(over="ignore"):
        s = _splitmix64(np.array([salt], dtype=np.uint64))
        h = _splitmix64(keys.astype(np.uint64) ^ s)
    return (h >> np.uint64(11)).astype(np.float64) / float(1 << 53)


def cell_keys(aoi: Iterable[str], ix: np.ndarray, iy: np.ndarray) -> np.ndarray:
    """セルの鍵: AOI の名の CRC32 と格子の位置。"""
    a = np.array([zlib.crc32(str(s).encode("utf-8")) for s in aoi], dtype=np.uint64)
    return (a << np.uint64(32)) | (ix.astype(np.uint64) << np.uint64(16)) | iy.astype(np.uint64)


def plasmode_salt(form: str, delta: float, rep: int) -> int:
    return zlib.crc32(f"B-6/plasmode/{form}/{delta:.4f}/{rep}".encode())


def inject_labels(y: np.ndarray, p: np.ndarray, target: np.ndarray, delta: float | np.ndarray,
                  u: np.ndarray) -> np.ndarray:
    """対象のセルのラベルに、地形を条件にしたオッズを e^δ 倍にする効果を入れる。

    δ > 0 ならラベル 0 のセルを確率 q = (p′ − p)/(1 − p) で 1 に、δ < 0（B-9 の C2）なら
    ラベル 1 のセルを確率 1 − p′/p で 0 にする。p′ = p e^δ / (1 − p + p e^δ)。
    y が Bernoulli(p) なら、注入後は Bernoulli(p′)。δ はセルごとの配列でもよい
    （B-9 の C3: 集落ごとに違う効果）。δ > 0 は B-6 の実装と同じ乱数・同じ式で、
    同じラベルを返す。
    """
    p = np.clip(p.astype(np.float64), 1e-6, 1 - 1e-6)
    d = np.asarray(delta, dtype=np.float64)
    e = np.exp(d)
    p2 = p * e / (1 - p + p * e)
    q = (p2 - p) / (1 - p)
    up = target & (y == 0) & (d > 0) & (u < q)
    down = target & (y == 1) & (d < 0) & (u < 1 - p2 / p)
    return np.where(up, 1, np.where(down, 0, y)).astype(np.int8)


def latent_names(gids: np.ndarray, wetside: np.ndarray, pi1: float, pi0: float,
                 salt: int) -> np.ndarray:
    """読み違いを入れた潜在の印 L（集落ごと）。

    辞書が水系側と読んだ集落は確率 π₁、それ以外は π₀ で L = 1。集落の gid と塩で決まる。
    """
    u = hashed_uniform(gids.astype(np.int64).astype(np.uint64), salt)
    return np.where(wetside, u < pi1, u < pi0) & (gids >= 0)


# %% 3: 事象のラベルの割合（漏れを承知の上限）
def label_share_by_settlement(cells_xy: np.ndarray, y: np.ndarray, settle_xy: np.ndarray, *,
                              settle_gid: np.ndarray, radius_m: float) -> pd.DataFrame:
    """集落の点から radius_m 以内の評価可能なセル（y ∈ {0, 1}）の浸水の和と数。"""
    tree = cKDTree(cells_xy)
    hits = tree.query_ball_point(settle_xy, radius_m)
    s = np.array([y[h].sum() if len(h) else 0 for h in hits], dtype=np.float64)
    n = np.array([len(h) for h in hits], dtype=np.float64)
    return pd.DataFrame({"gid": settle_gid.astype(np.int64), "sum": s, "n": n})


def loo_share(cell_y: np.ndarray, cell_xy: np.ndarray, gid: np.ndarray,  # noqa: PLR0913
              tab: pd.DataFrame, *,
              settle_xy: dict[int, tuple[float, float]], radius_m: float) -> np.ndarray:
    """最寄りの集落の割合。セル自身がその集落から radius_m 以内なら自分を除く（leave-one-out）。"""
    t = tab.drop_duplicates("gid").set_index("gid")
    s = t["sum"].reindex(gid).to_numpy(dtype=np.float64)
    n = t["n"].reindex(gid).to_numpy(dtype=np.float64)
    pts = np.array([settle_xy.get(int(g), (np.nan, np.nan)) for g in gid], dtype=np.float64)
    inside = np.hypot(cell_xy[:, 0] - pts[:, 0], cell_xy[:, 1] - pts[:, 1]) <= radius_m
    s = s - np.where(inside, cell_y, 0)
    n = n - inside
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, s / n, np.nan)


# %% 6(a): 地名辞典の絞り込み
def _unit(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lo, la = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def filter_gazetteer(gn: pd.DataFrame, exclude_fcodes: tuple[str, ...] = GN_EXCLUDE_FCODES,
                     dedup_m: float = DEDUP_RADIUS_M) -> tuple[pd.DataFrame, dict[str, int]]:
    """P 類の集落から地物コードを除き、同じ asciiname で dedup_m 以内の重複は gid の小さい方に。"""
    p = gn.loc[gn["fclass"] == "P"]
    keep = ~p["fcode"].isin(exclude_fcodes)
    q = p.loc[keep].reset_index(drop=True)
    tree = cKDTree(_unit(q["lon"].to_numpy(float), q["lat"].to_numpy(float)))
    pairs = tree.query_pairs(dedup_m / EARTH_R, output_type="ndarray")
    names = q["asciiname"].astype(str).str.lower().to_numpy()
    gids = q["gid"].to_numpy(np.int64)
    drop = {int(i) if gids[i] > gids[j] else int(j) for i, j in pairs if names[i] == names[j]}
    out = q.drop(index=sorted(drop)).reset_index(drop=True)
    rest = gn.loc[gn["fclass"] != "P"]
    info = {"n_p": len(p), "n_excluded_fcode": int((~keep).sum()), "n_duplicates": len(drop),
            "n_kept": len(out)}
    return pd.concat([out, rest], ignore_index=True), info
