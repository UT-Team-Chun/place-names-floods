"""事後 B-6 の部品の単体テスト（合成データのみ）。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from topsem.config import REPO_ROOT
from topsem.world.b6 import (
    SPRING_FORMS,
    cell_keys,
    filter_gazetteer,
    hashed_uniform,
    inject_labels,
    label_share_by_settlement,
    latent_names,
    loo_share,
    plasmode_salt,
    strict_lexicon,
)
from topsem.world.lexicon import Lexicon

LEXICON = REPO_ROOT / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"


def test_hashed_uniform_is_deterministic_and_uniform() -> None:
    keys = np.arange(200_000, dtype=np.uint64)
    u1, u2 = hashed_uniform(keys, 7), hashed_uniform(keys, 7)
    assert np.array_equal(u1, u2)
    assert not np.array_equal(u1, hashed_uniform(keys, 8))
    assert u1.min() >= 0.0
    assert u1.max() < 1.0
    assert abs(u1.mean() - 0.5) < 0.01


def test_cell_keys_distinguish_aoi_and_position() -> None:
    k = cell_keys(["A", "A", "B"], np.array([1, 2, 1]), np.array([5, 5, 5]))
    assert len(set(k.tolist())) == 3


def test_inject_labels_multiplies_odds_by_exp_delta() -> None:
    """y ~ Bernoulli(p) に注入すると、対象のセルの浸水の割合が p′ = p e^δ/(1 − p + p e^δ) に。"""
    rng = np.random.default_rng(0)
    n = 400_000
    p = np.full(n, 0.1)
    y = (rng.random(n) < p).astype(np.int8)
    target = np.ones(n, dtype=bool)
    target[n // 2:] = False
    delta = np.log(3.0)
    y2 = inject_labels(y, p, target, delta, hashed_uniform(np.arange(n, dtype=np.uint64), 1))
    p2 = 0.1 * 3 / (1 - 0.1 + 0.1 * 3)
    assert abs(y2[target].mean() - p2) < 0.005
    assert np.array_equal(y2[~target], y[~target])  # 対象外は変えない
    assert (y2 >= y).all()  # 1 を 0 にはしない


def test_latent_names_rates() -> None:
    g = np.arange(100_000, dtype=np.int64)
    wet = g % 2 == 0
    lat = latent_names(g, wet, 0.76, 0.18, plasmode_salt("CA", 1.0, 0))
    assert abs(lat[wet].mean() - 0.76) < 0.01
    assert abs(lat[~wet].mean() - 0.18) < 0.01
    assert not latent_names(np.array([-1]), np.array([True]), 1.0, 1.0, 0)[0]


def test_loo_share_excludes_the_cell_itself() -> None:
    """集落から 500 m 以内の 3 セル（浸水 1, 0, 1）。内側は自分を除き、外側は 3 セルの割合。"""
    cells = np.array([[0.0, 0.0], [100.0, 0.0], [200.0, 0.0], [5000.0, 0.0]])
    y = np.array([1.0, 0.0, 1.0, 1.0])
    tab = label_share_by_settlement(cells, y, np.array([[0.0, 0.0]]), settle_gid=np.array([9]),
                                    radius_m=500.0)
    assert tab["sum"].tolist() == [2.0]
    assert tab["n"].tolist() == [3.0]
    share = loo_share(y, cells, np.array([9, 9, 9, 9]), tab, settle_xy={9: (0.0, 0.0)},
                      radius_m=500.0)
    assert share.tolist() == [0.5, 1.0, 0.5, 2 / 3]


def test_strict_lexicon_drops_spring_forms_only() -> None:
    lex = Lexicon(LEXICON)
    ls = strict_lexicon(lex)
    assert lex.match("Holywell", "GB")["water"] == 1
    assert ls.match("Holywell", "GB")["water"] == 0  # 泉・井戸は厳密な読みでは水域でない
    assert ls.match("Riverside", "GB")["water"] == lex.match("Riverside", "GB")["water"]
    assert ls.match("Marshfield", "GB")["wet"] == lex.match("Marshfield", "GB")["wet"]
    # 付表の形は、その言語の水の区分の形に実際にある（綴りの誤りが無い）
    water_forms = {lang: {e.form for e in ents if e.cat == "water"}
                   for lang, ents in lex.by_lang.items()}
    for lang, forms in SPRING_FORMS.items():
        assert forms <= water_forms[lang], (lang, forms - water_forms[lang])


def test_filter_gazetteer_drops_fcodes_and_near_duplicates() -> None:
    gn = pd.DataFrame({
        "gid": [1, 2, 3, 4, 5],
        "asciiname": ["Alpha", "Alpha", "Beta", "Gamma", "Alpha"],
        "lat": [50.0, 50.0003, 50.0, 50.1, 51.0],  # 2 は 1 から約 33 m、5 は遠い
        "lon": [10.0, 10.0, 10.0, 10.1, 10.0],
        "fclass": ["P", "P", "P", "P", "P"],
        "fcode": ["PPL", "PPL", "PPLX", "PPL", "PPL"],
    })
    out, info = filter_gazetteer(gn)
    assert sorted(out["gid"].tolist()) == [1, 4, 5]
    assert info == {"n_p": 5, "n_excluded_fcode": 1, "n_duplicates": 1, "n_kept": 3}


def test_spring_table_matches_amendment() -> None:
    """付表の形が amendment B-6 の表と一致する（コードと文書の食い違いを防ぐ）。"""
    doc = (Path(REPO_ROOT) / "docs" / "protocol" / "amendments" / "B-6_blind_review2_analyses.md")
    rows = [line for line in doc.read_text(encoding="utf-8").splitlines()
            if line.startswith("| ") and not line.startswith(("| 符号", "|---"))]
    table = {r.split("|")[1].strip(): frozenset(f.strip() for f in r.split("|")[2].split(","))
             for r in rows}
    assert table == SPRING_FORMS
