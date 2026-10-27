# %% NHESS 版（paper_2_nhess）の図表・数値マクロ。数値は全て成果物由来（手打ち禁止）。
# 旧版（paper_2_global/figures/make_figures.py）の関数をライブラリとして読み込み、同じ入力から同じ数値を作る。
# その上に事後の解析（amendment B-3〜B-5: v41・v42・v43・v44・v45・v46・v47）の数値・図・付録の表を足す。
# 使い方（docs/paper/paper_2_nhess/ から）: make figures
#
# 入力が無い数値は \missing{...} の印にしてコンパイルは通す（投稿前に `make check` で印が無いことを確かめる）。
from __future__ import annotations

import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "0")  # PDF に作成日時を埋めない
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
OLD = REPO / "docs" / "paper" / "paper_2_global" / "figures" / "make_figures.py"
VER = REPO / "docs" / "verification"
DATA = REPO / "data"
SESOI = 0.05


def _load_old():
    spec = importlib.util.spec_from_file_location("p2g_figures", OLD)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["p2g_figures"] = mod
    spec.loader.exec_module(mod)
    mod.OUTDIR = HERE  # 旧版の図・表もこのフォルダに書く
    return mod


g = _load_old()
missing, rv, sci, pct, _n = g.missing, g.rv, g.sci, g.pct, g._n


def _json(name: str) -> dict | None:
    return g._json(VER / name)


# %% 旧版と同じ数値（確証解析・B-2 の事後解析）。旧版の main と同じ順に呼び、本文の図 1–2・梯子の図と表を作る
def base_macros() -> tuple[dict[str, str], dict | None, dict, dict]:
    m: dict[str, str] = {}
    m |= g.inventory_macros(g._json(REPO / "docs" / "data" / "ems_inventory_20261025.json"))
    inc = g._json(DATA / "global" / "incidence.json")
    if inc:
        m |= g.incidence_macros(inc) | g.sample_macros(inc)
        g.fig_map(inc)
    m |= g.lexicon_macros(REPO / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")
    g.fig_design(m["NLexLangs"])
    m |= g.design_macros() | g.channel_macros()
    m |= g.freeze_macros(REPO / "docs" / "protocol" / "preregistration_v2.md")
    m |= g.geonames_macros()
    m |= g.nread_macros(g._json(DATA / "global" / "nread" / "nread_meta.json"))
    m |= g.power_macros(g._json(DATA / "runs" / "power_v3_null" / "power_summary.json"),
                        g._json(DATA / "runs" / "power_v3_sesoi" / "power_summary.json"))
    run = DATA / "runs" / "world_v2"
    primary = g._json(run / "metrics_primary.json")
    m |= g.result_macros(primary)
    extra = {s: v for s in ("country", "temporal", "llmcut") if (v := g._json(run / f"metrics_{s}.json"))}
    post = {k: _json(f) for k, f in (("v36", "v36_yardstick_robustness.json"),
                                     ("v37", "v37_positive_control.json"),
                                     ("v38", "v38_ems_sensors.json"),
                                     ("v39", "v39_within_aoi_odds.json"))}
    if primary:
        m |= g.detail_macros(primary, post["v36"])
        m |= g.split_macros(extra)
        if primary.get("ladder"):
            g.fig_ladder(primary["ladder"], primary["primary"]["delta_bs_M0_M1"])
        g.write_si_tables(primary, extra, post)
        m |= g.robustness_macros(post["v36"])
        m |= g.positive_control_macros(post["v37"])
        m |= g.sensor_macros(post["v38"])
        m |= g.odds_macros(post["v39"])
        m |= g.power_split_macros(DATA / "runs")
        m |= g.hgb_macros(_json("v40_hgb_iterations.json"))
        m |= g.sens_macros(DATA / "runs")
        m |= g.contrast_macros(_json("v35_name_class_contrast.json"))
    return m, primary, extra, post


# %% 書式
def orv(x: float) -> str:
    return f"{x:.2f}"


def or_macros(prefix: str, t: dict | None) -> dict[str, str]:
    """オッズ比と 95% 区間（{prefix}, {prefix}Lo, {prefix}Hi）。"""
    if not t or t.get("odds_ratio") is None:
        return {f"{prefix}{s}": missing(prefix) for s in ("", "Lo", "Hi")}
    return {prefix: orv(t["odds_ratio"]), f"{prefix}Lo": orv(t["ci95_odds_ratio"][0]),
            f"{prefix}Hi": orv(t["ci95_odds_ratio"][1])}


def pp(x: float) -> str:
    return f"\\ensuremath{{{100 * x:+.1f}}}".replace("+", "{+}").replace("-", "{-}")


def pct1(x: float) -> str:
    return f"{100 * x:.1f}\\,\\%" if abs(100 * x) < 10 else f"{100 * x:.0f}\\,\\%"


def _word(n: int) -> str:
    return g._word(n)


# %% 凍結の記録（内容を固定したコミットと、凍結を記したコミット）
def freeze_record() -> dict[str, str]:
    text = (REPO / "docs" / "protocol" / "preregistration_v2.md").read_text(encoding="utf-8")
    mt = re.search(r"FROZEN\*\* (\S+) / commit ([0-9a-f]{7,40})", text)
    try:
        mark = subprocess.run(["git", "rev-parse", "--short=7", "prereg-v2-frozen^{commit}"],
                              cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        mark = missing("freeze tag")
    return {"FreezeContentCommit": mt.group(2)[:7] if mt else missing("freeze line"),
            "FreezeMarkCommit": mark}


# %% B-4（v46）
def norway_macros(v: dict | None) -> dict[str, str]:
    if not v:
        return {k: missing("v46") for k in ("NorwayCells", "NorwayMatched", "NorwayShare")}
    return {"NorwayCells": _n(v["n_cells_norway"]), "NorwayMatched": _n(v["n_cells_norway_matched"]),
            "NorwayShare": f"{100 * v['share_matched_of_all']:.2f}\\,\\%",
            "NCountriesNearest": str(v.get("n_countries_nearest", "")),
            "NUncoveredAfterFix": str(len(v.get("uncovered_countries_after_fix", [])))}


# %% B-3（v41・v45、C2 は v47）
def gfm_macros(v41: dict | None, v45: dict | None, v47: dict | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if v41:
        a0 = v41["A0"]["all"]
        out |= {"GFMNAoi": str(a0["n_aoi"]), "GFMNPlaces": str(a0["n_units"]), "GFMNCells": _n(a0["n_cells"]),
                "GFMAgreePct": pct1(a0["agreement"]), "GFMKappa": f"{a0['kappa']:.2f}",
                "GFMRateEms": pct1(a0["rate_ems"]), "GFMRateGfm": pct1(a0["rate_gfm"])}
        for key, name in (("A1_near_minus_far", "GFMAOne"), ("A2_wetside_le500m", "GFMATwo"),
                          ("A2_wetside_all", "GFMATwoAll")):
            r = v41["A1_A2"][key]
            out |= {name: pp(r["estimate"]), f"{name}Lo": pp(r["ci95"][0]), f"{name}Hi": pp(r["ci95"][1])}
        og, oe = v41["odds_gfm"], v41["odds_ems"]
        out |= or_macros("GFMNearFarGfm", og["bands_adjusted_M1"]["terms"]["band_gt1km"])
        out |= or_macros("GFMNearFarEms", oe["bands_adjusted_M1"]["terms"]["band_gt1km"])
        out |= or_macros("GFMORWetNear", og["wetside_adjusted_M1"]["terms"]["wetside:le500m"])
        out |= or_macros("GFMORWetNearRaw", og["wetside_raw"]["terms"]["wetside:le500m"])
        out |= or_macros("GFMORWetNearEms", oe["wetside_adjusted_M1"]["terms"]["wetside:le500m"])
        out |= or_macros("GFMORWetAll", og["wetside_allcells_adjusted_M1"]["terms"]["wetside:all"])
        out |= or_macros("GFMORWetAllEms", oe["wetside_allcells_adjusted_M1"]["terms"]["wetside:all"])
        out["GFMNOddsAoi"] = str(og["n_aoi"])
        oc = v41.get("odds_gfm_c2_M1")  # B の感度（B-3 で登録）: 調整に C2 の GFM で学習した M1
        if oc:
            out |= or_macros("GFMORWetNearCTwo", oc["wetside_adjusted_M1"]["terms"]["wetside:le500m"])
        else:
            out |= or_macros("GFMORWetNearCTwo", None)
        for key, name in (("C1_gfm", "GFMCOne"), ("C1_ems", "GFMCOneEms"), ("C1_gfm_lexfix", "GFMCOneLexfix")):
            c = v41.get(key)
            if c:
                out |= {f"{name}R": rv(c["R"]), f"{name}RLoNinety": rv(c["ci90_R"][0]),
                        f"{name}RHiNinety": rv(c["ci90_R"][1]), f"{name}Verdict": c["decision"]["label"]}
        sens = v41.get("sensitivity_independent_sensor", {}).get("odds_gfm")
        if sens:
            out |= or_macros("GFMORWetNearIndep", sens["wetside_adjusted_M1"]["terms"]["wetside:le500m"])
            out["GFMNIndepAoi"] = str(v41["sensitivity_independent_sensor"]["n_aoi"])
    else:
        out["GFMNAoi"] = missing("v41")
    if v45:
        for cut, name in (("v_ge_0.9", "GFMORWetNearVNine"), ("v_ge_1.0", "GFMORWetNearVOne")):
            out |= or_macros(name, v45[cut]["gfm"]["wetside_adjusted_M1"]["wetside:le500m"])
            out |= or_macros(f"{name}Ems", v45[cut]["ems"]["wetside_adjusted_M1"]["wetside:le500m"])
    if v41 and v41.get("C1_gfm"):
        c1 = v41["C1_gfm"]
        out |= {"GFMCOneYard": sci(c1["yardstick_M0_M1"]), "GFMCOneYardLo": sci(c1["ci95_yardstick"][0]),
                "GFMCOneYardHi": sci(c1["ci95_yardstick"][1])}
    c2 = (v47 or {}).get("gfm_c2")
    det = (v47 or {}).get("gfm_c2_detail")
    keys = ("GFMCTwoR", "GFMCTwoRLoNinety", "GFMCTwoRHiNinety", "GFMCTwoVerdict", "GFMCTwoNames", "GFMCTwoNamesLo",
            "GFMCTwoNamesHi", "GFMCTwoYard", "GFMCTwoYardLo", "GFMCTwoYardHi", "GFMCTwoPlacebo", "GFMCTwoBSMZero",
            "GFMCTwoBSMOne", "GFMCTwoBSMOnes", "GFMCTwoBSMOnec")
    if isinstance(c2, dict) and isinstance(det, dict):
        b = det["brier"]
        out |= {"GFMCTwoR": rv(c2["R"]), "GFMCTwoRLoNinety": rv(c2["ci90_R"][0]),
                "GFMCTwoRHiNinety": rv(c2["ci90_R"][1]), "GFMCTwoVerdict": c2["decision"],
                "GFMCTwoNames": sci(c2["delta_bs"]), "GFMCTwoNamesLo": sci(c2["ci95_delta"][0]),
                "GFMCTwoNamesHi": sci(c2["ci95_delta"][1]), "GFMCTwoYard": sci(c2["yardstick"]),
                "GFMCTwoYardLo": sci(c2["ci95_yardstick"][0]), "GFMCTwoYardHi": sci(c2["ci95_yardstick"][1]),
                "GFMCTwoPlacebo": sci(det["delta_bs_placebo"]), "GFMCTwoBSMZero": f"{b['M0']:.4f}",
                "GFMCTwoBSMOne": f"{b['M1']:.4f}", "GFMCTwoBSMOnes": f"{b['M1s']:.4f}",
                "GFMCTwoBSMOnec": f"{b['M1c']:.4f}"}
    else:
        out |= {k: missing("C2 (world_v2_gfm)") for k in keys}
    return out


# %% B-5（v42 と v47）
EG = (("wetland", "EGWetland"), ("island", "EGIsland"), ("spring_well", "EGSpring"), ("crossing", "EGCrossing"),
      ("lowland", "EGLowland"))


def b5_macros(v42: dict | None, v47: dict | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if v42:
        out |= {"TIORPooledEmp": orv(v42["inversion_area_wide"]["or_upper_q05"]),
                "TIORPooledNorm": orv(v42["inversion_area_wide"]["or_upper_q05_normal"]),
                "TIORLocalEmp": orv(v42["inversion_local_500m"]["or_upper_q05"]),
                "TIORLocalNorm": orv(v42["inversion_local_500m"]["or_upper_q05_normal"])}
        w = v42["within_area"]
        out |= {"AucDiffMedian": f"{w['auc_diff_median']:.4f}",
                "AucDiffQOne": f"\\ensuremath{{{w['auc_diff_iqr'][0]:.4f}}}",
                "AucDiffQThree": f"\\ensuremath{{{w['auc_diff_iqr'][1]:+.4f}}}", "NAucAoi": str(w["n_aoi_auc"]),
                "CondLL": sci(w["cond_ll_per_cell"]), "CondLLLo": sci(w["cond_ll_ci95"][0]),
                "CondLLHi": sci(w["cond_ll_ci95"][1]), "LooRMin": rv(w["loo_R_range"][0]),
                "LooRMax": rv(w["loo_R_range"][1])}
        gb = v42["groups"]["by_group"]
        ors = []
        for key, name in EG:
            t = gb[key].get("or_le500m_adj")
            out |= or_macros(name, t)
            if t:
                ors.append(t["odds_ratio"])
        out |= {"EGMin": orv(min(ors)), "EGMax": orv(max(ors)),
                "EGDikeCells": str(gb["dike_reclamation"]["n_cells_le500m"])}
        s = v42["settler_excluded"]
        out |= {"NSettlerAoi": str(s["n_aoi_dropped"]), "RSettler": rv(s["R"]),
                "RSettlerLoNinety": rv(s["ci90_R"][0]), "RSettlerHiNinety": rv(s["ci90_R"][1])}
    else:
        out["TIORPooledEmp"] = missing("v42")
    v47 = v47 or {}
    lf = (v47.get("lexfix") or {}).get("primary")
    if isinstance(lf, dict):
        out |= {"LexfixR": rv(lf["R"]), "LexfixRLoNinety": rv(lf["ci90_R"][0]),
                "LexfixRHiNinety": rv(lf["ci90_R"][1]), "LexfixVerdict": lf["decision"]}
        # B-4 の規則: 凍結したコードと直したコードの判定が全ての分割で同じなら、凍結したコードの実行を主とする
        pairs = [((v47.get("frozen") or {}).get(sp), (v47.get("lexfix") or {}).get(sp))
                 for sp in ("primary", "country", "temporal", "llmcut")]
        ok = [(a, b) for a, b in pairs if isinstance(a, dict) and isinstance(b, dict)]
        same = sum(a["decision"] == b["decision"] for a, b in ok)
        out["LexfixAgreement"] = (f"on all {_word(len(ok))} splits" if same == len(ok)
                                  else f"on {_word(same)} of the {_word(len(ok))} splits")
    else:
        out |= {k: missing("lexfix") for k in ("LexfixR", "LexfixRLoNinety", "LexfixRHiNinety", "LexfixVerdict",
                                               "LexfixAgreement")}
    fs = (v47.get("fold_seeds") or {}).get("fixed_code") or {}
    if fs.get("n_available", 0) >= 13 and "combined_ci90" in fs:
        out |= {"FSN": str(fs["n_available"]), "FSNWord": _word(fs["n_available"]), "FSRMin": rv(fs["R_min"]),
                "FSRMax": rv(fs["R_max"]), "FSCombLoNinety": rv(fs["combined_ci90"][0]),
                "FSCombHiNinety": rv(fs["combined_ci90"][1]), "FSNEquiv": str(fs["n_equivalent_null"]),
                "FSEquivPhrase": (f"all {fs['n_available']} assignments reached"
                                  if fs["n_equivalent_null"] == fs["n_available"]
                                  else f"{fs['n_equivalent_null']} of the {fs['n_available']} assignments reached")}
    else:
        out |= {k: missing(f"fold seeds ({fs.get('n_available', 0)} of 13)") for k in
                ("FSN", "FSNWord", "FSRMin", "FSRMax", "FSCombLoNinety", "FSCombHiNinety", "FSNEquiv", "FSEquivPhrase")}
    pm = v47.get("permutations")
    if isinstance(pm, dict):
        out |= {"PermN": str(pm["n_perms"]), "PermLo": rv(pm["R_perm_q025"]), "PermHi": rv(pm["R_perm_q975"]),
                "PermPct": f"{pm['percentile_of_real']:.0f}th", "PermRReal": rv(pm["R_real"]),
                "PermRMean": rv(pm["R_perm_mean"]),
                # 実の名前と置換の平均の差（R の尺度）と、置換が実以上だった割合（片側の順位）
                "PermDiff": rv(pm["R_real"] - pm["R_perm_mean"]), "PermShareGe": f"{pm['share_perm_ge_real']:.2f}"}
    else:
        out |= {k: missing("permutations") for k in ("PermN", "PermLo", "PermHi", "PermPct", "PermRReal",
                                                     "PermRMean", "PermDiff", "PermShareGe")}
    oc = v47.get("oracle_and_coefs")
    if isinstance(oc, dict):
        o = oc["oracle"]
        out |= {"OracleR": rv(o["R"]), "OracleRLo": rv(o["ci95_R"][0]), "OracleRHi": rv(o["ci95_R"][1]),
                "OracleRLoNinety": rv(o["ci90_R"][0]), "OracleRHiNinety": rv(o["ci90_R"][1]),
                "OracleVerdict": o["decision"]}
    else:
        out |= {k: missing("oracle") for k in ("OracleR", "OracleRLo", "OracleRHi", "OracleRLoNinety",
                                               "OracleRHiNinety", "OracleVerdict")}
    for key, name in (("notrace", "LVNoTrace"), ("nonwater", "LVNonWater")):
        lv = (v47.get("label_variants") or {}).get(key)
        if isinstance(lv, dict):
            out |= {f"{name}R": rv(lv["R"]), f"{name}Verdict": lv["decision"],
                    f"{name}RLoNinety": rv(lv["ci90_R"][0]), f"{name}RHiNinety": rv(lv["ci90_R"][1])}
        else:
            out |= {f"{name}{s}": missing(f"label variant {key}") for s in ("R", "Verdict", "RLoNinety",
                                                                            "RHiNinety")}
    return out


# %% B-5 の 10（v43・v44）
def _rate(r: dict | None) -> str:
    if not r or r.get("rate") is None:
        return "---"
    return f"{100 * r['rate']:.0f}\\,\\% ({100 * r['ci95'][0]:.0f}--{100 * r['ci95'][1]:.0f})"


def etym_macros(v43: dict | None, v44: dict | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if v43:
        out |= {"WDNItems": _n(v43["n_with_item"]), "WDNPOneThreeEight": _n(v43["n_with_p138"]),
                "WDNHydro": str(v43["n_hydro"]),
                "WDNHydroFR": str(v43["by_country"].get("FR", {}).get("n_hydro", 0)),
                "WDRecallLexPct": pct1(v43["recall_lexicon"]["rate"]),
                "WDRecallNreadPct": pct1(v43["recall_nread"]["rate"]),
                "WDRecallLex": _rate(v43["recall_lexicon"]), "WDRecallNread": _rate(v43["recall_nread"]),
                "WDPrecLex": _rate(v43["hydro_among_lexicon_wetside"]),
                "WDPrecNread": _rate(v43["hydro_among_nread_wetside"])}
    else:
        out["WDNItems"] = missing("v43")
    if v44:
        def w(blk: dict) -> str:
            return f"{100 * blk['est']:.0f}\\,\\% ({100 * blk['ci95'][0]:.0f}--{100 * blk['ci95'][1]:.0f})"

        sc, r1 = v44["status_counts"], v44.get("round1") or {}
        out |= {"WPNPopulation": _n(v44["n_population"]), "WPNStratumWet": _n(v44["stratum_sizes"]["wet"]),
                "WPNStratumNot": _n(v44["stratum_sizes"]["not"]), "WPNSampledNot": _n(v44["sample_sizes"]["not"]),
                "WPNSampled": _n(v44["n_sampled"]), "WPNWithSection": _n(v44["n_with_section"]),
                "WPNNoDeriv": _n(sc.get("no-derivation", 0)), "WPNDisputed": _n(sc.get("unclear", 0)),
                "WPRoundOneN": _n(r1["n_sampled"]) if r1 else missing("v44 round 1"),
                "WPRoundOneWithSection": _n(r1["n_with_section"]) if r1 else missing("v44 round 1"),
                "WPPrecNreadPct": pct1(v44["nread_wetside"]["precision"]["est"]),
                "WPPrecLex": w(v44["lex_wetside"]["precision"]), "WPRecallLex": w(v44["lex_wetside"]["recall"]),
                "WPPrecNread": w(v44["nread_wetside"]["precision"]),
                "WPRecallNread": w(v44["nread_wetside"]["recall"]),
                "WPPrecLexPct": pct1(v44["lex_wetside"]["precision"]["est"]),
                "WPPrecLexLo": pct1(v44["lex_wetside"]["precision"]["ci95"][0]),
                "WPPrecLexHi": pct1(v44["lex_wetside"]["precision"]["ci95"][1]),
                "WPRecallLexPct": pct1(v44["lex_wetside"]["recall"]["est"]),
                "WPRecallLexLo": pct1(v44["lex_wetside"]["recall"]["ci95"][0]),
                "WPRecallLexHi": pct1(v44["lex_wetside"]["recall"]["ci95"][1]),
                "WPRecallNreadPct": pct1(v44["nread_wetside"]["recall"]["est"]),
                "WPRateNotPct": pct1(v44["truth_rate_in_not"]["rate"]),
                "WPNCoded": _n(sc.get("coded", 0)),
                "EtymRetrievedUTC": v44.get("retrieved_utc") or missing("retrieval time")}
    else:
        out |= {k: missing("v44 (Wikipedia etymology)") for k in
                ("WPPrecLex", "WPRecallLex", "WPPrecNread", "WPRecallNread", "WPPrecLexPct", "WPPrecLexLo",
                 "WPPrecLexHi", "WPRecallLexPct", "WPRecallLexLo", "WPRecallLexHi", "WPRecallNreadPct",
                 "WPRateNotPct", "WPNCoded", "EtymRetrievedUTC", "WPNPopulation", "WPNStratumWet",
                 "WPNStratumNot", "WPNSampledNot", "WPNSampled", "WPNWithSection", "WPNNoDeriv", "WPNDisputed",
                 "WPRoundOneN", "WPRoundOneWithSection", "WPPrecNreadPct")}
    return out


# %% 図 4（本文）: 名前が加えるもの — 確証の比較・B-4・副の分割、事後の確認、fold の割り当てと置換の分布
# %% 図 1（本文）: 設計。旧版の図を、本文の幅で読める字の大きさと、改訂後の言い方（語のブロック、地名を含まない、
# T0 は標高と HAND が無い）で描き直す
def fig_design_nhess(n_blocks: str) -> None:
    from itertools import pairwise as _pw

    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig = plt.figure(figsize=(7.2, 3.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 48)
    ax.axis("off")
    blue, grey, light = "#1f6fb4", "0.4", "#e8f1f9"

    def box(x: float, y: float, w: float, h: float, text: str, *, fc: str = "white", ec: str = "0.3",
            ls: str = "-", size: float = 7.6) -> None:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.8", fc=fc, ec=ec, lw=0.8,
                                    ls=ls))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, linespacing=1.2)

    def arrow(x0: float, y0: float, x1: float, y1: float) -> None:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9, lw=0.9, color="0.3"))

    ax.text(0.5, 46, "a", fontsize=10, weight="bold")
    ax.text(3.5, 46, "Nested comparison within areas (primary, H2)", fontsize=8.4, weight="bold")
    y, h = 31, 10
    xs = [(1, 14, "Position\n(M0)"), (20, 23, "+ terrain, HAND,\nrivers, land cover,\nsettlement layout (M1)"),
          (49, 23, "+ regional composition\nof names: which names\noccur in the area (M1c)"),
          (78, 21.5, "+ the names of each\ncell's nearest\nsettlements (M3n)")]
    for x, w, txt in xs:
        box(x, y, w, h, txt, fc=light if x == 78 else "white", ec=blue if x == 78 else "0.3")
    for (x0, w0, _), (x1, _, _) in _pw(xs):
        arrow(x0 + w0 + 0.8, y + h / 2, x1 - 0.8, y + h / 2)
    ax.text(17.5, 28.5, "terrain increment $\\Delta_T$\n(the yardstick)", fontsize=7.4, ha="center", va="top",
            color=grey)
    ax.text(74.8, 28.5, "names increment $\\Delta_N$\n$R=\\Delta_N/\\Delta_T$", fontsize=7.4, ha="center", va="top",
            color=blue)
    box(78, 17, 21.5, 8, "placebo (M3np): the same\nnames reassigned among\nthe area's settlements", ec=grey, ls="--",
        size=7.0)
    ax.text(0.5, 14.2, "b", fontsize=10, weight="bold")
    ax.text(3.5, 14.2, "Terrain ladder (H3)", fontsize=8.4, weight="bold")
    rungs = ["T0\nno elev.\nor HAND", "T1\n1 km", "T2\n250 m", "T3\n30 m\nstats", "T4 = M1\n+ HAND"]
    for k, r in enumerate(rungs):
        box(1 + k * 6.9, 1.5, 6.0, 9, r, size=6.6)
    ax.text(36.5, 14.2, "c", fontsize=10, weight="bold")
    ax.text(39.5, 14.2, "Three readings of names", fontsize=8.4, weight="bold")
    reads = [((f"lexicon (N-lex):\nwater, wetland, riverine,\ncrossing, lowland,\nhighland elements;\n{n_blocks} "
               "language blocks"), 36.8),
             ("character model\n(N-learn):\nlearned on training\nareas only", 53.6),
             ("language model\n(N-read): name\nelements; tested on\nfloods after its\ntraining cut-off", 70.4)]
    for txt, x in reads:
        box(x, 1.5, 14.8, 10.5, txt, size=6.4)
    ax.text(87.5, 14.2, "d", fontsize=10, weight="bold")
    ax.text(90.5, 14.2, "Held out", fontsize=8.4, weight="bold")
    box(87.5, 1.5, 12, 10, "places (primary)\ncountries\nlater floods\nfloods after the\nLM cut-off", size=6.6)
    fig.savefig(HERE / "fig_design.pdf")
    plt.close(fig)


XMAX_FIG4 = 0.08


def _rows90(primary: dict, extra: dict, v36: dict | None) -> list[tuple[str, float | None, list | None, list | None, str]]:
    """図 4 の行: (表示名, R̂, 90% 区間, 95% 区間, 様式)。区間は固定マージン（チャネルと H1 は v36）。"""
    ch = primary.get("channels", {})

    def blk(b: dict | None, rk: str = "R", c95: str = "ci95_R", c90: str = "ci90_R"):
        if not b or b.get(rk) is None:
            return None
        return b[rk], b.get(c90), b.get(c95)

    def fm(comp: str, b: dict | None, rk: str, c95: str):
        f = g.fixed_margin(v36, "primary", comp)
        return (f["R"], f["ci90_R"], f["ci95_R"]) if f else blk(b, rk, c95, "ci90_R")

    spec = [("Names within areas (primary, H2)", blk(primary.get("primary")), "main"),
            ("   Lexicon only (N-lex)", fm("channel_lexicon", ch.get("M3nx"), "R", "ci95_R"), "main"),
            ("   Character model only (N-learn)", fm("channel_learned", ch.get("M3nl"), "R", "ci95_R"), "main"),
            ("   Placebo: names reassigned", blk(primary.get("placebo"), "R_placebo", "ci95_R_placebo",
                                                 "ci90_R_placebo"), "main"),
            ("Position only, within areas (H1)", fm("H1", primary.get("H1"), "ratio_to_terrain", "ci95_ratio"),
             "main")]
    if (met := extra.get("country")) and met.get("primary"):
        spec.append(("Countries held out", blk(met["primary"]), "main"))
    if (met := extra.get("temporal")) and met.get("primary"):
        und = met["primary"]["decision"] == "undefined-denominator"
        spec.append(("Time held out (floods from 2024)" + (" *" if und else ""), blk(met["primary"]),
                     "undefined" if und else "main"))
    if (met := extra.get("llmcut")) and met.get("H4"):
        spec.append(("Language model, floods after its cut-off (H4)", blk(met["H4"]), "main"))
        spec.append(("   Lexicon and character model, same floods", blk(met.get("primary")), "main"))
        if met.get("H4_total"):
            spec.append(("   Language model incl. its composition (secondary)", blk(met["H4_total"]), "secondary"))
    spec += [("Total incl. regional composition (secondary)", blk(primary.get("total")), "secondary"),
             ("Flexible model refitted with names (secondary)", blk(primary.get("flexible")), "secondary")]
    return [(n, v[0], v[1], v[2], st) for n, v, st in spec if v is not None]


def fig_main(primary: dict, extra: dict, v36: dict | None, v41: dict | None, v42: dict | None,
             v47: dict | None, v48: dict | None = None, v52: dict | None = None) -> None:
    rows = _rows90(primary, extra, v36)
    v47 = v47 or {}
    lf = (v47.get("lexfix") or {}).get("primary")
    if isinstance(lf, dict):
        rows.insert(1, ("   Same, Norwegian names read correctly", lf["R"], lf["ci90_R"], lf["ci95_R"], "main"))
    post: list[tuple[str, float | None, list | None, list | None, str]] = []
    if v41:
        c = v41["C1_gfm"]
        post.append(("Independent labels (GFM): rescored", c["R"], c["ci90_R"], c["ci95_R"], "post"))
    c2 = v47.get("gfm_c2")
    if isinstance(c2, dict):
        post.append(("Independent labels (GFM): refitted", None, None, None, "post"))
    for key, lab in (("notrace", "Flood traces counted as not flooded"),
                     ("nonwater", "Flooded share over land only")):
        lv = (v47.get("label_variants") or {}).get(key)
        if isinstance(lv, dict):
            post.append((lab, lv["R"], lv["ci90_R"], lv["ci95_R"], "post"))
    met = g._json(DATA / "runs" / "world_v2_actgrp" / "metrics_primary.json")
    if met:
        q = met["primary"]
        post.append(("Flood events held out whole", q["R"], q["ci90_R"], q["ci95_R"], "post"))
    if v42:
        st = v42["settler_excluded"]
        post.append(("Without the Americas and Oceania", st["R"], st["ci90_R"],
                     [v / st["yardstick_M0_M1"] for v in st["ci95_delta"]], "post"))
    if v48:
        a = v48["aoi_weighted"]["frozen"]
        post.append(("Every AOI weighted equally", a["R"], a["ci90_R"], a["ci95_R"], "post"))
    v52 = v52 or {}
    for key, cmp_, lab in (("strict", "strict", "Strict reading (water bodies, wetland, riverine)"),
                           ("nread", "H4", "Language-model reading, primary split"),
                           ("gn", "primary", "Gazetteer cleaned of city sections and duplicates"),
                           ("nohydro", "primary", "Without AOIs lacking a hydrography layer"),
                           ("typed", "primary", "Without AOIs lacking a flood type")):
        b = v52.get(key)
        c = b.get(cmp_) if isinstance(b, dict) else None
        if isinstance(c, dict):
            post.append((lab, c["R"], c["ci90_R"], c["ci95_R"], "post"))
    rows_all = rows + post
    fs = (v47.get("fold_seeds") or {}).get("fixed_code") or {}
    pm = v47.get("permutations")
    g.SOURCE["Fig4"] = [{"comparison": n.strip(), "R": r, "ci90_low": c9[0] if c9 else None,
                         "ci90_high": c9[1] if c9 else None, "ci95_low": c5[0] if c5 else None,
                         "ci95_high": c5[1] if c5 else None, "group": st} for n, r, c9, c5, st in rows_all]
    n = len(rows_all) + 2
    fig, ax = plt.subplots(figsize=(6.6, 0.155 * n + 0.85))  # 図と説明が原稿の 1 ページに収まる高さ
    ax.axvspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axvline(0, color="0.4", lw=0.8)
    ys = np.arange(n)[::-1]
    style = {"main": ("#1f6fb4", "-", True), "secondary": ("0.45", "-", False),
             "undefined": ("0.45", "--", False), "post": ("#8c5a2b", "-", False)}
    rows_y = ys[:-2]
    for y, (_, r, c9, c5, st) in zip(rows_y, rows_all, strict=True):
        col, ls, filled = style[st]
        if r is None:
            ax.text(0.0, y, "R undefined: terrain did not improve on position", ha="center", va="center",
                    fontsize=6, color=col, style="italic")
            continue
        if c5 is not None:  # 95% は細く、90%（判定に使う）は太く
            ax.plot(c5, [y, y], color=col, lw=0.7, ls=ls)
        if c9 is not None:
            ax.plot(c9, [y, y], color=col, lw=2.2, ls=ls, solid_capstyle="butt")
        ax.plot([r], [y], "o", ms=4.0, color=col, mfc=col if filled else "white", zorder=3)
        hi = max(r, *(c5[1:] if c5 else []), *(c9[1:] if c9 else []))
        if hi > XMAX_FIG4:  # 軸の外に出る区間は端に矢印と上限の値を
            ax.annotate(f"{hi:.2f}", xy=(XMAX_FIG4, y), xytext=(XMAX_FIG4 - 0.014, y + 0.42), fontsize=5.5,
                        color=col, ha="center", arrowprops={"arrowstyle": "->", "color": col, "lw": 0.6})
    y_fs, y_pm = 1, 0
    if fs.get("R"):
        rs = np.array(list(fs["R"].values()), dtype=float)
        ax.plot(rs, np.full(len(rs), y_fs), "o", ms=3.2, color="#8c5a2b", mfc="white", alpha=0.9)
    if isinstance(pm, dict):
        ax.plot([pm["R_perm_q025"], pm["R_perm_q975"]], [y_pm, y_pm], color="0.7", lw=5, solid_capstyle="butt")
        ax.plot([pm["R_real"]], [y_pm], "o", ms=4.2, color="#8c5a2b", mfc="#8c5a2b")
    n_perm = pm["n_perms"] if isinstance(pm, dict) else 0
    labels = [r[0] for r in rows_all] + [f"{fs.get('n_available', 0)} assignments of places to folds (range)",
                                         f"Names reassigned {n_perm} times, 95% range (grey); real names (dot)"]
    ax.set_yticks(list(rows_y) + [y_fs, y_pm])
    ax.set_yticklabels(labels, fontsize=6.4)
    for y, (nm, *_) in zip(rows_y, rows_all, strict=True):
        if nm.startswith(("Countries held out", "Independent labels (GFM): rescored")):
            ax.axhline(y + 0.5, color="0.8", lw=0.6)
    ax.axhline(y_fs + 0.5, color="0.8", lw=0.6)
    ax.set_xlabel("R (improvement from names / improvement from terrain)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlim(-XMAX_FIG4, XMAX_FIG4)
    ax.set_ylim(-0.8, ys[0] + 0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(HERE / "fig_main.pdf")
    plt.close(fig)


# %% 図 6（本文）: 地形情報の梯子。B-6 の 6(b) の T0−（川の距離と CCI の水域・湿地なし、事後）を左端に
def fig_ladder_b6(ladder: dict, terrain: float, v52: dict | None) -> None:
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    d = list(np.array(ladder["delta_names"]) / terrain)
    c = [list(np.array(x) / terrain) for x in ladder["ci95_delta_names"]]
    labels = ["T0\nno elevation", "T1\n1 km", "T2\n250 m", "T3\n30 m stats", "T4\n+HAND"]
    t0x = ((v52 or {}).get("t0x") or {}).get("T0x") if isinstance((v52 or {}).get("t0x"), dict) else None
    x = np.arange(len(d))
    if t0x:  # T0− は全体の物差しの単位（同じ実行の M0 → M1）で、確証と同じ固定マージン
        x = x + 1
        ax.errorbar([0], [t0x["R"]], yerr=[[t0x["R"] - t0x["ci95_R"][0]], [t0x["ci95_R"][1] - t0x["R"]]],
                    fmt="s", color="#8c5a2b", ms=4, lw=1.2, capsize=2, mfc="white")
        labels = ["T0$^-$\nno river,\nno water LC"] + labels
    g.SOURCE["Fig6"] = ([{"rung": "T0-", "R": t0x["R"], "ci95_low": t0x["ci95_R"][0], "ci95_high": t0x["ci95_R"][1]}]
                        if t0x else []) + [{"rung": rg, "R": float(d[k]), "ci95_low": float(c[k][0]),
                                            "ci95_high": float(c[k][1])} for k, rg in enumerate(ladder["rungs"])]
    ax.axhspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axhline(0, color="0.4", lw=0.8)
    dd, cc = np.array(d), np.array(c)
    ax.errorbar(x, dd, yerr=[dd - cc[:, 0], cc[:, 1] - dd], fmt="o-", color="#1f6fb4", ms=4, lw=1.2, capsize=2)
    ax.set_xticks(list([0] if t0x else []) + list(x))
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("names within areas\n(fraction of terrain increment)", fontsize=8)
    ax.set_ylim(-0.08, 0.08)
    ax.tick_params(axis="y", labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(HERE / "fig_ladder.pdf")
    plt.close(fig)


# %% 図 5（本文）: 名前が記述するものと、集落の近くの浸水
def fig_names(v37: dict, v39: dict | None, v41: dict | None, v42: dict | None,
              v50: dict | None = None, v52: dict | None = None) -> None:
    fig = plt.figure(figsize=(7.2, 5.6), layout="constrained")
    axs = fig.subplots(2, 2).ravel()
    classes = list(g.CLASS_LABELS)
    ys = np.arange(len(classes))[::-1]
    blue, brown, grey, green = "#1f6fb4", "#8c5a2b", "0.45", "#2b8c4a"
    pc = v37["positive_control"]
    # a: 集落の点での HAND ≤ 2 m の割合の差（塗り）と全セル（白抜き）
    ax = axs[0]
    for y, c in zip(ys, classes, strict=True):
        col = brown if c == "high" else blue
        for scale, dy, filled in (("settlements", 0.15, True), ("cells", -0.15, False)):
            r = pc[scale][c]["low_hand"]
            if not r:
                continue
            ax.plot(100 * np.array(r["ci95"]), [y + dy, y + dy], color=col if filled else grey, lw=1.4)
            ax.plot([100 * r["estimate"]], [y + dy], "o", ms=4, color=col if filled else grey,
                    mfc=(col if filled else "white"))
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([g.CLASS_LABELS[c] for c in classes], fontsize=7)
    ax.set_xlabel("difference in share of cells with ground ≤ 2 m above\ndrainage (pct. points; filled: at the "
                  "settlement, open: all cells)", fontsize=7)
    # b: 距離帯ごとの差（HAND・浸水）
    ax = axs[1]
    bands = list(g.BAND_LABELS)
    x = np.arange(len(bands))
    series = [("HAND ≤ 2 m (no flood labels)", [pc["bands"][b]["low_hand"] for b in bands], grey, "s", True),
              ("flooded (EMS)", [v37["flood_contrast"]["by_band"]["capped"][b]["y"] for b in bands], blue, "o",
               True),
              ("flooded, net of terrain model",
               [v37["flood_contrast"]["by_band"]["capped"][b]["resid_M1c"] for b in bands], blue, "o", False)]
    for k, (lab, rs, col, mk, filled) in enumerate(series):
        dx = (k - 1) * 0.18
        est = np.array([100 * r["estimate"] for r in rs])
        lo = np.array([100 * r["ci95"][0] for r in rs])
        hi = np.array([100 * r["ci95"][1] for r in rs])
        ax.errorbar(x + dx, est, yerr=[est - lo, hi - est], fmt=mk, ms=4, color=col, lw=1.2,
                    mfc=col if filled else "white", capsize=0, label=lab)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(["≤250 m", "250–500 m", "0.5–1 km", ">1 km"], fontsize=7)
    ax.set_xlabel("distance to the nearest settlement", fontsize=7)
    ax.set_ylabel("water-side − other names\n(pct. points)", fontsize=7)
    ax.legend(fontsize=6, frameon=False, loc="upper right")
    # c: 水系側の名の 500 m 以内の浸水のオッズ比（地形で調整）: セルの選び方を順に狭める（B-6 の 4）
    ax = axs[2]
    chain = []
    if v39:
        chain.append(("EMS, all analysed AOIs", (v39.get("wetside_adjusted_M1") or {}).get("terms", {})
                      .get("wetside:le500m"), blue, True))
    if v50:
        sel = v50["aoi_vs_cell_selection"]
        chain += [("EMS, all cells of the GFM AOIs", sel["ems_192aoi_all_cells"], blue, False),
                  ("EMS, GFM-assessable cells", sel["ems_gfm_cells"], grey, False),
                  ("GFM, same cells", sel["gfm_gfm_cells"], green, True)]
        v1 = v50.get("v1")
        if v1:
            chain += [("EMS, cells assessable in V1", v1["odds_ems"]["wetside_adjusted_M1"]["terms"]["wetside:le500m"],
                       grey, False),
                      ("GFM V1 (seasonal water assessable)",
                       v1["odds_v1"]["wetside_adjusted_M1"]["terms"]["wetside:le500m"], green, False)]
    yc = np.arange(len(chain))[::-1]
    for y, (lab, t, col, filled) in zip(yc, chain, strict=True):
        if not t:
            continue
        lo, hi = (max(v, 0.05) for v in t["ci95_odds_ratio"])
        ax.plot([lo, min(hi, 20)], [y, y], color=col, lw=1.4)
        ax.plot([t["odds_ratio"]], [y], "o", ms=4, color=col, mfc=col if filled else "white")
    g.SOURCE["Fig5c"] = [{"labels_and_cells": lab, "odds_ratio": t["odds_ratio"] if t else None,
                          "ci95": t["ci95_odds_ratio"] if t else None} for lab, t, *_ in chain]
    ax.set_xscale("log")
    ax.axvline(1, color="0.4", lw=0.8)
    ax.set_xlim(0.25, 2.5)
    ax.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.25", "0.5", "1", "2"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(yc)
    ax.set_yticklabels([c[0] for c in chain], fontsize=6.5)
    ax.set_xlabel("odds ratio of flooding within 500 m of settlements\nwith water-side names, given the terrain model",
                  fontsize=7)
    # d: plasmode（B-6 の 2）: 注入した効果のオッズ比ごとの増分（置換の平均 + 利得）と、観測の 90% 上限・余白
    ax = axs[3]
    pl = (v52 or {}).get("plasmode")
    if isinstance(pl, dict):
        ors = np.exp(np.concatenate([[0.0], np.asarray(pl["deltas"], dtype=float)]))
        for f, lab, col, ls in (("A", "all ground nearest the name", blue, "-"),
                                ("B", "within 500 m only", brown, "-"),
                                ("CA", "all ground, truly water-related", blue, ":"),
                                ("CB", "500 m, truly water-related", brown, ":")):
            r = pl["forms"][f]
            gain = np.concatenate([[0.0], np.asarray(r["gain"], dtype=float)])
            se = np.concatenate([[0.0], np.asarray(r["gain_se"], dtype=float)])
            yv = 1e4 * gain  # 利得そのもの（余白の線との交点 = δ_M、観測の上限 − 置換の平均の線との交点 = δ_A）
            ax.plot(ors, yv, ls, color=col, lw=1.2, marker="o" if ls == "-" else None, ms=3, label=lab)
            ax.fill_between(ors, yv - 1e4 * se, yv + 1e4 * se, color=col, alpha=0.15, lw=0)
        room = pl["U90"] - pl["delta_bs_null_perm_mean"]
        ax.axhline(1e4 * pl["margin"], color="0.2", lw=0.9, ls="--", label="margin (R = 0.05)")
        ax.axhline(1e4 * room, color="0.2", lw=0.9, label="observed limit above reassigned")
        # 交点（余白・上限のオッズ比）が見える範囲に絞る（大きなオッズ比の利得は上に出る）
        ax.set_ylim(-0.3, 2.0 * 1e4 * pl["margin"])
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20]))
        ax.xaxis.set_major_formatter(FixedFormatter(["1", "2", "5", "10", "20"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlabel("injected odds ratio of flooding", fontsize=7)
        ax.set_ylabel("gain in the names increment, ΔBS (×10$^{-4}$)", fontsize=7)
        ax.legend(fontsize=5.2, frameon=False, loc="center right", bbox_to_anchor=(1.04, 0.32))
        g.SOURCE["Fig5d"] = {"odds_ratio": ors.tolist(), "null": pl["delta_bs_null_perm_mean"], "U90": pl["U90"],
                             "U90_minus_null": room, "margin": pl["margin"],
                             "gain": {f: pl["forms"][f]["gain"] for f in pl["forms"]}}
    else:
        ax.text(0.5, 0.5, "plasmode results pending", ha="center", va="center", fontsize=7, transform=ax.transAxes)
    for a, lab in zip(axs, "abcd", strict=True):
        a.tick_params(labelsize=6.5)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.set_title(f"({lab})", fontsize=9, fontweight="bold", loc="left")
    fig.savefig(HERE / "fig_names.pdf")
    plt.close(fig)


# %% 図 3（本文）: 例の AOI（選び方は結果に依らない: GFM で評価でき、場所の上限で切られていない AOI のうち、
# 地図が読めるよう集落が 30 以下のもののうち、辞書が水系側と読む集落の最も多いもの。同数なら解析したセルの
# 多いもの。浸水ラベル・GFM で評価できる割合・名前の効果は選び方に使わない）
MAX_EXAMPLE_SETTLEMENTS = 30
def pick_example(t: pd.DataFrame, gfm_summary: dict) -> str | None:
    inc = {a for a, s in gfm_summary["aois"].items() if s.get("include")}
    n_tab = t.groupby("aoi").size()
    best, best_key = None, (-1, -1)
    gn = pd.read_parquet(DATA / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "asciiname", "cc", "lat", "lon", "fclass"]).drop_duplicates("gid")
    gn = gn[gn["fclass"] == "P"]  # 解析と同じ集落（P 類）だけ。以前は全ての地物を数えていた（2026-10-27 に修正）
    lex = g_lexicon()
    for a in sorted(inc):
        d = DATA / "global" / "labels" / a
        if not (d / "cells.parquet").is_file():
            continue
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        evaluable = int(n_tab.get(a, 0))
        if evaluable == 0 or evaluable < 0.99 * gfm_summary["aois"][a].get("n_eval", 0):
            continue  # 場所の上限で間引かれた AOI は予測の図が欠けるので除く
        sw = settlements_in(meta, gn)
        if len(sw) > MAX_EXAMPLE_SETTLEMENTS:
            continue
        n_wet = sum(any(lex.match(r.asciiname, r.cc)[c] for c in ("water", "wet", "riverine", "crossing", "low"))
                    for r in sw.itertuples())
        if (n_wet, evaluable) > best_key:
            best, best_key = a, (n_wet, evaluable)
    return best


def g_lexicon():
    sys.path.insert(0, str(REPO / "backend" / "packages" / "topsem"))
    from topsem.world.lexicon import Lexicon
    return Lexicon(REPO / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")


def settlements_in(meta: dict, gn: pd.DataFrame) -> pd.DataFrame:
    from pyproj import Transformer
    tr = Transformer.from_proj("+proj=longlat +datum=WGS84", meta["crs_proj4"], always_xy=True)
    gr = meta["grid"]
    x, y = tr.transform(gn["lon"].to_numpy(), gn["lat"].to_numpy())
    col = (x - gr["x0"]) / gr["cell_m"]
    row = (gr["y1"] - y) / gr["cell_m"]
    ok = (col >= 0) & (col < gr["nx"]) & (row >= 0) & (row < gr["ny"])
    return gn.loc[ok].assign(col=col[ok], row=row[ok])


def fig_example(aoi: str, t: pd.DataFrame) -> dict[str, str]:
    d = DATA / "global" / "labels" / aoi
    meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
    cells = pd.read_parquet(d / "cells.parquet", columns=["ix", "iy"])
    feat = pd.read_parquet(d / "features.parquet", columns=["hand_mean"])
    if len(feat) != len(cells):
        msg = f"{aoi}: features and cells differ in length"
        raise SystemExit(msg)
    gfm = pd.read_parquet(d / "gfm.parquet").merge(t.loc[t["aoi"] == aoi, ["ix", "iy"]], on=["ix", "iy"])
    ny, nx = meta["grid"]["ny"], meta["grid"]["nx"]

    def grid(vals: np.ndarray, ix: np.ndarray, iy: np.ndarray) -> np.ndarray:
        a = np.full((ny, nx), np.nan)
        a[iy.astype(int), ix.astype(int)] = vals
        return a

    s = t[t["aoi"] == aoi]
    ems = grid(s["y"].to_numpy(dtype=float), s["ix"].to_numpy(), s["iy"].to_numpy())  # 解析したセルの EMS
    gf = grid(np.where(gfm["gfm_eval"], gfm["gfm_y"], np.nan), gfm["ix"].to_numpy(), gfm["iy"].to_numpy())
    hand = grid(feat["hand_mean"].to_numpy(), cells["ix"].to_numpy(), cells["iy"].to_numpy())
    hand[np.isnan(ems)] = np.nan
    diff = grid((s["p_M3n"] - s["p_M1c"]).to_numpy(), s["ix"].to_numpy(), s["iy"].to_numpy())
    gn = pd.read_parquet(DATA / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "asciiname", "cc", "lat", "lon", "fclass"]).drop_duplicates("gid")
    st = settlements_in(meta, gn[gn["fclass"] == "P"])  # 解析と同じ集落（P 類）だけ
    lex = g_lexicon()
    st = st.assign(wet=[any(lex.match(r.asciiname, r.cc)[c] for c in ("water", "wet", "riverine", "crossing", "low"))
                        for r in st.itertuples()])
    aspect = ny / nx
    # 4 枚を横に並べる（2×2 だと縦長の AOI で図が原稿の 1 ページを超える）
    fig, axs = plt.subplots(1, 4, figsize=(7.0, min(4.2, 1.62 * aspect + 1.0)), layout="constrained")
    from matplotlib.colors import ListedColormap
    lab_cmap = ListedColormap(["0.82", "#1f6fb4"])
    titles = ["(a) EMS labels", "(b) GFM labels", "(c) HAND", "(d) change from names"]
    for ax, arr in ((axs[0], ems), (axs[1], gf)):
        ax.imshow(arr, cmap=lab_cmap, vmin=0, vmax=1, interpolation="nearest")
    hmax = float(np.nanquantile(hand, 0.98)) or 1.0
    im = axs[2].imshow(np.clip(hand, 0, hmax), cmap="YlOrBr_r", interpolation="nearest")
    cb = fig.colorbar(im, ax=axs[2], orientation="horizontal", shrink=0.9, pad=0.02)
    cb.set_label(f"HAND (m; clipped at {hmax:.0f} m)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    lim = float(np.nanquantile(np.abs(diff), 0.99)) or 1e-3
    im = axs[3].imshow(diff, cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    cb = fig.colorbar(im, ax=axs[3], orientation="horizontal", shrink=0.9, pad=0.02)
    cb.set_label("M3n − M1c (probability)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    cb.ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    for ax in axs.ravel():
        ax.plot(st.loc[~st["wet"], "col"], st.loc[~st["wet"], "row"], "k^", ms=3, mfc="white", mew=0.6)
        ax.plot(st.loc[st["wet"], "col"], st.loc[st["wet"], "row"], "^", ms=3.6, color="#e66101", mew=0.6)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax, title in zip(axs.ravel(), titles, strict=True):
        ax.set_title(title, fontsize=7.5, loc="left")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(color="#1f6fb4", label="flooded"), Patch(color="0.82", label="not flooded"),
               Line2D([], [], marker="^", ls="", color="#e66101", label="settlement, water-side name"),
               Line2D([], [], marker="^", ls="", color="k", mfc="white", label="other settlement")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    fig.savefig(HERE / "fig_example.pdf", bbox_inches="tight", pad_inches=0.04)  # 枠の外の集落の印で題が欠けないように
    plt.close(fig)
    act = meta["activation"]
    country = ", ".join(act.get("countries") or [])
    tag = aoi.replace("_", "\\_")
    name = act.get("name", "")
    name = name.replace("- ", ", ")
    where = name if country and country in name else f"{name}, {country}"
    # 補足の表: 例の AOI の集落名と辞書の読み（B-6 の 5）
    from topsem.world.lexicon import _SPLIT, _hit
    rows_ex = []
    for r in st.sort_values("asciiname").itertuples():
        toks = [t_ for t_ in _SPLIT.split(str(r.asciiname).lower()) if t_]
        ents = [e for e in lex.entries_for(r.cc) if any(_hit(t_, e) for t_ in toks)]
        forms = ", ".join(sorted({f"\\emph{{{e.form}}} ({e.pos}, {g.CLASS_LABELS.get(e.cat, e.cat)})" for e in ents}))
        rows_ex.append(f"{r.asciiname} & {r.cc} & {'yes' if r.wet else 'no'} & {forms or '---'} \\\\")
    (HERE / "app_table_b6_example.tex").write_text(_tbl(
        "tab:b6-example", f"Settlements of the example AOI of Fig.~3 ({tag}) and their readings by the lexicon "
        "(populated places of the gazetteer within the AOI's grid; water-side: water, wetland, riverine, crossing or "
        "lowland element).", "llcp{9cm}", "Name & Country & Water-side & Matched forms (position, class)",
        rows_ex), encoding="utf-8")
    out = {"ExampleAoiName": f"{tag}: {where}, {act.get('activation_time', '')[:4]}",
           "ExampleNSettlements": str(len(st)), "ExampleNWetSettlements": str(int(st['wet'].sum())),
           "ExampleHandClip": f"{hmax:.0f}"}
    # 図 3 の注: EMS で浸水した解析セルのうち GFM で評価できない割合と、その主な成分（B-6 の 4 の gfm_v1）
    ga = gfm.set_index(["ix", "iy"])
    fl = s.loc[s["y"] == 1, ["ix", "iy"]]
    ev_ = ga.reindex(pd.MultiIndex.from_frame(fl))["gfm_eval"].fillna(value=False).astype(bool).to_numpy()
    out["ExampleFloodedNotGfmPct"] = pct(1 - ev_.mean()) if len(ev_) else missing("example flooded cells")
    v1p = d / "gfm_v1.parquet"
    if v1p.is_file():
        v1 = pd.read_parquet(v1p).merge(fl, on=["ix", "iy"])
        v1 = v1.loc[~v1["gfm_eval"].astype(bool)]
        comp = {k: float(v1[f"c_{k}"].mean()) for k in ("excl", "perm", "seas", "nodata")} if len(v1) else {}
        names = {"excl": "the exclusion mask", "perm": "permanent reference water",
                 "seas": "seasonal reference water", "nodata": "no GFM image coverage"}
        if comp:
            k = max(comp, key=comp.get)
            out["ExampleMaskMain"] = names[k]
            out["ExampleMaskMainPct"] = pct(comp[k])
    else:
        out["ExampleMaskMain"] = missing("example GFM components")
        out["ExampleMaskMainPct"] = missing("example GFM components")
    return out


# %% 付録の表: 旧版の補足表を Copernicus の書式に（\tophline 等、Supplementary の参照を付録に）
def copernicus_tables() -> list[str]:
    names = []
    for f in sorted(HERE.glob("si_table_*.tex")):
        s = f.read_text(encoding="utf-8")
        wide = "\\begin{sidewaystable}" in s  # 旧版で横置きの表は横置きのまま（幅の上限は \\textheight）
        lim = "\\textheight" if wide else "\\textwidth"
        s = (s.replace("\\toprule", "\\tophline").replace("\\midrule", "\\middlehline")
             .replace("\\bottomrule", "\\bottomhline").replace("\\begin{table}[h!]", "\\begin{table*}[t]")
             .replace("\\end{table}", "\\end{table*}").replace("\\begin{sidewaystable}", "\\begin{sidewaystable*}")
             .replace("\\end{sidewaystable}", "\\end{sidewaystable*}")
             # はみ出す表だけを縮める（自然な幅が上限以下ならそのまま）
             .replace("\\begin{adjustbox}{max width=\\linewidth}",
                      f"\\resizebox{{\\ifdim\\width>{lim}{lim}\\else\\width\\fi}}{{!}}{{%")
             .replace("\\end{adjustbox}", "}").replace("Supplementary Tables~", "Tables~")
             .replace("Supplementary Table~", "Table~")
             .replace("Supplementary Note~10", "Appendix~\\ref{app:prereg}")
             .replace("(Methods)", "(Sect.~\\ref{sec:methods})").replace("\\centering\\footnotesize", "\\footnotesize"))
        s = re.sub(r"Supplementary Note~\d+", "Appendix~\\\\ref{app:results}", s)
        out = HERE / f.name.replace("si_table_", "app_table_")
        out.write_text(s, encoding="utf-8")
        names.append(out.name)
    return names


# %% 付録 A: 時系列（git のコミット時刻・ラン記録の作成時刻・イメージのタグ = ビルド時刻 ns）
def _git_time(*args: str) -> str | None:
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%cI", *args], cwd=REPO,
                             capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    if not out:
        return None
    import datetime as dt
    return dt.datetime.fromisoformat(out).astimezone(dt.UTC).strftime("%Y-%m-%dT%H:%MZ")


def _first_commit_time(path: str) -> str | None:
    try:
        out = subprocess.run(["git", "log", "--diff-filter=A", "--format=%cI", "--", path], cwd=REPO,
                             capture_output=True, text=True, check=True).stdout.strip().splitlines()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    if not out:
        return None
    import datetime as dt
    return dt.datetime.fromisoformat(out[-1]).astimezone(dt.UTC).strftime("%Y-%m-%dT%H:%MZ")


AMEND = {"BOne": "B-1_s6_group_buffer_fix.md", "BTwo": "B-2_posthoc_review_analyses.md",
         "BThree": "B-3_gfm_independent_labels.md", "BFour": "B-4_lexicon_norway_yaml_fix.md",
         "BFive": "B-5_blind_review_analyses.md", "BSix": "B-6_blind_review2_analyses.md"}
CONFIRM_IMAGE_NS = 1792961983094088000  # 確証解析のイメージ（amendment B-1 に記録したタグ）


def _sim_last_utc() -> str:
    """検出力の模擬の最後の反復の記録時刻（凍結の前。要約は凍結の後に同じ反復から集計し直した）。"""
    ts = []
    for tag in ("power_v3_null", "power_v3_sesoi"):
        for f in sorted((DATA / "runs" / tag).glob("sim_meta_*.json")):
            ts.append((g._json(f) or {}).get("created_utc") or "")
    ts = [t for t in ts if t]
    return max(ts)[:16] + "Z" if ts else missing("power simulation records")


def timeline(m: dict[str, str]) -> dict[str, str]:
    """付録 A の時系列（UTC、時刻の順）。改訂 B-k は本文では Amendment k と呼ぶ。"""
    import datetime as dt
    out: dict[str, str] = {}
    for k, f in AMEND.items():
        t = _first_commit_time(f"docs/protocol/amendments/{f}")
        out[f"AmendTime{k}"] = t or missing(f"amendment {f}")
    content = _git_time(m.get("FreezeContentCommit", "HEAD"))
    mark = _git_time(m.get("FreezeMarkCommit", "HEAD"))
    img = dt.datetime.fromtimestamp(CONFIRM_IMAGE_NS / 1e9, tz=dt.UTC).strftime("%Y-%m-%dT%H:%MZ")
    pw = g._json(DATA / "runs" / "power_v3_null" / "power_summary.json") or {}
    conf = g._json(DATA / "runs" / "world_v2" / "summary_primary.json") or {}
    gs = g._json(DATA / "global" / "gfm_summary.json") or {}
    blind1 = _first_commit_time("docs/research/reviews/2026-10-26_paper2_blind1_editor.md")
    blind2 = _first_commit_time("docs/research/reviews/2026-10-27_paper2_blind2_editor.md")
    addendum = _git_time("--grep=B-6 追記 1", "--", "docs/protocol/amendments/B-6_blind_review2_analyses.md")
    addenda23 = [_git_time(f"--grep=B-6 追記 {k}", "--", "docs/protocol/amendments/B-6_blind_review2_analyses.md")
                 for k in (2, 3)]
    out["PowerSimLastUTC"] = _sim_last_utc()
    out["PowerReaggUTC"] = (pw.get("created_utc") or "")[:16] + "Z"
    rows = [("EMS product inventory retrieved", m.get("InvRetrievedUTC", "")),
            ("GeoNames gazetteer retrieved", m.get("GeoNamesUTC", "")),
            ("Power simulations on synthetic labels completed (last replicate)", out["PowerSimLastUTC"]),
            ("Protocol content fixed, including the power results (commit " + m.get("FreezeContentCommit", "")
             + ")", content or ""),
            ("Protocol marked frozen (commit " + m.get("FreezeMarkCommit", "") + ", tag prereg-v2-frozen)",
             mark or ""),
            ("Container image of the confirmatory analysis built", img),
            ("Confirmatory analysis aggregated (primary split)", (conf.get("created_utc") or "")[:16] + "Z"),
            ("Amendment 1 (sensitivity analysis 6 corrected)", out["AmendTimeBOne"]),
            ("Amendment 2 (post hoc analyses after LLM-agent critiques with access to the records)",
             out["AmendTimeBTwo"]),
            ("Power summary re-aggregated from the same replicates (code update)", out["PowerReaggUTC"]),
            ("Amendment 3 (independent GFM labels)", out["AmendTimeBThree"]),
            ("First LLM-agent critiques of the manuscript recorded", blind1 or ""),
            ("Amendment 4 (Norwegian names; rerun)", out["AmendTimeBFour"]),
            ("Amendment 5 (post hoc analyses after those critiques)", out["AmendTimeBFive"]),
            ("GFM label summary written (AOIs selected for amendment 3)", (gs.get("built_utc") or "")[:16] + "Z"),
            ("Second LLM-agent critiques of the manuscript recorded", blind2 or ""),
            ("Amendment 6 (post hoc analyses after the second critiques)", out["AmendTimeBSix"]),
            ("Addendum 1 to amendment 6", addendum or missing("addendum 1 of amendment 6")),
            ("Addendum 2 to amendment 6 (implementation record)", addenda23[0] or missing("addendum 2 of amendment 6")),
            ("Addendum 3 to amendment 6 (correction of the implementation record)",
             addenda23[1] or missing("addendum 3 of amendment 6"))]
    rows.sort(key=lambda r: r[1])
    body = "\n".join(f"{a} & {b} \\\\" for a, b in rows)
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table}[t]\n"
           "\\caption{Timeline of the study (UTC), in order. Times are those of the version-controlled repository, of "
           "the run records and of the container image (whose tag is its build time). Amendment $k$ is the file "
           "B-$k$ of the protocol's amendments in the code repository.}\\label{tab:timeline}\n\\small\n"
           "\\begin{tabular}{lr}\n\\tophline\nEvent & Time \\\\\n\\middlehline\n" + body +
           "\n\\bottomhline\n\\end{tabular}\n\\end{table}\n")
    (HERE / "app_table_timeline.tex").write_text(tex, encoding="utf-8")
    return out


def gfm_items() -> dict[str, str]:
    """GFM の窓に入った画像（項目）の数と、層が読めず外した項目の数（B-3 追記 1）。"""
    n_items = n_skip = n_date = 0
    for f in sorted((DATA / "global" / "labels").glob("*/gfm_meta.json")):
        meta = json.loads(f.read_text(encoding="utf-8"))
        if meta.get("window"):
            n_date += 1
        for sc in meta.get("scenes", []):
            n_items += len(sc.get("items", []))
            n_skip += len(sc.get("skipped_unreadable", []))
    return {"GFMNWithDate": str(n_date), "GFMNItems": _n(n_items), "GFMNSkipped": str(n_skip)}


# %% 付録 C の新しい表（B-3・B-4・B-5）
def _tbl(label: str, caption: str, cols: str, head: str, rows: list[str]) -> str:
    return ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[t]\n"
            f"\\caption{{{caption}}}\\label{{{label}}}\n\\footnotesize\n\\begin{{tabular}}{{{cols}}}\n\\tophline\n"
            f"{head} \\\\\n\\middlehline\n" + "\n".join(rows) + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")


def _orc(t: dict | None) -> tuple[str, str]:
    if not t or t.get("odds_ratio") is None:
        return "---", "---"
    lo, hi = t["ci95_odds_ratio"]
    return orv(t["odds_ratio"]), f"{orv(lo)}--{orv(hi)}"


def _ppc(r: dict) -> tuple[str, str]:
    return f"${100 * r['estimate']:+.1f}$", f"${100 * r['ci95'][0]:+.1f}$ to ${100 * r['ci95'][1]:+.1f}$"


def _rr(c: dict | str | None, key90: str = "ci90_R") -> tuple[str, str, str]:
    if not isinstance(c, dict):
        return "---", "---", "---"
    lab = c["decision"]["label"] if isinstance(c.get("decision"), dict) else c.get("decision", "---")
    return g._r(c["R"]), g._ci(c.get(key90)), lab


def table_gfm(v41: dict | None, v45: dict | None, v47: dict | None) -> str:
    if not v41:
        return "\\missing{v41}"
    rows = []
    for key, lab in (("all", "all AOIs"), ("sensor_s1", "EMS image includes Sentinel-1"),
                     ("sensor_sar", "EMS image from other radar"), ("sensor_optical", "EMS image optical")):
        a = v41["A0"].get(key)
        if a:
            rows.append(f"Agreement, {lab} & both & {100 * a['agreement']:.1f}\\,\\% ($\\kappa={a['kappa']:.2f}$) & "
                        f"--- & {_n(a['n_cells'])} cells, {a['n_aoi']} AOIs \\\\")
    a1 = v41["A1_A2"]
    for key, lab in (("A1_near_minus_far", "EMS $-$ GFM, $\\leq500$\\,m minus $>1$\\,km (pct. points)"),
                     ("A2_wetside_le500m", "EMS $-$ GFM, water-side minus other names, $\\leq500$\\,m"),
                     ("A2_wetside_all", "EMS $-$ GFM, water-side minus other names, all cells")):
        e, c = _ppc(a1[key])
        rows.append(f"{lab} & both & {e} & {c} & {a1[key]['n_aoi']} AOIs \\\\")
    for lab, src in (("GFM", v41["odds_gfm"]), ("EMS", v41["odds_ems"])):
        e, c = _orc(src["bands_adjusted_M1"]["terms"]["band_gt1km"])
        rows.append(f"Odds ratio, $>1$\\,km against $\\leq500$\\,m, given terrain & {lab} & {e} & {c} & "
                    f"{src['n_aoi']} AOIs \\\\")
    rows.append("\\middlehline")
    for lab, src in (("GFM", v41["odds_gfm"]), ("EMS", v41["odds_ems"])):
        for key, term, what in (("wetside_adjusted_M1", "wetside:le500m", "water-side names, $\\leq500$\\,m"),
                                ("wetside_raw", "wetside:le500m", "water-side names, $\\leq500$\\,m, unadjusted"),
                                ("wetside_allcells_adjusted_M1", "wetside:all", "water-side names, all cells")):
            e, c = _orc(src[key]["terms"].get(term))
            rows.append(f"Odds ratio, {what} & {lab} & {e} & {c} & --- \\\\")
    for cls in ("water", "wet", "riverine", "crossing", "low", "high"):
        e, c = _orc(v41["odds_gfm"]["categories_adjusted_M1"]["terms"].get(f"{cls}:le500m"))
        rows.append(f"\\quad {g.CLASS_LABELS[cls]} names, $\\leq500$\\,m & GFM & {e} & {c} & --- \\\\")
    if v45:
        for cut, what in (("v_ge_0.9", "cells $\\geq90$\\,\\% assessable"), ("v_ge_1.0", "cells fully assessable")):
            for lab in ("gfm", "ems"):
                e, c = _orc(v45[cut][lab]["wetside_adjusted_M1"]["wetside:le500m"])
                rows.append(f"Odds ratio, water-side, $\\leq500$\\,m, {what}$^{{a}}$ & {lab.upper()} & {e} & {c} & "
                            f"{_n(v45[cut][lab]['n_cells'])} cells \\\\")
    for lab in ("gfm", "ems"):
        oc = v41.get(f"odds_{lab}_c2_M1")
        if oc:
            e, c = _orc(oc["wetside_adjusted_M1"]["terms"].get("wetside:le500m"))
            rows.append(f"Odds ratio, water-side, $\\leq500$\\,m, adjusted for the terrain model refitted on GFM "
                        f"& {lab.upper()} & {e} & {c} & --- \\\\")
    sens = v41.get("sensitivity_independent_sensor")
    if sens:
        for lab in ("gfm", "ems"):
            e, c = _orc(sens[f"odds_{lab}"]["wetside_adjusted_M1"]["terms"].get("wetside:le500m"))
            rows.append(f"Odds ratio, water-side, $\\leq500$\\,m, EMS image not Sentinel-1 & {lab.upper()} & {e} & {c} & "
                        f"{sens['n_aoi']} AOIs \\\\")
    rows.append("\\middlehline")
    for key, what in (("C1_gfm", "$R$, held-out predictions rescored"), ("C1_ems", "$R$, same cells"),
                      ("C1_gfm_lexfix", "$R$, rescored, corrected Norwegian reading"),
                      ("C1_ems_lexfix", "$R$, same cells, corrected Norwegian reading")):
        r, c, lab = _rr(v41.get(key))
        rows.append(f"{what} & {'GFM' if 'gfm' in key else 'EMS'} & {r} & {c}$^{{b}}$ & {lab} \\\\")
    c2 = (v47 or {}).get("gfm_c2")
    if isinstance(c2, dict) and c2["decision"] == "undefined-denominator":
        y = c2
        rows.append(f"Terrain increment $\\Delta$BS(M0$\\to$M1), every model refitted on GFM labels & GFM & "
                    f"{sci(y['yardstick'])} & {sci(y['ci95_yardstick'][0])} to {sci(y['ci95_yardstick'][1])} & "
                    f"$R$ undefined \\\\")
        rows.append(f"Names increment $\\Delta$BS(M1c$\\to$M3n), refitted & GFM & {sci(y['delta_bs'])} & "
                    f"{sci(y['ci95_delta'][0])} to {sci(y['ci95_delta'][1])} & --- \\\\")
    else:
        r, c, lab = _rr(c2 if isinstance(c2, dict) else None)
        rows.append(f"$R$, every model refitted on GFM labels & GFM & {r} & {c}$^{{b}}$ & {lab} \\\\")
    cap = ("Independent labels from GFM (post hoc, amendment B-3) on the \\GFMNAoi{} AOIs assessable on both label "
           "sets. Odds ratios come from logistic regressions with an intercept for each AOI and the logit of the "
           "held-out terrain prediction, with standard errors clustered by place; pct.\\ points are within-AOI "
           "contrasts. $^{a}$Added after the result above had been seen (B-3 addendum~2). $^{b}$90\\,\\% interval.")
    return _tbl("tab:gfm", cap, "p{7.6cm}lccl", "Quantity & Labels & Estimate & 95\\,\\% interval & Note", rows)


def table_b5(v42: dict | None, v47: dict | None) -> str:
    rows = []
    if v42:
        for key, what in (("inversion_area_wide", "effect on all ground nearest a water-side name"),
                          ("inversion_local_500m", "effect within 500\\,m of the settlement")):
            inv = v42[key]
            rows.append(f"Largest odds ratio consistent with the observed increment, {what} & "
                        f"{orv(inv['or_upper_q05'])} & --- & normal approximation "
                        f"{orv(inv['or_upper_q05_normal'])} \\\\")
        w = v42["within_area"]
        rows.append(f"Change in AUC within AOIs, median & ${w['auc_diff_median']:+.4f}$ & "
                    f"${w['auc_diff_iqr'][0]:+.4f}$ to ${w['auc_diff_iqr'][1]:+.4f}$ (IQR) & {w['n_aoi_auc']} AOIs \\\\")
        rows.append(f"Change in PR-AUC within AOIs, median & ${w['prauc_diff_median']:+.4f}$ & "
                    f"${w['prauc_diff_iqr'][0]:+.4f}$ to ${w['prauc_diff_iqr'][1]:+.4f}$ (IQR) & --- \\\\")
        rows.append(f"Change in log-likelihood per cell, AOI intercepts refitted & {sci(w['cond_ll_per_cell'])} & "
                    f"{sci(w['cond_ll_ci95'][0])} to {sci(w['cond_ll_ci95'][1])} & --- \\\\")
        rows.append(f"$R$ leaving out one place at a time & {g._r(w['loo_R_range'][0])} to "
                    f"{g._r(w['loo_R_range'][1])} & --- & range \\\\")
        s = v42["settler_excluded"]
        rows.append(f"$R$ without the Americas and Oceania & {g._r(s['R'])} & {g._ci(s['ci90_R'])}$^{{a}}$ & "
                    f"{s['n_aoi_dropped']} AOIs removed \\\\")
    v47 = v47 or {}
    pm = v47.get("permutations")
    if isinstance(pm, dict):
        rows.append(f"$R$ with names reassigned {pm['n_perms']} times & {g._r(pm['R_perm_mean'])} (mean) & "
                    f"{g._r(pm['R_perm_q025'])} to {g._r(pm['R_perm_q975'])}$^{{b}}$ & real names "
                    f"{g._r(pm['R_real'])}, {pm['percentile_of_real']:.0f}th percentile \\\\")
    fs = (v47.get("fold_seeds") or {}).get("fixed_code") or {}
    if "combined_ci90" in fs:
        rows.append(f"$R$ over {fs['n_available']} fold assignments & {g._r(fs['R_min'])} to {g._r(fs['R_max'])} & "
                    f"{g._ci(fs['combined_ci90'])}$^{{c}}$ & {fs['n_equivalent_null']} equivalent-null \\\\")
    oc = v47.get("oracle_and_coefs")
    if isinstance(oc, dict):
        o = oc["oracle"]
        rows.append(f"$R$ of the oracle (terrain at the settlement's site) & {g._r(o['R'])} & "
                    f"{g._ci(o['ci95_R'])}$^{{d}}$ & {o['decision']} \\\\")
    for key, what in (("notrace", "flood traces counted as not flooded"), ("nonwater", "flooded share over land")):
        lv = (v47.get("label_variants") or {}).get(key)
        if isinstance(lv, dict):
            rows.append(f"$R$, labels with {what} & {g._r(lv['R'])} & {g._ci(lv['ci90_R'])}$^{{a}}$ & "
                        f"{lv['decision']} \\\\")
    cap = ("Post hoc analyses of amendment B-5 on the primary split. The largest odds ratios are upper bounds by test "
           "inversion under the simulation model (the largest effect whose simulated increments place the observed "
           "increment above their 5th percentile). $^{a}$90\\,\\% interval. $^{b}$2.5th to 97.5th percentile of the "
           "reassignments. $^{c}$90\\,\\% interval combining the bootstrap error of the first assignment with the "
           "spread over assignments. $^{d}$95\\,\\% interval.")
    return _tbl("tab:b5", cap, "p{8.2cm}ccl", "Quantity & Estimate & Interval & Note", rows or ["\\missing{B-5}"])


def table_groups(v42: dict | None) -> str:
    if not v42:
        return "\\missing{v42}"
    rows = []
    names = {"wetland": "wetland", "island": "island and river island$^{a}$", "spring_well": "spring and well",
             "crossing": "crossing", "dike_reclamation": "dike and polder", "lowland": "lowland (benchmark)"}
    expect = {"wetland": "more flooding, lower", "island": "site higher", "spring_well": "unrelated to flooding",
              "crossing": "none", "dike_reclamation": "none", "lowland": "none"}  # amendment 5 で固定した向き
    for key, lab in names.items():
        gb = v42["groups"]["by_group"][key]
        e1, c1 = _orc(gb.get("or_le500m_adj"))
        e2, c2 = _orc(gb.get("or_all_adj"))
        hm = gb.get("pc_hand_mean_all")
        hand = (f"${hm['estimate']:+.1f}$ (${hm['ci95'][0]:+.1f}$ to ${hm['ci95'][1]:+.1f}$)"
                if isinstance(hm, dict) else "---")
        ph = gb.get("p_le500m_holm")
        rows.append(f"{lab} & {expect[key]} & {gb['n_settlements']} & {_n(gb['n_cells_le500m'])} & {gb['n_flooded_le500m']} & "
                    f"{e1} ({c1}) & {e2} ({c2}) & {hand} & {'---' if ph is None else f'{ph:.2f}'} \\\\")
    cap = ("Element groups of the lexicon (post hoc, amendment B-5), with the direction that the amendment expected "
           "for each group: odds ratios of flooding given the held-out terrain prediction within 500\\,m of the "
           "settlement and over all cells, the difference in mean HAND (m) between the cells of the group and other "
           "cells of the same AOI, and Holm-adjusted $p$ values. Each group is compared with all other cells of the "
           "same AOI; the lowland group is a benchmark, not the reference category of the regression. No cell within "
           "500\\,m of a dike or polder name was flooded, so that odds ratio cannot be estimated. $^{a}$Includes German "
           "\\emph{Au} and \\emph{Aue}, which mean water meadow rather than island.")
    return _tbl("tab:groups", cap, "llrrrcccc",
                "Group & Expected & Settlements & Cells $\\leq500$\\,m & Flooded & OR $\\leq500$\\,m & OR all cells & HAND & "
                "Holm $p$", rows)


def table_lexfix(v47: dict | None) -> str:
    v47 = v47 or {}
    rows = []
    for split, lab in (("primary", "places held out (primary)"), ("country", "countries held out"),
                       ("temporal", "time held out"), ("llmcut", "after the language model's cut-off")):
        a, b = (v47.get("frozen") or {}).get(split), (v47.get("lexfix") or {}).get(split)
        ra, ca, la = _rr(a, "ci95_R")
        rb, cb, lb = _rr(b, "ci95_R")
        rows.append(f"{lab} & {ra} & {ca} & {la} & {rb} & {cb} & {lb} \\\\")
    cap = ("The confirmatory analysis with the frozen code (lexicon not applied to Norwegian names) and with the "
           "corrected reading (amendment B-4): $R$ of the within-area comparison, 95\\,\\% interval and verdict.")
    return _tbl("tab:lexfix", cap, "lcclccl",
                "Split & $R$ (frozen) & 95\\,\\% & Verdict & $R$ (corrected) & 95\\,\\% & Verdict", rows)


CLS = ("water", "wet", "riverine", "crossing", "low", "high")
FEAT_LABEL = {**{f"lex_{c}_n1": f"{g.CLASS_LABELS[c]}, nearest name" for c in CLS},
              **{f"lex_{c}_k": f"{g.CLASS_LABELS[c]}, three nearest (weighted)" for c in CLS},
              "lex_wetside_n1": "water-side, nearest name", "lex_covered_n1": "lexicon covers the language",
              "nlearn": "character-model score"}


def table_coefs(v47: dict | None) -> str:
    """B-5 の 2: 二段目の係数（標準化した尺度、10 fold の平均と範囲）。"""
    oc = (v47 or {}).get("oracle_and_coefs")
    if not isinstance(oc, dict) or not oc.get("coefs"):
        return "\\missing{coefficients (world\\_v2\\_b5coef)}"
    cf = oc["coefs"]

    def cell(stack: str, name: str) -> str:
        c = (cf.get(stack) or {}).get(name)
        return "---" if c is None else f"${c['mean']:+.3f}$ (${c['min']:+.3f}$, ${c['max']:+.3f}$)"

    first = (f"logit of the terrain prediction & {cell('M3n', 'logit_M1')} & {cell('M3np', 'logit_M1')} & "
             f"{cell('M3d', 'logit_M1')} & {cell('M1c', 'logit_M1')} \\\\")
    rows = [first, "\\middlehline"]
    for f, lab in FEAT_LABEL.items():
        pf = f"p_{f}" if f != "nlearn" else "nlearn_p"
        rows.append(f"{lab} & {cell('M3n', f)} & {cell('M3np', pf)} & {cell('M3d', f'd_{f}')} & "
                    f"{cell('M1c', f'c_{f}')} \\\\")
    cap = ("Second-stage coefficients on the standardised scale (post hoc, amendment B-5): mean over the ten folds "
           "of the primary split, with the smallest and largest fold in parentheses. M3n: the name features; M3np: "
           "the same features of names reassigned among settlements; M3d: the features centred within each AOI; M1c: "
           "the regional composition (AOI means). The reading of Norwegian names is the corrected one (B-4).")
    return _tbl("tab:coefs", cap, "lcccc", "Feature & M3n & M3np (placebo) & M3d & M1c (composition)", rows)


LANG_NAMES = {"en": "English", "fr": "French", "de": "German", "pl": "Polish", "sl": "Slovene", "cs": "Czech",
              "hu": "Hungarian", "it": "Italian", "es": "Spanish", "pt": "Portuguese", "nl": "Dutch", "ro": "Romanian"}


def _kn(r: dict) -> str:
    return f"{r['k']}/{r['n']}"


def table_etym(v44: dict | None) -> str:
    if not v44:
        return "\\missing{v44}"
    rows = []
    for c, r in v44["by_lexicon_category"].items():
        rows.append(f"{g.CLASS_LABELS[c]} & {r['same_category']['n']} & "
                    f"{_rate(r['same_category'])} & {_rate(r['any_wetside'])} \\\\")
    rows += ["\\middlehline",
             "Language of the name & Coded & Water-related, read as water-side & Water-related, others \\\\",
             "\\middlehline"]
    for lang, r in v44["by_language"].items():
        rows.append(f"{LANG_NAMES.get(lang, lang)} & {r['n']} & "
                    f"{_rate(r['precision_lex'])}; {_kn(r['precision_lex'])} & "
                    f"{_rate(r['truth_rate_not_flagged'])}; {_kn(r['truth_rate_not_flagged'])} \\\\")
    cap = ("Recorded etymologies on Wikipedia by lexicon class and by language of the name (post hoc, amendment B-5). "
           "Upper part: among the names in each class of the lexicon whose Wikipedia section states a derivation, "
           "the share whose recorded derivation is in the same class and the share whose derivation is "
           "water-related in any class (the highland class is not water-side). Lower part, for languages with at "
           "least 30 coded sections: the share of water-related derivations among the names that the lexicon read "
           "as water-side (its precision) and among the other names, with counts. Shares are unweighted; intervals "
           "are 95\\,\\% (Wilson).")
    return _tbl("tab:etym-detail", cap, "lccc",
                "Lexicon class & Coded & Same class & Water-related, any class", rows)


def write_new_tables(v41: dict | None, v42: dict | None, v44: dict | None, v45: dict | None,
                     v47: dict | None) -> None:
    for name, text in (("gfm", table_gfm(v41, v45, v47)), ("b5", table_b5(v42, v47)),
                       ("groups", table_groups(v42)), ("lexfix", table_lexfix(v47)),
                       ("coefs", table_coefs(v47)), ("etym", table_etym(v44))):
        (HERE / f"app_table_{name}.tex").write_text(text, encoding="utf-8")


def write_numbers(m: dict[str, str]) -> None:
    lines = ["% 自動生成（figures/make_figures.py）。手で編集しない。",
             "\\providecommand{\\missing}[1]{\\textbf{[??\\,#1]}}"]
    lines += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(m.items())]
    (HERE / "numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


# %% B-6（初見の査読 第 2 回への事後の解析）のマクロ
def _ci2(prefix: str, est: float | None, ci: list | None, fmt=None) -> dict[str, str]:
    """値と区間（{prefix}, {prefix}Lo, {prefix}Hi）。fmt の既定は rv。"""
    f = fmt or rv
    if est is None or ci is None:
        return {f"{prefix}{s}": missing(prefix) for s in ("", "Lo", "Hi")}
    return {prefix: f(est), f"{prefix}Lo": f(ci[0]), f"{prefix}Hi": f(ci[1])}


def _pct0(x: float) -> str:
    return f"{100 * x:.0f}\\,\\%"


def _verdict(c: dict | None) -> str:
    return c["decision"] if isinstance(c, dict) and c.get("decision") else missing("verdict")


def _cmp_macros(prefix: str, c: dict | None) -> dict[str, str]:
    """B-6 の比較（v52 の verdict の形）: R、90%・95% 区間、判定。"""
    if not isinstance(c, dict) or c.get("R") is None:
        return {f"{prefix}{s}": missing(prefix) for s in ("R", "RLoNinety", "RHiNinety", "RLo", "RHi", "Verdict")}
    m = {f"{prefix}R": rv(c["R"]), f"{prefix}RLoNinety": rv(c["ci90_R"][0]), f"{prefix}RHiNinety": rv(c["ci90_R"][1]),
         f"{prefix}RLo": rv(c["ci95_R"][0]), f"{prefix}RHi": rv(c["ci95_R"][1]),
         f"{prefix}Verdict": c["decision"]}
    pl = c.get("placebo")
    if isinstance(pl, dict) and pl.get("R_placebo") is not None:  # 名前を置き換えた 1 回（placebo）の R
        m[f"{prefix}Rp"] = rv(pl["R_placebo"])
    return m


def actgrp_macros() -> dict[str, str]:
    """事象をまとめて保持した分割（B-2）の R の 95% 区間（not-meaning の判定に使う区間）。"""
    m = g._json(DATA / "runs" / "world_v2_actgrp" / "metrics_primary.json")
    if not m:
        return {"RActGrpLo": missing("actgrp"), "RActGrpHi": missing("actgrp")}
    return {"RActGrpLo": rv(m["primary"]["ci95_R"][0]), "RActGrpHi": rv(m["primary"]["ci95_R"][1])}


def ti_spread_macros(v42: dict | None, v47: dict | None) -> dict[str, str]:
    """B-5 の 4 の検定の反転が実データの不確かさを小さく見た理由（B-6 の 2 の報告の規則）: 合成の null の平均・標準偏差と、
    置換の null の平均・確証の増分の bootstrap の標準誤差（直したコード）。"""
    lex = g._json(DATA / "runs" / "world_v2_lexfix" / "metrics_primary.json")
    pm = (v47 or {}).get("permutations")
    if not (v42 and lex and isinstance(pm, dict)):
        return {k: missing("test-inversion spread") for k in ("TISimSd", "TISimMean", "TIRealSe", "TIPermMean")}
    z = v42["inversion_area_wide"]["grid"]["0.0"]
    c = lex["primary"]["ci95_delta"]
    return {"TISimSd": sci(z["sd"]), "TISimMean": sci(z["mean"]), "TIRealSe": sci((c[1] - c[0]) / 3.92),
            "TIPermMean": sci(pm["R_perm_mean"] * pm["yardstick"])}


def fs_frozen_macros(v47: dict | None) -> dict[str, str]:
    """B-5 追記 3: 凍結したコードの fold の割り当て（修正前に回した 3 つ）と、修正による R の動き（同じ割り当て）。"""
    fs = (v47 or {}).get("fold_seeds") or {}
    fz, fx = (fs.get("frozen_code") or {}).get("R") or {}, (fs.get("fixed_code") or {}).get("R") or {}
    if not fz:
        return {k: missing("fold seeds") for k in ("FSFrozenMin", "FSFrozenMax", "FSFixShiftMax")}
    shift = max(abs(fz[k] - fx[k]) for k in fz if k in fx)
    return {"FSFrozenMin": rv(min(fz.values())), "FSFrozenMax": rv(max(fz.values())), "FSFixShiftMax": rv(shift),
            "FSFrozenN": _word(len(fz))}


def b6_ceiling_macros(v48: dict | None) -> dict[str, str]:
    m: dict[str, str] = {}
    if not v48:
        return {k: missing("v48") for k in ("CeilRFiveHundred", "CeilRTwoFifty", "CeilROneKm", "CeilRAll",
                                            "CeilShareFiveHundred", "AoiWR", "AoiWRLoNinety", "AoiWRHiNinety",
                                            "RoundedSharePct")}
    c = v48["ceiling"]["lexfix"]
    m |= {"CeilRFiveHundred": rv(c["le500m"]["R_max"]), "CeilRTwoFifty": rv(c["le250m"]["R_max"]),
          "CeilROneKm": rv(c["le1km"]["R_max"]), "CeilRAll": rv(c["all_wetside"]["R_max"]),
          "CeilShareFiveHundred": pct(c["le500m"]["share_cells"]),
          "CeilFloodedFiveHundred": pct(c["le500m"]["flooded_share"]),
          "CeilRFiveHundredFrozen": rv(v48["ceiling"]["frozen"]["le500m"]["R_max"])}
    a = v48["aoi_weighted"]["frozen"]
    m |= {"AoiWR": rv(a["R"]), "AoiWRLoNinety": rv(a["ci90_R"][0]), "AoiWRHiNinety": rv(a["ci90_R"][1]),
          "AoiWRLo": rv(a["ci95_R"][0]), "AoiWRHi": rv(a["ci95_R"][1]),
          "AoiWVerdictEquiv": "held" if a["equivalent_TOST"] else "did not hold"}
    r = v48["rounded"]
    m["RoundedSharePct"] = pct(r["share_cells_minute_rounded"])
    m |= or_macros("RoundEms", r["ems_all"]["not_rounded"]["le500m"])
    m |= or_macros("RoundGfm", r["gfm_cells"]["gfm"]["not_rounded"]["le500m"])
    m |= or_macros("RoundEmsGfmCells", r["gfm_cells"]["ems"]["not_rounded"]["le500m"])
    cov = v48["coverage"]
    nc = cov.get("null_coverage") or {}
    dc = cov.get("null_decision_counts") or {}
    m |= {"CovNinety": _pct0(nc["coverage90"]) if nc else missing("coverage"),
          "CovNinetyFive": _pct0(nc["coverage95"]) if nc else missing("coverage"),
          "CovNinetyZero": _pct0(nc["coverage90_zero"]) if nc else missing("coverage"),
          "CovNinetyFiveZero": _pct0(nc["coverage95_zero"]) if nc else missing("coverage"),
          "CovSdRatio": f"{nc['sd_R_between_reps'] / nc['median_bootstrap_se_R']:.1f}" if nc else missing("coverage"),
          "NullNEquiv": str(dc.get("equivalent-null", 0)), "NullNNegative": str(dc.get("negative", 0)),
          "NullNInconclusive": str(dc.get("inconclusive", 0)), "NullNNotMeaning": str(dc.get("not-meaning", 0)),
          "NullNReps": str(cov.get("n_reps") or missing("null reps"))}
    return m


def b6_gfm_macros(v50: dict | None) -> dict[str, str]:
    if not v50:
        return {"GFMCSI": missing("v50")}
    ag = v50["agreement"]
    m = {"GFMCSI": f"{ag['all']['csi']:.2f}", "GFMCSIOptical": f"{ag['sensor_optical']['csi']:.2f}",
         "GFMCSISOne": f"{ag['sensor_s1']['csi']:.2f}", "GFMCSISar": f"{ag['sensor_sar']['csi']:.2f}"}
    sel = v50["aoi_vs_cell_selection"]
    m |= or_macros("SelEmsAll", sel["ems_192aoi_all_cells"]) | or_macros("SelEmsGfm", sel["ems_gfm_cells"])
    m |= or_macros("SelGfm", sel["gfm_gfm_cells"])
    for key, pre in (("gfm_gfm_cells", "CbGfm"), ("ems_gfm_cells", "CbEmsGfm"), ("ems_192aoi_all_cells", "CbEmsSel"),
                     ("ems_all_aoi", "CbEmsAll")):
        c = v50["contributing"][key]
        m |= {f"{pre}Flooded": _n(c["n_flooded_group"]), f"{pre}Cells": _n(c["n_cells_group"]),
              f"{pre}Aoi": str(c["n_aoi_with_flooded"]), f"{pre}Units": str(c["n_units_with_flooded"])}
    for key, pre in (("gfm_gfm_cells", "RbGfm"), ("ems_gfm_cells", "RbEmsGfm"), ("ems_all_aoi", "RbEmsAll")):
        r = v50["robustness"][key]
        m |= {f"{pre}PairsLo": orv(r["pairs_bootstrap"]["ci95_odds_ratio"][0]),
              f"{pre}PairsHi": orv(r["pairs_bootstrap"]["ci95_odds_ratio"][1]),
              f"{pre}LooMin": orv(r["leave_one_place_out"]["min_odds_ratio"]),
              f"{pre}LooMax": orv(r["leave_one_place_out"]["max_odds_ratio"])}
    c2 = v50.get("c2") or {}
    if c2:
        b = c2["brier"]
        m |= {"CTwoBSConst": f"{b['const_train_rate']:.4f}", "CTwoBSConstEval": f"{b['const_eval_rate_pooled']:.4f}",
              "CTwoBSMZero": f"{b['M0']:.4f}", "CTwoBSMOne": f"{b['M1']:.4f}", "CTwoBSMOnes": f"{b['M1s']:.4f}",
              "CTwoBSMOnec": f"{b['M1c']:.4f}", "CTwoBSMThreed": f"{b['M3d']:.4f}"}
        m |= _ci2("CTwoDm", c2["delta_bs_M1s_M3d"], c2["ci95_delta_M1s_M3d"], sci)
    v1 = v50.get("v1")
    if v1:
        m |= or_macros("VOneOR", v1["odds_v1"]["wetside_adjusted_M1"]["terms"]["wetside:le500m"])
        m |= or_macros("VOneOREms", v1["odds_ems"]["wetside_adjusted_M1"]["terms"]["wetside:le500m"])
        m |= or_macros("VOneORAll", v1["odds_v1"]["wetside_allcells_adjusted_M1"]["terms"]["wetside:all"])
        m |= {"VOneNAoi": str(v1["n_aoi"]), "VOneNCells": _n(v1["n_cells"]), "VOneRate": pct(v1["rate_v1"]),
              "VOneCSI": f"{v1['agreement']['csi']:.2f}"}
        c1 = v1["C1_v1"]
        m |= {"VOneCOneR": rv(c1["R"]) if c1.get("R") is not None else missing("V1 C1"),
              "VOneCOneRLoNinety": rv(c1["ci90_R"][0]) if c1.get("R") is not None else missing("V1 C1"),
              "VOneCOneRHiNinety": rv(c1["ci90_R"][1]) if c1.get("R") is not None else missing("V1 C1"),
              "VOneCOneVerdict": (c1["decision"]["label"] if isinstance(c1.get("decision"), dict)
                                  else str(c1.get("decision")))}
        cb = v1["contributing_v1"]
        m |= {"VOneCbFlooded": _n(cb["n_flooded_group"]), "VOneCbUnits": str(cb["n_units_with_flooded"])}
        rb = (v50.get("robustness_v1") or {}).get("v1_v1_cells")
        if rb:
            m |= {"VOnePairsLo": orv(rb["pairs_bootstrap"]["ci95_odds_ratio"][0]),
                  "VOnePairsHi": orv(rb["pairs_bootstrap"]["ci95_odds_ratio"][1])}
    else:
        m |= {k: missing("GFM V1") for k in ("VOneOR", "VOneORLo", "VOneORHi", "VOneOREms", "VOneOREmsLo",
                                             "VOneOREmsHi", "VOneNAoi", "VOneCOneR", "VOneCOneVerdict")}
    rep = v50.get("v1_reproduces_b3") or {}
    if "n_aoi_identical" in rep:
        m |= {"VOneRepAoi": str(rep["n_aoi_identical"]), "VOneRepNAoi": str(rep["n_aoi"]),
              "VOneRepCells": _n(rep["n_cells_all"]), "VOneRepDiff": str(rep["n_eval_differs"] + rep["n_label_differs"]),
              "VOneRepItems": str(rep["n_items_unreadable_on_rebuild"])}
    comp = v50.get("components")
    if comp:
        a, w, o = comp["all"], comp["wetside_le500m"], comp["other_le500m"]
        m |= {"MaskShareAll": pct(a["share_not_gfm_evaluable"]), "MaskShareWet": pct(w["share_not_gfm_evaluable"]),
              "MaskShareOther": pct(o["share_not_gfm_evaluable"]),
              "MaskEmsRateNotEval": pct(a["ems_rate_not_evaluable"]), "MaskEmsRateEval": pct(a["ems_rate_evaluable"]),
              "MaskExclShare": pct(comp["share_not_evaluable_with_excl_majority"]),
              "MaskExclHandHigh": pct(comp["excluded_majority"]["share_hand_ge10m"]),
              "MaskEvalHandHigh": pct(comp["gfm_evaluable"]["share_hand_ge10m"]),
              "MaskExclBuilt": pct(comp["excluded_majority"]["share_built_ge05"]),
              "MaskEvalBuilt": pct(comp["gfm_evaluable"]["share_built_ge05"]),
              "MaskExclTree": pct(comp["excluded_majority"]["share_tree_ge05"]),
              "MaskEvalTree": pct(comp["gfm_evaluable"]["share_tree_ge05"]),
              "MaskSeasWet": pct(w["mean_components_not_evaluable"]["c_seas"]),
              "MaskSeasAll": pct(a["mean_components_not_evaluable"]["c_seas"]),
              "MaskExclAll": pct(a["mean_components_not_evaluable"]["c_excl"]),
              "MaskNoDataAll": pct(a["mean_components_not_evaluable"]["c_nodata"]),
              "MaskPermAll": pct(a["mean_components_not_evaluable"]["c_perm"]),
              "MaskSeasOther": pct(o["mean_components_not_evaluable"]["c_seas"])}
    req = [*(f"VOne{k}" for k in ("OR", "ORLo", "ORHi", "OREms", "OREmsLo", "OREmsHi", "ORAll", "ORAllLo", "ORAllHi",
                                  "NAoi", "NCells", "Rate", "CSI", "COneR", "COneRLoNinety", "COneRHiNinety",
                                  "COneVerdict", "CbFlooded", "CbUnits", "PairsLo", "PairsHi")),
           *(f"Mask{k}" for k in ("ShareAll", "ShareWet", "ShareOther", "EmsRateNotEval", "EmsRateEval", "ExclShare",
                                  "ExclHandHigh", "EvalHandHigh", "ExclBuilt", "EvalBuilt", "ExclTree", "EvalTree",
                                  "SeasWet", "SeasOther"))]
    m |= {k: missing("GFM V1 / components") for k in req if k not in m}
    return m


def b6_readings_macros(v51: dict | None) -> dict[str, str]:
    if not v51:
        return {"MatchAnyPct": missing("v51")}
    mt = v51["match"]
    m = {"LexNCodes": str(v51["lexicon"]["n_codes"]),
         "LexNForms": _n(sum(v51["lexicon"]["n_forms_by_category"].values())),
         "MatchAnyPct": pct(mt["share_cells_any"]), "MatchWetPct": pct(mt["share_cells_wetside"])}
    pb = v51["precision_by_position"]
    for pos, pre in (("word", "PrecWord"), ("prefix", "PrecPrefix"), ("suffix", "PrecSuffix")):
        r = pb[pos]["any"]
        m |= {pre: _pct0(r["rate"]) if r["rate"] is not None else "--", f"{pre}K": str(r["k"]), f"{pre}N": str(r["n"])}
    st = v51["strict"]
    m |= {"StrictPrec": _pct0(st["precision"]["est"]), "StrictPrecLo": _pct0(st["precision"]["ci95"][0]),
          "StrictPrecHi": _pct0(st["precision"]["ci95"][1]), "StrictRecall": _pct0(st["recall"]["est"]),
          "StrictRecallLo": _pct0(st["recall"]["ci95"][0]), "StrictRecallHi": _pct0(st["recall"]["ci95"][1]),
          "StrictNSpringOnly": str(st["n_water_spring_only"]),
          "StrictPrecCell": f"{100 * st['precision']['est']:.0f}\\,\\% ({100 * st['precision']['ci95'][0]:.0f}--"
                            f"{100 * st['precision']['ci95'][1]:.0f})",
          "StrictRecallCell": f"{100 * st['recall']['est']:.0f}\\,\\% ({100 * st['recall']['ci95'][0]:.0f}--"
                              f"{100 * st['recall']['ci95'][1]:.0f})"}
    for c, pre in (("water", "Water"), ("wet", "Wet"), ("riverine", "Riverine"), ("crossing", "Crossing"),
                   ("low", "Low"), ("high", "High")):
        r = st["by_description_category"][c]["precision"]
        m |= {f"CatPrec{pre}K": str(r["k"]), f"CatPrec{pre}N": str(r["n"])}
    m["AuditN"] = str(v51["audit"]["n"])
    return m


def _or_tex(x: dict | None) -> str | None:
    """plasmode のオッズ比の書式: 格子の外は「> 20」（数式の >）、10 以上は小数 1 桁、それ未満は 2 桁。"""
    if not x:
        return None
    if x.get("odds_ratio") is None:
        return r"\ensuremath{>}\,20" if x.get("how") == "above grid" else "undefined"
    v = float(x["odds_ratio"])
    return f"{v:.1f}" if v >= 10 else f"{v:.2f}"


def b6_cluster_macros(v52: dict | None) -> dict[str, str]:
    m: dict[str, str] = {}
    v52 = v52 or {}

    def blk(key: str) -> dict:
        x = v52.get(key)
        return x if isinstance(x, dict) else {}

    m |= _cmp_macros("Strict", blk("strict").get("strict"))
    m |= _cmp_macros("TZeroX", blk("t0x").get("T0x")) | _cmp_macros("TZeroB", blk("t0x").get("T0"))
    m |= _cmp_macros("Gn", blk("gn").get("primary")) | _cmp_macros("NoHydro", blk("nohydro").get("primary"))
    m |= _cmp_macros("Typed", blk("typed").get("primary")) | _cmp_macros("NReadPrim", blk("nread").get("H4"))
    m |= _cmp_macros("LabOracle", blk("label_oracle").get("oracle_label"))
    m |= _cmp_macros("GswOracle", blk("gsw_oracle").get("oracle_gsw"))
    m |= _cmp_macros("GswNpOracle", blk("gsw_oracle").get("oracle_gsw_nonperm"))
    gz = (blk("gn").get("b6_record") or {}).get("gazetteer")
    m |= ({"GnNP": _n(gz["n_p"]), "GnNExclCode": _n(gz["n_excluded_fcode"]), "GnNDup": _n(gz["n_duplicates"])}
          if gz else {k: missing("gazetteer filter") for k in ("GnNP", "GnNExclCode", "GnNDup")})
    n_all = blk("strict").get("n_aoi")  # 全 AOI（同じ表の実行）
    for key, pre in (("gn", "Gn"), ("nohydro", "NoHydro"), ("typed", "Typed")):
        b = blk(key)
        m[f"{pre}NAoi"] = str(b["n_aoi"]) if b else missing(key)
        m[f"{pre}NCells"] = _n(b["n_cells"]) if b else missing(key)
        m[f"{pre}NExcl"] = str(n_all - b["n_aoi"]) if b and n_all else missing(key)
    ap = v52.get("actperm")
    act_keys = ("ActPermPct", "ActPermN", "ActPermPlacPct", "ActPermRReal", "ActPermRMean", "ActPermRLo", "ActPermRHi")
    if isinstance(ap, dict) and ap.get("yardstick"):
        y = float(ap["yardstick"])
        m |= {"ActPermPct": f"{ap['percentile_of_real']:.0f}", "ActPermN": str(ap["n_perms"]),
              "ActPermPlacPct": f"{ap['percentile_of_single_placebo']:.0f}",
              # 直したコードの実の増分と、置換の平均・95% の範囲（R の尺度、物差しは B-2 の同じ分割）
              "ActPermRReal": rv(ap["R_real"]), "ActPermRMean": rv(ap["R_perm_mean"]),
              "ActPermRLo": rv(ap["delta_bs_perm_q025"] / y), "ActPermRHi": rv(ap["delta_bs_perm_q975"] / y)}
    else:
        m |= {k: missing("actperm") for k in act_keys}
    pl = v52.get("plasmode")
    forms = {"A": "A", "B": "B", "CA": "CA", "CB": "CB"}
    if isinstance(pl, dict):
        m |= {"PlasReps": str(pl["reps"]), "PlasNullDelta": sci(pl["delta_bs_null_perm_mean"]),
              "PlasUNinety": sci(pl["U90"])}
        for f, pre in forms.items():
            r = pl["forms"][f]
            for key, mac in (("delta_A", "Bound"), ("delta_B", "BoundZero"), ("delta_A_U95", "BoundNinetyFive"),
                             ("delta_M", "Margin"), ("delta_M_prime", "MarginNull"), ("mde80", "Mde")):
                x = r.get(key)
                m[f"Plas{mac}{pre}"] = _or_tex(x) or missing(f"plasmode {f} {key}")
    else:
        m |= {k: missing("plasmode") for k in ("PlasReps", "PlasNullDelta", "PlasUNinety")}
        for pre in forms.values():
            for mac in ("Bound", "BoundZero", "BoundNinetyFive", "Margin", "MarginNull", "Mde"):
                m[f"Plas{mac}{pre}"] = missing("plasmode")
    return m


def llmcut_macros() -> dict[str, str]:
    """言語モデルの分割の位置の模型と、最良の定数の予測（評価側の浸水率 p で p(1 − p)）の Brier。"""
    m = g._json(DATA / "runs" / "world_v2" / "metrics_llmcut.json")
    if not m:
        return {"LlmcutBSMZero": missing("llmcut metrics"), "LlmcutBSConst": missing("llmcut metrics")}
    w = m["wet_rate"]
    return {"LlmcutBSMZero": f"{m['brier']['M0']:.4f}", "LlmcutBSConst": f"{w * (1 - w):.4f}"}


def b6_gsw_macros(v53: dict | None) -> dict[str, str]:
    if not v53:
        return {"GswOccWet": missing("v53")}
    m: dict[str, str] = {}
    for g_, pre in (("wetside", "Wet"), ("strict", "Strict"), ("riverine", "Riverine"), ("crossing", "Crossing"),
                    ("wet", "Wetland"), ("water", "Water"), ("high", "High")):
        for v, vp in (("occ", "Occ"), ("seas", "Seas")):
            c = v53["all_named"][g_][v]
            fmt = (lambda x: f"\\ensuremath{{{x:+.1f}}}".replace("+", "{+}").replace("-", "{-}")) if v == "occ" else \
                (lambda x: f"\\ensuremath{{{x:+.2f}}}".replace("+", "{+}").replace("-", "{-}"))
            m |= _ci2(f"Gsw{vp}{pre}", c["estimate"], c["ci95"], fmt) if c else {f"Gsw{vp}{pre}": missing("gsw")}
    m["GswNSettlements"] = _n(v53["info"]["n_settlements"])
    return m


def b6_desc_macros(v54: dict | None) -> dict[str, str]:
    if not v54:
        return {"RefCSIHand": missing("v54")}
    rp = v54["reference_points"]
    m = {"RefCSIHand": f"{rp['hand_threshold_map']['csi']:.2f}", "RefCSIMOne": f"{rp['M1_map']['csi']:.2f}",
         "RefCSIMThreen": f"{rp['M3n_map']['csi']:.2f}", "RefBSConst": f"{rp['brier']['const_train_rate']:.4f}",
         "RefBSMZero": f"{rp['brier']['M0']:.4f}", "RefBSMOne": f"{rp['brier']['M1']:.4f}"}
    rc = v54["residual_correlogram"]["corr"]
    for i, nm in enumerate(("A", "B", "C", "D", "E", "F")):
        m[f"ResCorr{nm}"] = f"{rc[i]:.2f}"
    rl = v54["regional_languages"]
    reg = sum(v["n_cells"] for v in rl["regions"].values())
    m |= {"RegLangPct": pct(reg / rl["n_cells"]), "AuPct": pct(rl["countries"]["AU"]["share"]),
          "IePct": pct(rl["countries"]["IE"]["share"]), "RiversidePct": pct(v54["riverside"]["share_cells"])}
    nr = v54["nread"]
    by = nr["wetside_by_continent_cells"]
    m |= {"NReadWetPct": pct(nr["wetside_all_cells"]["llm"]), "LexWetPct": pct(nr["wetside_all_cells"]["lexicon"]),
          "NReadWetEU": pct(by["EU"]["llm"]), "LexWetEU": pct(by["EU"]["lexicon"]),
          "NReadWetAF": pct(by["AF"]["llm"]), "LexWetAF": pct(by["AF"]["lexicon"]),
          "NReadUnknownPct": pct(nr["category_share_cells"]["unknown"]),
          "NReadNonEnvPct": pct(nr["category_share_cells"]["non_environmental"])}
    fl = v54["inclusion_flow"]
    m |= {"NAoiLabelledFlow": str(fl["n_aoi_with_labels"]), "NExclEval": str(fl["reasons"].get("eval<400", 0)),
          "NExclWet": str(fl["reasons"].get("wet<20", 0)),
          "NActLabelled": str(fl["included"]["n_activations"] + fl["activations_only_in_excluded"]),
          "NExclOtherType": str(fl["excluded"]["type"].get("other", 0))}
    pv = v54["p_values"]["world_v2"]
    tost = pv["H2_tost"]
    m["HTwoTostP"] = ("< " + f"{tost['resolution']:.5f}".rstrip("0")) if tost["p_tost"] == 0 else f"{tost['p_tost']:.4f}"
    if "H4_tost" in pv:
        t4 = pv["H4_tost"]
        m["HFourTostP"] = ("< " + f"{t4['resolution']:.5f}".rstrip("0")) if t4["p_tost"] == 0 else f"{t4['p_tost']:.4f}"
        m["HFourImproveP"] = f"{pv['H4_improvement_p']:.2f}"
    m |= {"HThreeP": f"{pv['H3']['p_delta0_gt_delta4']:.2f}", "HThreeSpearman": f"{pv['H3']['spearman']:.1f}",
          "HThreeSpearmanP": f"{pv['H3']['p_spearman_exact']:.2f}",
          "HolmHThree": f"{pv['holm']['H3']:.2f}", "HolmHFour": f"{pv['holm']['H4']:.2f}"}
    return m


# %% B-6 の補足の表と、付録の GFM の短い表（数値は JSON から。手打ちしない）
CONT_NAMES = {"EU": "Europe", "AS": "Asia", "AF": "Africa", "OC": "Oceania", "NA": "North America",
              "SA": "South America"}


def _tex(sx: str) -> str:
    return (str(sx).replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%").replace("_", "\\_")
            .replace("#", "\\#").replace("$", "\\$"))


def _r3(x: float | None) -> str:
    return "---" if x is None else rv(x)


def _ci90(c: dict | None) -> str:
    return "---" if not isinstance(c, dict) or c.get("ci90_R") is None else f"{rv(c['ci90_R'][0])} to {rv(c['ci90_R'][1])}"


def table_gfm_compact(v41: dict | None, v47: dict | None, v50: dict | None) -> str:
    """付録 C の短い GFM の表（B-3 と B-6 の 4 の主な量）。詳しい表は補足。"""
    rows: list[str] = []
    if v41:
        a = v41["A0"]["all"]
        rows.append(f"Agreement with EMS & {100 * a['agreement']:.1f}\\,\\% ($\\kappa={a['kappa']:.2f}$) & --- & "
                    f"{_n(a['n_cells'])} cells, {a['n_aoi']} AOIs \\\\")
    if v50:
        rows.append(f"Critical success index for flooded cells & {v50['agreement']['all']['csi']:.2f} & --- & "
                    f"S1 {v50['agreement']['sensor_s1']['csi']:.2f}, other radar {v50['agreement']['sensor_sar']['csi']:.2f}, "
                    f"optical {v50['agreement']['sensor_optical']['csi']:.2f} \\\\")
        c2 = v50.get("c2") or {}
        if c2:
            b = c2["brier"]
            rows.append("Refitted on GFM: Brier score, constant / M0 / M1 / M1s / M1c & "
                        f"{b['const_train_rate']:.4f} / {b['M0']:.4f} / {b['M1']:.4f} / {b['M1s']:.4f} / {b['M1c']:.4f} "
                        "& --- & $R$ undefined \\\\")
            rows.append(f"Refitted on GFM: $\\Delta$BS(M1s$\\to$M3d) & {sci(c2['delta_bs_M1s_M3d'])} & "
                        f"{sci(c2['ci95_delta_M1s_M3d'][0])} to {sci(c2['ci95_delta_M1s_M3d'][1])} & no composition \\\\")
    if v41:
        c = v41["C1_gfm"]
        rows.append(f"Rescored held-out predictions, $R$ & {rv(c['R'])} & {rv(c['ci90_R'][0])} to {rv(c['ci90_R'][1])}$^{{a}}$ "
                    f"& {c['decision']['label']} \\\\")
    rows.append("\\middlehline")
    chain = []
    v39 = _json("v39_within_aoi_odds.json")
    if v39:
        chain.append(("EMS, all analysed AOIs", (v39.get("wetside_adjusted_M1") or {}).get("terms", {}).get("wetside:le500m"),
                      (v50 or {}).get("robustness", {}).get("ems_all_aoi")))
    if v50:
        sel = v50["aoi_vs_cell_selection"]
        chain += [("EMS, all cells of the GFM AOIs", sel["ems_192aoi_all_cells"], None),
                  ("EMS, GFM-assessable cells", sel["ems_gfm_cells"], v50["robustness"].get("ems_gfm_cells")),
                  ("GFM, same cells", sel["gfm_gfm_cells"], v50["robustness"].get("gfm_gfm_cells"))]
        v1 = v50.get("v1")
        if v1:
            chain += [("GFM V1, seasonal water assessable", v1["odds_v1"]["wetside_adjusted_M1"]["terms"]["wetside:le500m"],
                       (v50.get("robustness_v1") or {}).get("v1_v1_cells")),
                      ("EMS, cells assessable in V1", v1["odds_ems"]["wetside_adjusted_M1"]["terms"]["wetside:le500m"],
                       (v50.get("robustness_v1") or {}).get("ems_v1_cells"))]
    for lab, t, rb in chain:
        o, ci = _orc(t)
        extra = "---"
        if rb:
            pb, lo = rb["pairs_bootstrap"], rb["leave_one_place_out"]
            extra = (f"pairs {orv(pb['ci95_odds_ratio'][0])}--{orv(pb['ci95_odds_ratio'][1])}; "
                     f"LOO {orv(lo['min_odds_ratio'])}--{orv(lo['max_odds_ratio'])}")
        rows.append(f"Odds ratio, water-side names, $\\leq$500\\,m: {lab} & {o} & {ci} & {extra} \\\\")
    return _tbl("tab:gfm", "Independent labels from GFM (post hoc; amendments 3 and 6), condensed; the full set of "
                "quantities is in the Supplement. Odds ratios: logistic regression with an intercept for each AOI and the "
                "logit of the held-out terrain prediction, standard errors clustered by place (95\\,\\% interval); "
                "pairs: pairs cluster bootstrap over places (1000 replicates); LOO: range when leaving out one place. "
                "$^{a}$90\\,\\% interval.",
                "p{5.8cm}p{3.2cm}p{2.9cm}p{3.8cm}", "Quantity & Estimate & Interval & Note", rows)


def b6_supp_tables(v48: dict | None, v50: dict | None, v51: dict | None, v52: dict | None, v53: dict | None,
                   v54: dict | None) -> None:
    """B-6 の補足の表（app_table_b6_*.tex）。"""
    out: dict[str, str] = {}
    # 辞書の言語の符号ごと（longtable）
    lt = REPO / "docs" / "verification" / "v51_lexicon_table.csv"
    if lt.is_file():
        d = pd.read_csv(lt)
        body = "\n".join(
            f"{r.code} & {_tex(r.languages)} & {_tex(r.countries)} & {r.n_water} & {r.n_wet} & {r.n_riverine} & "
            f"{r.n_crossing} & {r.n_low} & {r.n_high} & {_tex(r.positions)} & {r.n_spring_forms_strict} \\\\"
            for r in d.itertuples())
        out["b6_lexicon"] = (
            "% generated by figures/make_figures.py — do not edit\n{\\scriptsize\n"
            "\\begin{longtable}{lp{2.4cm}p{5.6cm}rrrrrrp{1.4cm}r}\n"
            "\\caption{The multilingual lexicon by block: languages covered, countries mapped to the block, number of "
            "distinct forms per class (water, wetland, riverine, crossing, lowland, highland), match positions used and "
            "number of spring and well forms removed in the strict reading. Countries are mapped to the blocks of their "
            "main languages; regional languages are not read.}\\label{tab:b6-lexicon}\\\\\n\\tophline\n"
            "Block & Languages & Countries & Wa & We & Ri & Cr & Lo & Hi & Positions & Spr \\\\\n\\middlehline\n"
            "\\endfirsthead\n\\tophline\nBlock & Languages & Countries & Wa & We & Ri & Cr & Lo & Hi & Positions & Spr "
            "\\\\\n\\middlehline\n\\endhead\n" + body + "\n\\bottomhline\n\\end{longtable}\n}\n")
    tf = REPO / "docs" / "verification" / "v51_top_forms.csv"
    if tf.is_file():
        d = pd.read_csv(tf)
        rows = [f"{r.lang} & \\emph{{{_tex(r.form)}}} & {r.pos} & {g.CLASS_LABELS.get(r.cat, r.cat)} & {_n(int(r.n_cells))} & "
                f"{int(r.n_settlements)} \\\\" for r in d.itertuples()]
        out["b6_forms"] = _tbl("tab:b6-forms", "The 20 water-side forms that matched the nearest name of the most analysed "
                               "cells (amendment 6): language block, form, position, class, number of cells and of "
                               "settlements.", "llllrr", "Block & Form & Position & Class & Cells & Settlements", rows)
    ex = REPO / "docs" / "verification" / "v51_examples.csv"
    if ex.is_file():
        d = pd.read_csv(ex)
        rows = [f"{g.CLASS_LABELS.get(r.lexicon_category, r.lexicon_category)} & {r.kind} & {_tex(r.name)} ({r.country}) & "
                f"{_tex(str(r.derivation_codes).replace(';', ', '))} & {_tex(str(r.quote)[:150])}"
                f"{'\\ldots' if len(str(r.quote)) > 150 else ''} \\\\" for r in d.itertuples()]
        out["b6_examples"] = _tbl("tab:b6-examples", "Examples of names read as water-side by the lexicon, with the "
                                  "derivation recorded on Wikipedia (amendment 6): for each class, the two lowest-numbered "
                                  "coded names whose recorded derivation matches the class and the lowest-numbered one whose "
                                  "derivation does not. Quotations abridged; sources in the code repository.",
                                  "llp{3.0cm}p{2.0cm}p{7.6cm}", "Class & Agreement & Name & Recorded derivation & Quotation",
                                  rows)
    if v54:
        fl = v54["inclusion_flow"]
        rows = []
        for key, lab in (("included", "Analysed"), ("excluded", "Excluded, all")):
            x = fl[key]
            rows.append(f"{lab} & {x['n_aoi']} & {x['n_activations']} & {x['type'].get('riverine', 0)} / "
                        f"{x['type'].get('unspecified-as-riverine', 0)} / {x['type'].get('unspecified', 0)} / "
                        f"{x['type'].get('other', 0)} & {', '.join(f'{CONT_NAMES.get(k, k)} {v}' for k, v in x['continent'].items())}"
                        f" & {_n(int(x['median_n_eval']))} \\\\")
        for r, x in fl["excluded_by_reason"].items():
            lab = {"wet<20": "\\quad fewer than 20 flooded cells", "eval<400": "\\quad fewer than 400 evaluable cells"}.get(
                r, r)
            rows.append(f"{lab} & {x['n_aoi']} & {x['n_activations']} & {x['type'].get('riverine', 0)} / "
                        f"{x['type'].get('unspecified-as-riverine', 0)} / {x['type'].get('unspecified', 0)} / "
                        f"{x['type'].get('other', 0)} & {', '.join(f'{CONT_NAMES.get(k, k)} {v}' for k, v in x['continent'].items())}"
                        f" & {_n(int(x['median_n_eval']))} \\\\")
        out["b6_flow"] = _tbl("tab:b6-flow", f"Inclusion of the {fl['n_aoi_with_labels']} labelled AOIs (amendment 6). "
                              "Flood type of the AOI: riverine / untyped but inland (treated as riverine) / untyped / other "
                              "types. Activations are counted within each row and overlap between rows.",
                              "lrrlp{5.2cm}r", "AOIs & Number & Activations & Type & Continent & Median evaluable cells", rows)
    if v51 and v54:
        pb = v51["precision_by_position"]
        st = v51["strict"]
        rows = [(f"Precision, matched as a whole word & {pb['word']['any']['k']}/{pb['word']['any']['n']} & "
                 f"{_pct0(pb['word']['any']['rate'])} \\\\"),
                (f"Precision, matched as an ending & {pb['suffix']['any']['k']}/{pb['suffix']['any']['n']} & "
                 f"{_pct0(pb['suffix']['any']['rate'])} \\\\"),
                (f"Precision, matched as a beginning & {pb['prefix']['any']['k']}/{pb['prefix']['any']['n']} & "
                 f"{_pct0(pb['prefix']['any']['rate'])} \\\\")]
        for c in ("water", "wet", "riverine", "crossing", "low", "high"):
            r, rc = st["by_description_category"][c]["precision"], st["by_description_category"][c]["recall_unweighted"]
            rows.append(f"Class {g.CLASS_LABELS.get(c, c)}: precision / recall (unweighted) & {r['k']}/{r['n']} / "
                        f"{rc['k']}/{rc['n']} & {_pct0(r['rate'])} / {_pct0(rc['rate'])} \\\\")
        rows.append(f"Strict reading: precision / recall (weighted) & --- & {_pct0(st['precision']['est'])} / "
                    f"{_pct0(st['recall']['est'])} \\\\")
        rows.append(f"Water codes reclassified as springs or wells only & {st['n_water_spring_only']} & --- \\\\")
        rows.append("\\middlehline")
        nr = v54["nread"]["wetside_by_continent_cells"]
        for k, x in nr.items():
            rows.append(f"Cells whose nearest name is water-side, {CONT_NAMES.get(k, k)}: language model / lexicon & "
                        f"{_n(x['n_cells'])} cells & {pct(x['llm'])} / {pct(x['lexicon'])} \\\\")
        rl = v54["regional_languages"]
        for k, x in {**rl["regions"], **rl["countries"]}.items():
            if x["n_cells"]:
                rows.append(f"Cells nearest a settlement in {_tex(k)} & {_n(x['n_cells'])} & {pct(x['share'])} \\\\")
        rows.append(f"Cells whose nearest name has a riverside construction & --- & {pct(v54['riverside']['share_cells'])} \\\\")
        out["b6_readings"] = _tbl("tab:b6-readings", "The readings in more detail (amendment 6). Precision against "
                                  "derivations recorded on Wikipedia by the position of the matched water-side form and by "
                                  "class; recall by class is unweighted and therefore not comparable with the weighted "
                                  "recall reported in the paper. Lower part: readings by continent and languages that are "
                                  "not read.", "p{9.4cm}rr", "Quantity & Count & Share", rows)
    if v50:
        rows = []
        for k, lab in (("all", "all"), ("sensor_s1", "EMS image includes Sentinel-1"), ("sensor_sar", "other radar"),
                       ("sensor_optical", "optical")):
            a = v50["agreement"][k]
            rows.append(f"CSI / P(GFM=1$\\mid$EMS=0) / P(EMS=1$\\mid$GFM=0), {lab} & {a['csi']:.2f} / "
                        f"{100 * a['p_gfm1_given_ems0']:.1f}\\,\\% / {100 * a['p_ems1_given_gfm0']:.1f}\\,\\% \\\\")
        for k, lab in (("gfm_gfm_cells", "GFM labels"), ("ems_gfm_cells", "EMS labels, same cells"),
                       ("ems_192aoi_all_cells", "EMS labels, all cells of the GFM AOIs"), ("ems_all_aoi", "EMS labels, all AOIs")):
            c = v50["contributing"][k]
            rows.append(f"Flooded cells within 500\\,m of water-side names, {lab} & {_n(c['n_flooded_group'])} of "
                        f"{_n(c['n_cells_group'])} cells; {c['n_aoi_with_flooded']} AOIs, {c['n_units_with_flooded']} places \\\\")
        comp = v50.get("components")
        if comp:
            for k, lab in (("all", "all cells"), ("wetside_le500m", "within 500\\,m of water-side names"),
                           ("other_le500m", "within 500\\,m of other names")):
                x = comp[k]
                mc = x["mean_components_not_evaluable"]
                rows.append(f"Not assessable on GFM, {lab}: share; mean exclusion / permanent / seasonal / no image & "
                            f"{pct(x['share_not_gfm_evaluable'])}; {pct(mc['c_excl'])} / {pct(mc['c_perm'])} / "
                            f"{pct(mc['c_seas'])} / {pct(mc['c_nodata'])} \\\\")
            e, v = comp["excluded_majority"], comp["gfm_evaluable"]
            rows.append(f"Cells mostly excluded vs assessable: HAND $\\geq$10\\,m; built-up $\\geq$0.5; trees $\\geq$0.5; "
                        f"EMS flooded & {pct(e['share_hand_ge10m'])} vs {pct(v['share_hand_ge10m'])}; {pct(e['share_built_ge05'])} vs "
                        f"{pct(v['share_built_ge05'])}; {pct(e['share_tree_ge05'])} vs {pct(v['share_tree_ge05'])}; "
                        f"{pct(e['ems_rate'])} vs {pct(v['ems_rate'])} \\\\")
        v1 = v50.get("v1")
        if v1:
            for lab, blk in (("V1 labels", v1["odds_v1"]), ("EMS, V1 cells", v1["odds_ems"])):
                t1 = blk["wetside_allcells_adjusted_M1"]["terms"]["wetside:all"]
                o, ci = _orc(t1)
                rows.append(f"Odds ratio, water-side names, all cells, {lab} & {o} ({ci}) \\\\")
            c1 = v1["C1_v1"]
            if c1.get("R") is not None:
                rows.append(f"$R$, held-out predictions rescored on V1 labels & {rv(c1['R'])} (90\\,\\% {rv(c1['ci90_R'][0])} "
                            f"to {rv(c1['ci90_R'][1])}), {c1['decision']['label']} \\\\")
        if v48:
            r = v48["rounded"]
            for lab, blk in (("EMS, all AOIs", r["ems_all"]), ("GFM", r["gfm_cells"]["gfm"]), ("EMS, GFM cells", r["gfm_cells"]["ems"])):
                o1, c1_ = _orc(blk["all_cells"]["le500m"])
                o2, c2_ = _orc(blk["not_rounded"]["le500m"])
                rows.append(f"Odds ratio, water-side, $\\leq$500\\,m, {lab}: all settlements / without rounded coordinates & "
                            f"{o1} ({c1_}) / {o2} ({c2_}) \\\\")
        out["b6_gfm"] = _tbl("tab:b6-gfm", "Diagnostics of the GFM labels and of the near-settlement odds ratios "
                             "(amendment 6). CSI: critical success index for flooded cells. Rounded coordinates: settlements "
                             "whose GeoNames coordinates are whole arc-minutes.", "p{10.6cm}p{4.6cm}", "Quantity & Value", rows)
    pl = (v52 or {}).get("plasmode")
    if isinstance(pl, dict):
        rows = []
        for f, lab in (("A", "all ground nearest the name"), ("B", "within 500\\,m"), ("CA", "all ground, latent"),
                       ("CB", "within 500\\,m, latent")):
            r = pl["forms"][f]
            gains = " & ".join(f"{1e5 * x:.1f}" for x in r["gain"])
            rows.append(f"{lab} & {gains} & {_or_tex(r['delta_M'])} & {_or_tex(r['delta_A'])} & "
                        f"{_or_tex(r['delta_B'])} & {_or_tex(r['mde80'])} \\\\")
        head = ("Target & " + " & ".join(f"{np.exp(x):.1f}" for x in pl["deltas"])
                + " & Margin & Bound & Bound$_0$ & MDE")
        out["b6_plasmode"] = _tbl("tab:b6-plasmode", "Effects injected into the observed labels (plasmode, amendment 6): "
                                  "mean gain in $\\Delta$BS ($\\times10^{-5}$) over the uninjected run for each injected odds "
                                  f"ratio ({pl['reps']} repetitions each), and the odds ratios read from these curves: margin "
                                  "($R=0.05$); upper bound with the mean increment of reassigned names as the no-effect "
                                  "reference (Bound) or zero (Bound$_0$), against the upper 90\\,\\% limit of the observed "
                                  "increment; detectable with 80\\,\\% power (MDE, normal approximation). Latent: effect on a "
                                  "latent indicator of truly water-related names.",
                                  "l" + "r" * len(pl["deltas"]) + "llll", head, rows)
    if v52:
        rows = []
        for key, cmp_, lab in (("label_oracle", "oracle_label", "Flood-label oracle (same event, leaks the outcome)"),
                               ("gsw_oracle", "oracle_gsw", "Recurrent surface water (JRC GSW)"),
                               ("gsw_oracle", "oracle_gsw_nonperm", "Recurrent surface water without permanent water"),
                               ("t0x", "T0x", "Names on T0$^-$ (no river distances, no water land cover)"),
                               ("t0x", "T0", "Names on T0 (same run)"),
                               ("strict", "strict", "Strict reading of the names"),
                               ("nread", "H4", "Language-model reading, primary split"),
                               ("gn", "primary", "Gazetteer without city sections, historical places, duplicates"),
                               ("nohydro", "primary", "Without AOIs lacking a hydrography layer"),
                               ("typed", "primary", "Without AOIs lacking a flood type")):
            b = v52.get(key)
            c = b.get(cmp_) if isinstance(b, dict) else None
            nn = f"{b['n_aoi']} AOIs" if isinstance(b, dict) else "---"
            rows.append(f"{lab} & {_r3(c['R']) if isinstance(c, dict) else '---'} & {_ci90(c)} & "
                        f"{c['decision'] if isinstance(c, dict) else 'pending'} & {nn} \\\\")
        ap = v52.get("actperm")
        if isinstance(ap, dict):
            rows.append(f"Events held out whole: real increment among {ap['n_perms']} reassignments & --- & "
                        f"percentile {ap['percentile_of_real']:.0f} & --- & --- \\\\")
        if v48:
            a = v48["aoi_weighted"]["frozen"]
            rows.append(f"Every AOI weighted equally & {rv(a['R'])} & {rv(a['ci90_R'][0])} to {rv(a['ci90_R'][1])} & "
                        f"{'equivalent' if a['equivalent_TOST'] else 'not equivalent'} & --- \\\\")
        out["b6_runs"] = _tbl("tab:b6-runs", "Post hoc runs and controls of amendment 6, primary split, corrected code. "
                              "$R$ is the within-area increment as a fraction of the terrain increment of the same run, with "
                              "the fixed-margin 90\\,\\% interval and the verdict of the decision table.",
                              "p{8.4cm}rlll", "Run & $R$ & 90\\,\\% interval & Verdict & Sample", rows)
    if v53:
        rows = []
        for gname in ("wetside", "strict", "water", "wet", "riverine", "crossing", "low", "high"):
            cells = []
            for v in ("occ", "seas", "occ_np"):
                c = v53["all_named"][gname][v]
                cells.append(f"{c['estimate']:+.2f} ({c['ci95'][0]:+.2f} to {c['ci95'][1]:+.2f})" if c else "---")
            lab = {"wetside": "water-side (any)", "strict": "strict reading"}.get(gname, g.CLASS_LABELS.get(gname, gname))
            rows.append(f"{lab} & " + " & ".join(cells) + " \\\\")
        out["b6_gsw"] = _tbl("tab:b6-gsw", "Positive control without flood labels (amendment 6): within-AOI difference "
                             "between settlements whose name is in each class and other settlements in the mean occurrence "
                             "of surface water (percentage points of time), its mean seasonality (months) and the occurrence "
                             "without permanent water, within 500\\,m of the settlement point (JRC Global Surface Water "
                             "1984--2021); 95\\,\\% cluster-bootstrap intervals over places.",
                             "lccc", "Names & Occurrence & Seasonality & Occurrence, non-permanent", rows)
    if v48 and v54:
        rows = []
        nc = v48["coverage"]["null_coverage"]
        dc = v48["coverage"]["null_decision_counts"]
        rows.append(f"Null simulations: coverage of the 90\\,\\% / 95\\,\\% interval (mean $R$) & {_pct0(nc['coverage90'])} / "
                    f"{_pct0(nc['coverage95'])} \\\\")
        rows.append(f"Null simulations: coverage of zero, 90\\,\\% / 95\\,\\% & {_pct0(nc['coverage90_zero'])} / "
                    f"{_pct0(nc['coverage95_zero'])} \\\\")
        rows.append(f"Null simulations: between-replicate SD / median bootstrap SE of $R$ & "
                    f"{nc['sd_R_between_reps'] / nc['median_bootstrap_se_R']:.1f} \\\\")
        rows.append("Null simulations: verdicts & " + ", ".join(f"{k} {v}" for k, v in dc.items()) + " \\\\")
        pv = v54["p_values"]["world_v2"]
        rows.append(f"H2: TOST $p$ & $<${pv['H2_tost']['resolution']:.5f} \\\\")
        rows.append(f"H3: one-sided $p$ ($\\Delta_0>\\Delta_4$); Spearman $\\rho$ and exact one-sided $p$ & "
                    f"{pv['H3']['p_delta0_gt_delta4']:.2f}; {pv['H3']['spearman']:.1f}, {pv['H3']['p_spearman_exact']:.2f} \\\\")
        if "H4_tost" in pv:
            rows.append(f"H4: TOST $p$; one-sided $p$ for an improvement & $<${pv['H4_tost']['resolution']:.5f}; "
                        f"{pv['H4_improvement_p']:.2f} \\\\")
        rows.append(f"Holm-adjusted $p$, H3 / H4 & {pv['holm']['H3']:.2f} / {pv['holm']['H4']:.2f} \\\\")
        rc = v54["residual_correlogram"]
        rows.append("Correlation of held-out M1 residuals within AOIs by distance (km) & " + "; ".join(
            f"{b}: {c:.2f}" for b, c in zip(rc["bins_km"], rc["corr"], strict=True)) + " \\\\")
        rp = v54["reference_points"]
        rows.append(f"CSI of maps flagging the observed flooded share per AOI: lowest HAND / M1 / M3n & "
                    f"{rp['hand_threshold_map']['csi']:.2f} / {rp['M1_map']['csi']:.2f} / {rp['M3n_map']['csi']:.2f} \\\\")
        rows.append(f"Brier score: constant (training rate) / M0 / M1 & {rp['brier']['const_train_rate']:.4f} / "
                    f"{rp['brier']['M0']:.4f} / {rp['brier']['M1']:.4f} \\\\")
        out["b6_stats"] = _tbl("tab:b6-stats", "Interval coverage in the pre-specified null simulations (amendment 2), "
                               "the $p$ values required by the protocol (amendment 6, addendum 1), the spatial correlation of "
                               "held-out residuals and reference points for the terrain baseline (primary split). The CSI "
                               "maps use the observed flooded count of each AOI and are therefore optimistic.",
                               "p{9.6cm}p{5.6cm}", "Quantity & Value", rows)
    for name, tex in out.items():
        (HERE / f"app_table_{name}.tex").write_text(tex, encoding="utf-8")


# %% 補足（Supplement）: 付録に残さない表を S 番号で。本文はマクロ \STab<名前>（例 \STabSiYardsticks = S2）で引く
APPENDIX_TABLES = ("comparisons", "splits", "sensitivity", "gfm_compact")
SUPP_TABLES = ("splits_secondary", "yardsticks", "registered", "ratio", "posthoc", "power_splits", "strata", "metrics",
               "lexfix", "coefs", "positive_control", "odds", "class", "band_contrast", "composition", "groups", "b5",
               "gfm", "etym", "b6_lexicon", "b6_forms", "b6_examples", "b6_example", "b6_readings", "b6_flow", "b6_gfm",
               "b6_plasmode", "b6_runs", "b6_gsw", "b6_stats")
_DIGITS = dict(zip("0123456789", ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"),
                   strict=True))


def _label_macro(label: str) -> str:
    parts = re.split(r"[-_:]", label.removeprefix("tab:"))
    name = "".join(p[:1].upper() + p[1:] for p in parts if p)
    return "STab" + "".join(_DIGITS.get(ch, ch) for ch in name)


def _rename_amendments(sx: str) -> str:
    """改訂を番号だけで呼び（付録 B と衝突しないように）、事後の解析に registered を使わない（第三者の登録と区別）。"""
    for a, b in (("Registered quantities", "Quantities planned in advance"), ("Not registered", "Not planned"),
                 ("not registered", "not planned"), ("were registered", "were planned"),
                 ("the registered analyses", "the planned analyses"), ("analysis registered in", "analysis planned in"),
                 ("registered in amendment", "planned in amendment"), ("reporting registered in", "reporting planned in")):
        sx = sx.replace(a, b)
    sx = re.sub(r"B-(\d) addendum~?(\d)", r"amendment \1, addendum~\2", sx)
    sx = re.sub(r"amendment B-(\d)", r"amendment \1", sx)
    sx = re.sub(r"amendments B-(\d) (and|to) B-(\d)", r"amendments \1 \2 \3", sx)
    return re.sub(r"\bB-(\d)\b", r"amendment \1", sx)


def supplement() -> dict[str, str]:
    """補足の表を supp_table_*.tex に書き、付録の表の改訂の名前を直し、S 番号のマクロを返す。"""
    m: dict[str, str] = {}
    k = 0
    lines = []
    for name in SUPP_TABLES:
        f = HERE / f"app_table_{name}.tex"
        if not f.is_file():
            continue
        sx = _rename_amendments(f.read_text(encoding="utf-8"))
        sx = (sx.replace("Appendix~\\ref{app:prereg}", "Appendix~A of the paper")
              .replace("Appendix~\\ref{app:results}", "Appendix~C of the paper")
              .replace("Sect.~\\ref{sec:methods}", "Sect.~2 of the paper")
              .replace("Table~\\ref{tab:etym}", "the table of etymologies in the paper"))
        if "\\resizebox" not in sx and "\\begin{longtable}" not in sx:  # 補足の版面より広い表だけを縮める
            lim = "\\textheight" if "sidewaystable" in sx else "\\textwidth"
            sx = re.sub(r"(\\begin\{tabular\}\{.*?\}\n)", lambda mm, lim=lim: f"\\resizebox{{\\ifdim\\width>{lim}{lim}"
                        f"\\else\\width\\fi}}{{!}}{{%\n{mm.group(1)}", sx, count=1)
            sx = sx.replace("\\end{tabular}", "\\end{tabular}}", 1)
        sx = sx.replace("}\\label{tab:si-posthoc}", " Values of $R$ for subsets with a small terrain increment "
                        "(for example single classes of names) are unstable; the increments by class on the Brier scale "
                        "are in the table of increments by class.}\\label{tab:si-posthoc}")
        sx = sx.replace("}\\label{tab:si-strata}", " Language strata are labelled by the first language of the country "
                        "in GeoNames (for example ur for Pakistan), whereas the lexicon reads the names of a country with "
                        "its blocks (fa and hi for Pakistan).}\\label{tab:si-strata}")
        # 補足は一段組の article: 表をその場に置き（番号の順に並べる）、横置きの表は自分の頁に
        sx = (sx.replace("\\begin{table*}[t]", "\\begin{table}[H]").replace("\\end{table*}", "\\end{table}")
              .replace("\\begin{sidewaystable*}", "\\begin{sidewaystable}")
              .replace("\\end{sidewaystable*}", "\\end{sidewaystable}"))
        if "\\begin{sidewaystable}" in sx:  # 横置きの表の前後で頁を改め、番号の順を保つ
            sx = "\\clearpage\n" + sx + "\\clearpage\n"
        for lab in re.findall(r"\\label\{(tab:[^}]*)\}", sx):
            k += 1
            m[_label_macro(lab)] = f"S{k}"
        (HERE / f"supp_table_{name}.tex").write_text(sx, encoding="utf-8")
        lines.append(f"\\input{{figures/supp_table_{name}}}")
    (HERE / "supp_tables.tex").write_text("% generated by figures/make_figures.py — do not edit\n" + "\n".join(lines)
                                          + "\n", encoding="utf-8")
    supp_labels = {lab for name in SUPP_TABLES if (HERE / f"supp_table_{name}.tex").is_file()
                   for lab in re.findall(r"\\label\{(tab:[^}]*)\}", (HERE / f"supp_table_{name}.tex").read_text())}
    for name in APPENDIX_TABLES:  # 付録に残す表から補足の表への参照は S 番号のマクロに
        f = HERE / f"app_table_{name}.tex"
        if f.is_file():
            sx = _rename_amendments(f.read_text(encoding="utf-8"))
            sx = re.sub(r"Tables?~\\ref\{(tab:[^}]*)\}",
                        lambda mm: (f"Table~\\{_label_macro(mm.group(1))}{{}} of the Supplement"
                                    if mm.group(1) in supp_labels else mm.group(0)), sx)
            f.write_text(sx, encoding="utf-8")
    return m


def table_hazard() -> None:
    """本文の表 3: 確証の結果の翻訳と、検定が見逃しうるもの（数値は全てマクロ）。"""
    rows = [
        r"\multicolumn{3}{l}{\emph{The pre-specified test (within-area increment $R$, 90\,\% interval)}} \\",
        r"Observed $R$ & \RHat{} (\RLoNinety{} to \RHiNinety) & \RWetNear{} (\RWetNearLoNinety{} to \RWetNearHiNinety)$^{a}$ \\",
        r"Largest $R$ if names predicted this ground perfectly & \CeilRAll & \CeilRFiveHundred \\",
        r"\multicolumn{3}{l}{\emph{Odds ratios from effects injected into the observed labels (plasmode)}} \\",
        r"Odds ratio that the margin ($R=0.05$) corresponds to & \PlasMarginA & \PlasMarginB \\",
        r"Upper bound consistent with the observed increment & \PlasBoundA & \PlasBoundB \\",
        r"\quad same, for truly water-related names$^{b}$ & \PlasBoundCA & \PlasBoundCB \\",
        r"\quad same, with zero as the no-effect reference & \PlasBoundZeroA & \PlasBoundZeroB \\",
        r"Detectable with 80\,\% power (normal approximation) & \PlasMdeA & \PlasMdeB \\",
        (r"\multicolumn{3}{l}{\emph{Controls: other information in place of the names, assigned to every cell through "
         r"its nearest settlement ($R$, 90\,\% interval)}} \\"),
        (r"Flooded share around the settlement in the same event$^{d}$ & \multicolumn{2}{l}{\LabOracleR{} "
         r"(\LabOracleRLoNinety{} to \LabOracleRHiNinety)} \\"),
        (r"Recurrent surface water around the settlement, 1984--2021 & \multicolumn{2}{l}{\GswOracleR{} "
         r"(\GswOracleRLoNinety{} to \GswOracleRHiNinety)} \\"),
        r"HAND of the settlement's site & \multicolumn{2}{l}{\OracleR{} (\OracleRLoNinety{} to \OracleRHiNinety)} \\",
        r"\multicolumn{3}{l}{\emph{Direct estimates: odds ratio of flooding, water-side against other names (95\,\% interval)}} \\",
        r"EMS labels, all AOIs, unconditional & \ORWetAllRaw{} (\ORWetAllRawLo--\ORWetAllRawHi) & \ORWetNearRaw{} (\ORWetNearRawLo--\ORWetNearRawHi) \\",
        r"EMS labels, all AOIs, given terrain & \ORWetAllAdj{} (\ORWetAllAdjLo--\ORWetAllAdjHi) & \ORWetNearAdj{} (\ORWetNearAdjLo--\ORWetNearAdjHi) \\",
        r"GFM labels, given terrain$^{c}$ & \GFMORWetAll{} (\GFMORWetAllLo--\GFMORWetAllHi) & \GFMORWetNear{} (\GFMORWetNearLo--\GFMORWetNearHi) \\",
        r"GFM labels with seasonal water assessable (V1), given terrain & \VOneORAll{} (\VOneORAllLo--\VOneORAllHi) & \VOneOR{} (\VOneORLo--\VOneORHi) \\",
    ]
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[t]\n"
           "\\caption{What the test shows in hazard terms and what it could have missed, for all ground nearest a "
           "settlement with a water-side name (\\MatchWetPct{} of the cells) and for the ground within 500\\,m of such a "
           "settlement. Only the first row is pre-specified; the rest are post hoc (amendments 2, 3, 5 and 6). "
           "$^{a}$~Cells within 500\\,m of a settlement with a water-side name (amendment 2). $^{b}$~Effect injected "
           "on a latent indicator of truly water-related names, drawn with the misreading rates measured against "
           "recorded etymologies. $^{c}$~AOIs and cells that radar can assess; built-up land and seasonal water are "
           "excluded. $^{d}$~Leaks the outcome by construction (the cell itself is left out): an upper bound, not a "
           "possible reading of names.}\\label{tab:hazard}\n\\small\n\\begin{tabular}{lll}\n\\tophline\n"
           " & All ground nearest the name & Within 500\\,m of the settlement \\\\\n\\middlehline\n"
           + "\n".join(rows) + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    (HERE / "tab_hazard.tex").write_text(tex, encoding="utf-8")


def main() -> None:
    m, primary, extra, post = base_macros()
    fig_design_nhess(m.get("NLexLangs", "47"))
    v41, v42, v43 = _json("v41_gfm_independent_labels.json"), _json("v42_blind_review_analyses.json"), \
        _json("v43_wikidata_etymology.json")
    v44, v45, v46, v47 = (_json("v44_wikipedia_etymology.json"), _json("v45_gfm_masking_check.json"),
                          _json("v46_norway_impact.json"), _json("v47_b5_cluster.json"))
    m |= freeze_record() | norway_macros(v46) | gfm_macros(v41, v45, v47) | b5_macros(v42, v47)
    m |= etym_macros(v43, v44)
    v48, v50, v51, v52, v53, v54 = (_json("v48_b6_local.json"), _json("v50_gfm_diagnostics.json"),
                                    _json("v51_readings.json"), _json("v52_b6_cluster.json"),
                                    _json("v53_gsw_positive_control.json"), _json("v54_b6_descriptives.json"))
    m |= b6_ceiling_macros(v48) | b6_gfm_macros(v50) | b6_readings_macros(v51) | b6_cluster_macros(v52)
    m |= b6_gsw_macros(v53) | b6_desc_macros(v54) | llmcut_macros() | fs_frozen_macros(v47)
    m |= ti_spread_macros(v42, v47) | actgrp_macros()
    if primary:
        fig_main(primary, extra, post["v36"], v41, v42, v47, v48, v52)
        if primary.get("ladder"):
            fig_ladder_b6(primary["ladder"], primary["primary"]["delta_bs_M0_M1"], v52)
    if post["v37"]:
        fig_names(post["v37"], post["v39"], v41, v42, v50, v52)
    t = pd.read_parquet(DATA / "verification" / "v41_table.parquet",
                        columns=["aoi", "y", "ix", "iy", "p_M1c", "p_M3n"])
    gs = g._json(DATA / "global" / "gfm_summary.json")
    ex = pick_example(t, gs) if gs else None
    m |= fig_example(ex, t) if ex else {"ExampleAoiName": missing("example AOI")}
    m |= timeline(m) | gfm_items()
    write_new_tables(v41, v42, v44, v45, v47)
    table_hazard()
    (HERE / "app_table_gfm_compact.tex").write_text(table_gfm_compact(v41, v47, v50), encoding="utf-8")
    b6_supp_tables(v48, v50, v51, v52, v53, v54)
    tables = copernicus_tables()
    m |= supplement()
    write_numbers(m)
    print(f"macros: {len(m)}; tables: {len(tables)}; missing: {len(g.MISSING)}")
    for w in sorted(set(g.MISSING)):
        print("  missing:", w)


if __name__ == "__main__":
    main()
