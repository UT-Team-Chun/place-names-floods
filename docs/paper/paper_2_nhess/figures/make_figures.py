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
SESOI_DELTA = 1.2  # 凍結したプロトコル §5: 主分割で R ≈ 0.05 になるように較正した合成の効果（power_v3_sesoi）


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
        return {k: missing("v46") for k in ("NorwayCells", "NorwayMatched", "NorwayShare", "NorwayCellsPct",
                                            "NorwayNAoi")}
    return {"NorwayCells": _n(v["n_cells_norway"]), "NorwayMatched": _n(v["n_cells_norway_matched"]),
            "NorwayShare": f"{100 * v['share_matched_of_all']:.2f}\\,\\%",
            # 辞書が読まなかったのはノルウェーの集落に最も近いセル全て（初見の査読 第 5 回: 0.34% は直した辞書が
            # 照合するセルの割合で、影響の範囲ではない）
            "NorwayCellsPct": f"{100 * v['n_cells_norway'] / v['n_cells']:.1f}\\,\\%",
            "NorwayNAoi": str(v["n_aoi_norway"]),
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
    ax.text(17.5, 28.5, "terrain increment\n(the yardstick)", fontsize=7.4, ha="center", va="top",
            color=grey)
    ax.text(74.8, 28.5, "names increment\nR = names / terrain increment", fontsize=7.4, ha="center", va="top",
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
        # 分割の名前を一行ずつ点で区切る（初見の査読 第 5 回: 続けて並んでいて区切りが読めない）
    box(87.5, 1.5, 12, 10, "• location groups\n  (primary)\n• countries\n• floods from 2024\n• floods after the\n  LM cut-off",
        size=6.4)
    fig.savefig(HERE / "fig_design.pdf")
    plt.close(fig)


XMAX_FIG4 = 0.08


def _rows90(primary: dict, extra: dict, v36: dict | None) -> list[tuple[str, float | None, list | None, list | None, str]]:
    """図 4a の行（事前に固定した比較、凍結したコード）: (表示名, R̂, 90% 区間, 95% 区間, 様式)。主分割の行を先に。"""
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
             "main"),
            ("   Position only, flexible models (H1, descriptive)",
             fm("H1_flexible", primary.get("H1_flexible"), "ratio_to_terrain", "ci95_ratio"), "secondary"),
            ("Total incl. regional composition (secondary)", blk(primary.get("total")), "secondary"),
            ("Flexible model refitted with names (secondary)", blk(primary.get("flexible")), "secondary")]
    if (met := extra.get("country")) and met.get("primary"):
        spec.append(("Countries held out", blk(met["primary"]), "main"))
    if (met := extra.get("temporal")) and met.get("primary"):
        und = met["primary"]["decision"] == "undefined-denominator"
        spec.append(("Time held out (floods from 2024)", blk(met["primary"]), "undefined" if und else "main"))
    if (met := extra.get("llmcut")) and met.get("H4"):
        # 星印: 判定はプロトコルの物差しでだけ equivalent（位置のモデルが定数より悪く物差しが膨らむ、表 C1 の注）
        spec.append(("Language model, floods after its cut-off (H4) *", blk(met["H4"]), "main"))
        spec.append(("   Lexicon and character model, same floods", blk(met.get("primary")), "main"))
        spec.append(("   Same, flexible models refitted (secondary)", blk(met.get("flexible")), "secondary"))
        if met.get("H4_total"):
            spec.append(("   Language model incl. its composition (secondary)", blk(met["H4_total"]), "secondary"))
    return [(n, v[0], v[1], v[2], st) for n, v, st in spec if v is not None]


def _draw_rows(ax, rows: list, style: dict) -> None:
    ys = np.arange(len(rows))[::-1]
    ax.axvspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axvline(0, color="0.4", lw=0.8)
    for y, (_, r, c9, c5, st) in zip(ys, rows, strict=True):
        col, ls, filled = style[st]
        if r is None:
            ax.text(0.0, y, "R undefined: terrain did not improve on position", ha="center", va="center",
                    fontsize=6, color=col, style="italic")
            continue
        if c5 is not None:  # 95% は細く、90%（判定に使う）は太く
            ax.plot(np.clip(c5, -XMAX_FIG4, XMAX_FIG4), [y, y], color=col, lw=0.7, ls=ls)
        if c9 is not None:
            ax.plot(np.clip(c9, -XMAX_FIG4, XMAX_FIG4), [y, y], color=col, lw=2.2, ls=ls, solid_capstyle="butt")
        ax.plot([np.clip(r, -XMAX_FIG4, XMAX_FIG4)], [y], "o", ms=4.0, color=col, mfc=col if filled else "white",
                zorder=3)
        hi = max(r, *(c5[1:] if c5 else []), *(c9[1:] if c9 else []))
        if hi > XMAX_FIG4:  # 軸の外に出る区間は端に矢じりと、同じ行の軸の外に上限の値を（重ならないように）
            txt = f"R {r:.2f}, to {hi:.2f}" if r > XMAX_FIG4 else f"to {hi:.2f}"
            ax.plot([XMAX_FIG4], [y], marker=">", ms=3.5, color=col, clip_on=False, zorder=4)
            ax.text(XMAX_FIG4 + 0.004, y, txt, fontsize=5.5, color=col, ha="left", va="center", clip_on=False)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.4)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(-XMAX_FIG4, XMAX_FIG4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def fig_main(primary: dict, extra: dict, v36: dict | None, v41: dict | None, v42: dict | None,
             v47: dict | None, v48: dict | None, v52: dict | None, checks: list, v56: dict | None,
             v59: dict | None = None) -> None:
    """図 4: (a) 事前に固定した比較（凍結したコード）、(b) 事後の確認、(c) 区間ではなく範囲。"""
    v47 = v47 or {}
    rows_a = _rows90(primary, extra, v36)
    lf = (v47.get("lexfix") or {}).get("primary")
    if isinstance(lf, dict):  # 確証の比較の直したコードでの再実行（amendment 4）
        rows_a.insert(1, ("   Same, corrected code (Norwegian names read)", lf["R"], lf["ci90_R"], lf["ci95_R"],
                          "main"))
    rows_b = []
    c2 = (v59 or {}).get("C2")
    rf = (c2.get("corrected") or {}).get("refit") if isinstance(c2, dict) else None
    if isinstance(rf, dict):  # B-8 の C2: 二段目を bootstrap の中で当て直した区間（直したコード）
        rows_b.append(("Second stage refitted in each bootstrap replicate", rf["R"], rf["ci90_R"], rf["ci95_R"], "post"))
    for lab, _am, code, c in checks:
        if not isinstance(c, dict) or c.get("R") is None:
            continue
        c95 = c.get("ci95_R")
        if c95 is None and c.get("ci95_delta") and c.get("yardstick_M0_M1"):
            c95 = [v / c["yardstick_M0_M1"] for v in c["ci95_delta"]]
        rows_b.append((lab + (" (frozen code)" if code == "frozen" else ""), c["R"], c.get("ci90_R"), c95, "post"))
    fs = (v47.get("fold_seeds") or {}).get("fixed_code") or {}
    pm5 = v47.get("permutations")
    pm7 = (v56 or {}).get("permutations") if isinstance((v56 or {}).get("permutations"), dict) else None
    rows_c = []
    if fs.get("R"):
        rows_c.append(("fs", f"{fs.get('n_available', len(fs['R']))} assignments of location groups to folds"))
    if isinstance(pm5, dict):
        rows_c.append(("pm5", f"Names reassigned {pm5['n_perms']} times, real composition"))
    if pm7:
        rows_c.append(("pm7", f"Names reassigned {pm7['n_perms']} times, composition recomputed"))
    style = {"main": ("#1f6fb4", "-", True), "secondary": ("0.45", "-", False),
             "undefined": ("0.45", "--", False), "post": ("#8c5a2b", "-", False)}
    na, nb, nc = len(rows_a), len(rows_b), max(len(rows_c), 1)
    fig, axs = plt.subplots(3, 1, figsize=(6.6, 0.122 * (na + nb + nc) + 1.3), sharex=True,
                            gridspec_kw={"height_ratios": [na + 0.6, nb + 0.6, nc + 0.6]})
    _draw_rows(axs[0], rows_a, style)
    for y, (nm, *_) in zip(np.arange(na)[::-1], rows_a, strict=True):
        if nm.startswith("Countries held out"):
            axs[0].axhline(y + 0.5, color="0.8", lw=0.6)
    _draw_rows(axs[1], rows_b, style)
    ax = axs[2]
    ax.axvspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axvline(0, color="0.4", lw=0.8)
    ys = np.arange(len(rows_c))[::-1]
    for y, (key, _lab) in zip(ys, rows_c, strict=True):
        if key == "fs":
            rs = np.array(list(fs["R"].values()), dtype=float)
            ax.plot(rs, np.full(len(rs), y), "o", ms=3.2, color="#8c5a2b", mfc="white", alpha=0.9)
        else:
            d = pm5 if key == "pm5" else pm7["recomputed_composition"]
            ax.plot([d["R_perm_q025"], d["R_perm_q975"]], [y, y], color="0.7", lw=5, solid_capstyle="butt")
            ax.plot([d["R_real"]], [y], "o", ms=4.2, color="#8c5a2b", mfc="#8c5a2b")
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _, lab in rows_c], fontsize=6.4)
    ax.set_ylim(-0.7, len(rows_c) - 0.3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for a, (lab, title) in zip(axs, (("a", "Pre-specified (frozen code)"), ("b", "Exploratory checks"),
                                     ("c", "Ranges, not intervals")), strict=True):
        a.set_title(f"({lab}) {title}", fontsize=7.5, loc="left", fontweight="bold")
    axs[2].set_xlabel("R (improvement from names / improvement from terrain)", fontsize=8)
    axs[2].tick_params(axis="x", labelsize=7)
    g.SOURCE["Fig4"] = [{"panel": p, "comparison": n.strip(), "R": r, "ci90_low": c9[0] if c9 else None,
                         "ci90_high": c9[1] if c9 else None, "ci95_low": c5[0] if c5 else None,
                         "ci95_high": c5[1] if c5 else None, "group": st}
                        for p, rows in (("a", rows_a), ("b", rows_b)) for n, r, c9, c5, st in rows]
    fig.tight_layout(h_pad=0.6)
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
              v50: dict | None = None, v52: dict | None = None, v56: dict | None = None) -> None:
    fig = plt.figure(figsize=(7.2, 2.9), layout="constrained")
    gs = fig.add_gridspec(1, 2)
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    classes = list(g.CLASS_LABELS)
    ys = np.arange(len(classes))[::-1]
    blue, brown, grey = "#1f6fb4", "#8c5a2b", "0.45"
    pc = v37["positive_control"]
    # a: 集落の点での平均 HAND の差（塗り）と全セル（白抜き）。連続の量で（初見の査読 第 4 回: 2 m 以下の割合は
    # どの群でも 7〜9 割で区別しにくい）
    ax = axs[0]
    for y, c in zip(ys, classes, strict=True):
        col = brown if c == "high" else blue
        for scale, dy, filled in (("settlements", 0.15, True), ("cells", -0.15, False)):
            r = pc[scale][c]["hand_mean"]
            if not r:
                continue
            ax.plot(np.array(r["ci95"]), [y + dy, y + dy], color=col if filled else grey, lw=1.4)
            ax.plot([r["estimate"]], [y + dy], "o", ms=4, color=col if filled else grey,
                    mfc=(col if filled else "white"))
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([g.CLASS_LABELS[c] for c in classes], fontsize=7)
    ax.set_xlabel("difference in mean height above the nearest drainage\n(m; filled: at the settlement, open: all "
                  "cells nearest the name)", fontsize=7)
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
    for a, lab in zip(axs, "ab", strict=True):
        a.tick_params(labelsize=6.5)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.set_title(f"({lab})", fontsize=9, fontweight="bold", loc="left")
    fig.savefig(HERE / "fig_names.pdf")
    plt.close(fig)


# %% 図（第 5 節）: 集落から 500 m 以内の浸水のオッズ比（地形の予測を所与）。区分・版・川の大きさ・同じセルの二つのラベル
def fig_near(v55: dict | None, v57: dict | None, v61: dict | None = None) -> None:
    v57, v61 = v57 or {}, v61 or {}
    b5, b6, b7, b12 = (v57.get(k) or {} for k in ("B5", "B6", "B7", "B12"))
    blue, grey, green, black = "#1f6fb4", "0.45", "#2b8c4a", "0.15"
    groups: list[tuple[str, list]] = []
    ca = b6.get("class_adjusted") or {}
    rows = [(lab, ca[c]["odds_ratio"], ca[c]["bonferroni_ci95"], blue, True) for c, lab in CLASS_LABELS if c in ca]
    fp = ((v55 or {}).get("odds_ems") or {}).get("floodplain")
    if fp:
        rows.append(("Floodplain reading (defined after these)", fp["odds_ratio"], fp["ci95"], blue, False))
    groups.append(("By element class (Bonferroni intervals)", rows))
    rows = []
    for key, lab in (("main", "All AOIs"), ("typed_aois", "AOIs with a recorded flood type"),
                     ("without_no_hydro", "Without AOIs lacking hydrography"),
                     ("without_oceania_americas", "Without Oceania and the Americas"),
                     ("cleaned_reading", "Cleaned reading"), ("uncapped", "Every evaluable cell (no cap)")):
        d = b6.get(key)
        if isinstance(d, dict) and d.get("pairs_ci95"):
            rows.append((lab, d["odds_ratio"], d["pairs_ci95"], black, key == "main"))
    d = v61.get("without_pakistan_bangladesh")
    if isinstance(d, dict) and d.get("pairs_ci95"):
        rows.append(("Without Pakistan and Bangladesh", d["odds_ratio"], d["pairs_ci95"], black, False))
        v64 = _json("v64_b9_aggregate.json") or {}  # B-9 の C5・C7（本文 5.2 節と同じ版、pairs bootstrap）
    for key, lab in (("C5_near", "One definition of permanent water"), ("C7_near", "Bare-earth terrain model")):
        d = v64.get(key)
        if isinstance(d, dict) and d.get("pairs_ci95"):
            rows.append((lab, d["odds_ratio"], d["pairs_ci95"], black, False))
    mc = b6.get("misclassification_corrected")
    if isinstance(mc, dict):
        rows.append(("Truly water-related names (corrected)", mc["median"], mc["ci95"], black, False))
    groups.append(("Pooled water-side class and its variants (pairs bootstrap)", rows))
    rs = b7.get("odds_by_river_size") or {}
    rows = [(lab, rs[k]["odds_ratio"], rs[k]["ci95"], grey, True)
            for k, lab in (("lt100", "Nearest river < 100 km$^2$"), ("100_1000", "Nearest river 100–1000 km$^2$"),
                           ("ge1000", "Nearest river ≥ 1000 km$^2$")) if k in rs]
    if b12.get("near_500m"):
        rows.append(("Name contains a nearby river or lake name", b12["near_500m"]["odds_ratio"],
                     b12["near_500m"]["ci95"], grey, False))
    groups.append(("By river size; reading by river names", rows))
    sc = b5.get("same_cells") or {}
    rows = [(lab, sc[k]["odds_ratio"], sc[k]["pairs_ci95"], col, True)
            for k, lab, col in (("ems", "EMS labels", blue), ("gfm", "GFM labels", green)) if k in sc]
    groups.append(("The two label sets on the same cells (pairs bootstrap)", rows))
    n = sum(len(r) + 1 for _, r in groups)
    fig, ax = plt.subplots(figsize=(6.6, 0.165 * n + 0.75), layout="constrained")
    y, ticks, labels, src = n - 1, [], [], []
    for title, rows in groups:
        ax.text(0.205, y, title, fontsize=7, fontweight="bold", va="center", ha="left")
        y -= 1
        for lab, orr, ci, col, filled in rows:
            lo, hi = max(ci[0], 0.2), min(ci[1], 5.0)
            ax.plot([lo, hi], [y, y], color=col, lw=1.3)
            ax.plot([orr], [y], "o", ms=3.8, color=col, mfc=col if filled else "white", zorder=3)
            if ci[1] > 5.0:
                ax.annotate(f"{ci[1]:.1f}", xy=(5.0, y), xytext=(4.1, y + 0.45), fontsize=5.5, color=col,
                            arrowprops={"arrowstyle": "->", "color": col, "lw": 0.6})
            ticks.append(y)
            labels.append(lab)
            src.append({"group": title, "estimate": lab, "odds_ratio": orr, "ci95": list(ci)})
            y -= 1
    ax.axvline(1, color="0.4", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 5.0)
    ax.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.25", "0.5", "1", "2", "4"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel("odds ratio of flooding within 500 m, given the terrain model's prediction", fontsize=7.5)
    ax.tick_params(axis="x", labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    g.SOURCE["FigNear"] = src
    fig.savefig(HERE / "fig_near.pdf")
    plt.close(fig)


def fig_plasmode(v52: dict | None, v56: dict | None, v59: dict | None = None, v64: dict | None = None) -> None:
    """図 6（第 4.3 節）: 観測のラベルに注入した効果のオッズ比ごとの利得と、余白・観測の上限。

    B-9 の結果があれば、その格子（形 A・CA は低い δ を足す）と、浸水を減らす向き（OR < 1、C2）・集落ごとに違う効果
    （C3）を描く。模擬した効果には全ての曲線に印を付ける（初見の査読 第 5 回、査読者 2 の細かい指摘 11）。
    無効果の基準は B-7 の作り直した置換の平均。本文の幅（12 cm）で作る。
    """
    blue, brown = "#1f6fb4", "#8c5a2b"
    fig, ax = plt.subplots(figsize=(4.7, 3.9), layout="constrained")
    pl = (v52 or {}).get("plasmode")
    b1 = (v59 or {}).get("B1") if isinstance((v59 or {}).get("B1"), dict) else None
    v = v64 or {}
    b9 = v.get("B1_B2") if isinstance(v.get("B1_B2"), dict) and v["B1_B2"].get("gain_mean") else None
    if not isinstance(pl, dict):
        ax.text(0.5, 0.5, "plasmode results pending", ha="center", va="center", fontsize=7, transform=ax.transAxes)
        fig.savefig(HERE / "fig_plasmode.pdf")
        plt.close(fig)
        return
    spec = (("A", "all ground", blue, "-", "o"), ("B", "within 500 m", brown, "-", "o"),
            ("CA", "all ground, truly water-related", blue, ":", "s"),
            ("CB", "500 m, truly water-related", brown, ":", "s"))
    src: dict = {}
    for f, lab, col, ls, mk in spec:
        if b9:
            dl = np.asarray(b9["grid"][f], dtype=float)
            gm, gs, nrep = b9["gain_mean"][f], b9["gain_mc_se"][f], b9["reps"].get(f, "?")
        elif b1:
            dl = np.asarray(pl["deltas"], dtype=float)
            gm, gs, nrep = b1["gain_mean"][f], b1["gain_mc_se"][f], b1["reps"].get(f, "?")
        else:
            dl = np.asarray(pl["deltas"], dtype=float)
            gm, gs, nrep = pl["forms"][f]["gain"], pl["forms"][f]["gain_se"], pl.get("reps", 10)
        ors = np.exp(np.concatenate([[0.0], dl]))
        yv = 1e4 * np.concatenate([[0.0], np.asarray(gm, dtype=float)])
        se = 1e4 * np.concatenate([[0.0], np.asarray(gs, dtype=float)])
        ax.plot(ors, yv, ls, color=col, lw=1.1, marker=mk, ms=2.8, mfc=col if ls == "-" else "white", mew=0.7,
                label=f"{lab} ({nrep} replicates)")
        ax.fill_between(ors, yv - se, yv + se, color=col, alpha=0.15, lw=0)
        src[f] = {"odds_ratio": ors.tolist(), "gain": (yv / 1e4).tolist()}
    c2 = v.get("C2") if isinstance(v.get("C2"), dict) else None
    if c2:  # 浸水を減らす向き（OR < 1）: 同じ色、下向きの三角
        for f, col in (("A", blue), ("B", brown)):
            ors = np.exp(-np.concatenate([[0.0], np.asarray(c2["abs_deltas"], dtype=float)]))
            yv = 1e4 * np.concatenate([[0.0], np.asarray(c2["gain_mean"][f], dtype=float)])
            se = 1e4 * np.concatenate([[0.0], np.asarray(c2["gain_mc_se"][f], dtype=float)])
            ax.plot(ors, yv, "-", color=col, lw=1.1, marker="v", ms=3.0,
                    label=("effects that remove flooding, all ground and 500 m "
                           f"({c2['reps'].get(f, '?')} replicates)") if f == "A" else None)
            ax.fill_between(ors, yv - se, yv + se, color=col, alpha=0.15, lw=0)
            src[f"{f}_negative"] = {"odds_ratio": ors.tolist(), "gain": (yv / 1e4).tolist()}
    c3 = v.get("C3") if isinstance(v.get("C3"), dict) else None
    if c3:  # 集落ごとに違う効果（平均の OR で）
        ors = np.exp(np.concatenate([[0.0], np.asarray(c3["deltas"], dtype=float)]))
        yv = 1e4 * np.concatenate([[0.0], np.asarray(c3["gain_mean"], dtype=float)])
        ax.plot(ors, yv, "-.", color=blue, lw=1.1, marker="D", ms=2.4, mfc="white", mew=0.7,
                label=f"all ground, effect differing between settlements ({c3['reps']} replicates)")
        src["AH"] = {"odds_ratio": ors.tolist(), "gain": (yv / 1e4).tolist()}
    pr = (v56 or {}).get("plasmode_refit")
    if isinstance(pr, dict):  # B-7 の 2: 全てのモデルを当て直した利得（形 A、3 回の平均）
        ors_r = np.exp(np.concatenate([[0.0], np.asarray(pr["deltas"], dtype=float)]))
        gr = np.concatenate([[0.0], np.asarray(pr["gain_refit"], dtype=float)])
        ax.plot(ors_r, 1e4 * gr, "x", ms=4.5, color=blue, mew=1.2, ls="none",
                label=f"all ground, models refitted ({pr.get('reps', 3)} replicates)")
    if b9:
        d_null, u90 = b9["references"]["b7_recomputed"], b9["intervals"]["fixed"]["U90"]
    else:
        d_null, u90 = (b1["references"]["b7_recomputed"] if b1 else pl["delta_bs_null_perm_mean"]), pl["U90"]
    room = u90 - d_null
    ax.axhline(1e4 * pl["margin"], color="0.2", lw=0.9, ls="--", label="margin (R = 0.05)")
    ax.axhline(1e4 * room, color="0.2", lw=0.9, label="observed upper limit, above the no-effect reference")
    ax.axvline(1.0, color="0.6", lw=0.6)
    ax.set_ylim(-0.3, 2.0 * 1e4 * pl["margin"])
    ax.set_xscale("log")
    ticks = ([0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20] if c2 else [1, 2, 5, 10, 20])
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("injected odds ratio of flooding", fontsize=7.5)
    ax.set_ylabel("gain in the names increment (ΔBS, 10$^{-4}$)", fontsize=7.5)
    fig.legend(fontsize=6.5, frameon=False, loc="outside lower center", ncol=1, handlelength=2.2)
    g.SOURCE["Fig6"] = {"null_reference": d_null, "U90": u90, "U90_minus_null": room, "margin": pl["margin"],
                        "curves": src}
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.savefig(HERE / "fig_plasmode.pdf")
    plt.close(fig)


# %% 図 3（本文）: 例の AOI（選び方は結果に依らない: GFM で評価でき、場所の上限で切られていない AOI のうち、
# 地図が読めるよう集落が 30 以下のもののうち、辞書が水系側と読む集落の最も多いもの。同数なら解析したセルの
# 多いもの。浸水ラベル・GFM で評価できる割合・名前の効果は選び方に使わない）
MAX_EXAMPLE_SETTLEMENTS = 30
HAND_CLIP_M = 20.0  # 例の図の HAND の色の上限（m）
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


def pick_example_near(t: pd.DataFrame, gfm_summary: dict) -> str | None:
    """図 3 の 2 つ目の AOI（第 3 回の査読を受けて、候補を見る前に固定した規則）: GFM でも解析でき、場所の上限で間引かれて
    いない AOI のうち、水系側の名の集落から 500 m 以内の解析セルが最も多いもの（同数なら解析セルの多いもの）。"""
    inc = {a for a, s in gfm_summary["aois"].items() if s.get("include")}
    n_tab = t.groupby("aoi").size()
    near = t[(t["P_d1"] <= 500) & (t["lex_wetside_n1"] > 0)].groupby("aoi").size()
    best, best_key = None, (-1, -1)
    for a in sorted(inc):
        evaluable = int(n_tab.get(a, 0))
        if evaluable == 0 or evaluable < 0.99 * gfm_summary["aois"][a].get("n_eval", 0):
            continue
        key = (int(near.get(a, 0)), evaluable)
        if key > best_key:
            best, best_key = a, key
    return best


def _example_data(aoi: str, t: pd.DataFrame, lex, gn: pd.DataFrame) -> dict:
    d = DATA / "global" / "labels" / aoi
    meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
    cells = pd.read_parquet(d / "cells.parquet", columns=["ix", "iy", "f_flood_any"])
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
    cover = grid(cells["f_flood_any"].to_numpy(dtype=float), cells["ix"].to_numpy(), cells["iy"].to_numpy())
    cover[np.isnan(ems)] = np.nan  # 25 m の地図で浸水した割合（解析したセルだけ）
    prob = grid(s["p_M1c"].to_numpy(dtype=float), s["ix"].to_numpy(), s["iy"].to_numpy())
    diff = grid((s["p_M3n"] - s["p_M1c"]).to_numpy(), s["ix"].to_numpy(), s["iy"].to_numpy())
    st = settlements_in(meta, gn[gn["fclass"] == "P"])  # 解析と同じ集落（P 類）だけ
    st = st.assign(wet=[any(lex.match(r.asciiname, r.cc)[c] for c in ("water", "wet", "riverine", "crossing", "low"))
                        for r in st.itertuples()])
    return {"aoi": aoi, "d": d, "meta": meta, "s": s, "gfm": gfm, "ems": ems, "gf": gf, "hand": hand, "diff": diff,
            "cover": cover, "prob": prob, "st": st, "aspect": ny / nx}


def _scale_bar(ax, meta: dict) -> None:
    """格子のセル（meta の cell_m）で描いた図に、幅のおよそ 4 分の 1 の切りのよい長さの縮尺を左下に。"""
    cell_km = float(meta["grid"]["cell_m"]) / 1000
    width_km = meta["grid"]["nx"] * cell_km
    km = next((v for v in (20, 10, 5, 2, 1) if v <= width_km / 3), 1)
    n = km / cell_km
    x0, y0 = 0.06 * meta["grid"]["nx"], 0.94 * meta["grid"]["ny"]
    ax.plot([x0, x0 + n], [y0, y0], color="k", lw=2.2, solid_capstyle="butt")
    ax.text(x0 + n / 2, y0 - 0.025 * meta["grid"]["ny"], f"{km} km", ha="center", va="bottom", fontsize=6)


def _image_dates(aoi: str) -> str:
    """EMS の事後画像の日付（GFM を探した範囲 ± WINDOW_DAYS から戻す、gfm_meta.json）。"""
    import datetime as dt

    from topsem.world.gfm import WINDOW_DAYS
    gm = g._json(DATA / "global" / "labels" / aoi / "gfm_meta.json") or {}
    w = gm.get("window")
    if not w:
        return missing(f"image date {aoi}")
    a = dt.date.fromisoformat(w[0]) + dt.timedelta(days=WINDOW_DAYS)
    b = dt.date.fromisoformat(w[1]) - dt.timedelta(days=WINDOW_DAYS)
    return a.isoformat() if a == b else f"{a.isoformat()} to {b.isoformat()}"


def fig_example(aois: list[str], t: pd.DataFrame) -> dict[str, str]:
    gn = pd.read_parquet(DATA / "global" / "geonames" / "aoi_subset.parquet",
                         columns=["gid", "asciiname", "cc", "lat", "lon", "fclass"]).drop_duplicates("gid")
    lex = g_lexicon()
    data = [_example_data(a, t, lex, gn) for a in aois]
    from matplotlib.colors import ListedColormap
    lab_cmap = ListedColormap(["0.82", "#1f6fb4"])
    heights = [min(3.0, 1.32 * x["aspect"] + 0.95) for x in data]
    # 5 枚を横に並べ、AOI ごとに 1 行。HAND は 20 m で切る（初見の査読 第 4 回: 0〜10 m が見えるように）
    fig, axs = plt.subplots(len(data), 5, figsize=(7.2, sum(heights)), layout="constrained", squeeze=False,
                            gridspec_kw={"height_ratios": heights})
    letters = "abcdefghij"
    hmaxes = []

    def bar(im, ax, label: str) -> None:
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.9, pad=0.02)
        cb.set_label(label, fontsize=6.5)
        cb.ax.tick_params(labelsize=6)
        cb.ax.xaxis.set_major_locator(plt.MaxNLocator(3))

    for i, x in enumerate(data):
        row = axs[i]
        im = row[0].imshow(x["cover"], cmap="Blues", vmin=0, vmax=1, interpolation="nearest")
        lab = np.where(np.isnan(x["ems"]), 0.0, x["ems"])
        row[0].contour(lab, levels=[0.5], colors="k", linewidths=0.35)
        bar(im, row[0], "share flooded at 25 m (outline: label)")
        _scale_bar(row[0], x["meta"])  # 縮尺（査読者 3 の細かい指摘 14、第 5 回）
        row[1].imshow(x["gf"], cmap=lab_cmap, vmin=0, vmax=1, interpolation="nearest")
        hmax = HAND_CLIP_M
        hmaxes.append(hmax)
        im = row[2].imshow(np.clip(x["hand"], 0, hmax), cmap="YlOrBr_r", interpolation="nearest")
        bar(im, row[2], f"HAND (m; clipped at {hmax:.0f} m)")
        im = row[3].imshow(x["prob"], cmap="Blues", vmin=0, vmax=1, interpolation="nearest")
        bar(im, row[3], "held-out probability (M1c)")
        lim = float(np.nanquantile(np.abs(x["diff"]), 0.99)) or 1e-3
        im = row[4].imshow(x["diff"], cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
        bar(im, row[4], "M3n − M1c (probability)")
        st = x["st"]
        ms = 2.6 if len(st) <= 60 else 1.8
        for j, (ax, title) in enumerate(zip(row, ("EMS map", "GFM labels", "HAND", "terrain model", "change from names"),
                                            strict=True)):
            ax.plot(st.loc[~st["wet"], "col"], st.loc[~st["wet"], "row"], "k^", ms=ms, mfc="white", mew=0.5)
            ax.plot(st.loc[st["wet"], "col"], st.loc[st["wet"], "row"], "^", ms=ms + 0.6, color="#e66101", mew=0.5)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"({letters[5 * i + j]}) {title}", fontsize=7.5, loc="left")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(color="#1f6fb4", label="flooded (GFM)"), Patch(color="0.82", label="not flooded (GFM)"),
               Line2D([], [], marker="^", ls="", color="#e66101", label="settlement, water-side name"),
               Line2D([], [], marker="^", ls="", color="k", mfc="white", label="other settlement")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    fig.savefig(HERE / "fig_example.pdf", bbox_inches="tight", pad_inches=0.04)  # 枠の外の集落の印で題が欠けないように
    plt.close(fig)
    out: dict[str, str] = {}
    if len(data) > 1:  # 2 つ目の AOI の名前と HAND の上限
        x = data[1]
        act = x["meta"]["activation"]
        country = ", ".join(act.get("countries") or [])
        name = act.get("name", "").replace("- ", ", ")
        where = name if country and country in name else f"{name}, {country}"
        out |= {"ExampleTwoAoiName": f"{x['aoi'].replace('_', chr(92) + '_')}: {where}, "
                                     f"{act.get('activation_time', '')[:4]}",
                "ExampleTwoHandClip": f"{hmaxes[1]:.0f}", "ExampleTwoNSettlements": str(len(x["st"])),
                "ExampleTwoImageDate": _image_dates(x["aoi"]),
                "ExampleTwoNWetSettlements": str(int(x["st"]["wet"].sum()))}
    x = data[0]
    aoi, d, meta, s, gfm, st = x["aoi"], x["d"], x["meta"], x["s"], x["gfm"], x["st"]
    hmax = hmaxes[0]
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
    # 同じ要素の綴りの揺れ（ASCII の ae と a、査読者 3 の細かい指摘 14）: この AOI にあるときだけ例を書く
    nm = set(st["asciiname"].astype(str))
    pair = [n for n in ("Kittilae", "Ala-Kittila", "Yli-Kittila") if n in nm]
    translit = (" The gazetteer's ASCII forms can spell one element differently within an AOI ("
                + ", ".join(f"\\emph{{{n}}}" for n in pair) + "), which the matching does not reconcile."
                if len(pair) >= 2 else "")  # noqa: PLR2004
    (HERE / "app_table_b6_example.tex").write_text(_tbl(
        "tab:b6-example", f"Settlements of the example AOI of Fig.~3 ({tag}) and their readings by the lexicon "
        "(populated places of the gazetteer within the AOI's grid; water-side: water, wetland, riverine, crossing or "
        "lowland element)." + translit, "llcp{9cm}", "Name & Country & Water-side & Matched forms (position, class)",
        rows_ex), encoding="utf-8")
    out |= {"ExampleAoiName": f"{tag}: {where}, {act.get('activation_time', '')[:4]}",
            "ExampleNSettlements": str(len(st)), "ExampleNWetSettlements": str(int(st['wet'].sum())),
            "ExampleHandClip": f"{hmax:.0f}", "ExampleImageDate": _image_dates(aoi)}
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
             .replace("\\bottomrule", "\\bottomhline").replace("\\begin{table}[h!]", "\\begin{table*}[tp]")
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
        if f.name == "si_table_splits.tex":  # H4 の TOST は物差しの膨らみに依る（初見の査読 第 4 回、査読者 2）
            s = s.replace("Language-model cut-off: H4, M1rc$\\to$M3rn", "Language-model cut-off: H4, M1rc$\\to$M3rn$^{\\dagger}$")
            s = s.replace("}\\label{tab:si-splits}", (
                " $^{\\dagger}$Equivalent only with the protocol's yardstick, which a position model that predicts worse "
                "than the best constant inflates on this split; with three of the four alternative yardsticks the 90\\,\\% "
                "interval exceeds the margins (Table~\\ref{tab:si-yardsticks}).}\\label{tab:si-splits}"), 1)
        if f.name == "si_table_comparisons.tex":  # 本文（第 3 節）に置く表: 版面の高さに収める
            s = s.replace("\\footnotesize\\setlength{\\tabcolsep}{4pt}", "\\scriptsize\\setlength{\\tabcolsep}{4pt}")
            s = _add_h3_rows(s)
        out = HERE / f.name.replace("si_table_", "app_table_")
        out.write_text(s, encoding="utf-8")
        names.append(out.name)
    return names


def _add_h3_rows(s: str) -> str:
    """表 2 に H3 の事前に決めた量（T0 と T4 の差と片側の p、T0 で取り戻した地形の情報の割合）を足す（編集者の指摘 8）。"""
    met = g._json(DATA / "runs" / "world_v2" / "metrics_primary.json") or {}
    lad = met.get("ladder")
    key = "H3 ladder, within areas, T4 = M1 (M1c$\\to$M3n)"
    if not lad or key not in s:
        return s
    d01 = met["primary"]["delta_bs_M0_M1"]
    j = s.index("\\\\", s.index(key)) + 2
    rows = (f"\n\\quad T0 minus T4, one-sided $p={lad['p_one_sided_delta0_gt_delta4']:.2f}$ & "
            f"{g._r(lad['delta0_minus_delta4'] / d01)} & {g._ci([x / d01 for x in lad['ci95_delta0_minus_delta4']])} & "
            "--- & --- & --- & --- \\\\"
            f"\n\\quad share of terrain information recovered at T0$^{{\\ddagger}}$ & "
            f"{g._r(lad['terrain_info_recovered_at_T0'])} & --- & --- & {g._ci(lad['ci95_recovered_at_T0'])} & --- & --- "
            "\\\\")
    s = s[:j] + rows + s[j:]
    return s.replace("}\\label{tab:si-comparisons}", " $^{\\ddagger}$The names increment at T0 as a share of what "
                     "elevation and HAND add (T0 to T4), with a ratio-bootstrap interval.}\\label{tab:si-comparisons}", 1)


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
         "BFive": "B-5_blind_review_analyses.md", "BSix": "B-6_blind_review2_analyses.md",
         "BSeven": "B-7_blind_review3_analyses.md", "BEight": "B-8_blind_review4_analyses.md",
         "BNine": "B-9_blind_review5_analyses.md"}
# 改訂のファイルを変えた後のコミット（追記・記録）→ 年表の英語の記述（開発リポジトリのコミット）
AMEND_UPDATES = {
    "3debd72": "Amendment 1: result of the corrected run recorded",
    "c8859c0": "Addendum 1 to amendment 2", "031d7bd": "Addendum 2 to amendment 2",
    "b7b5ae0": "Addendum 3 to amendment 2", "7f8a551": "Addendum 4 to amendment 2",
    "108627f": "Addendum 1 to amendment 3 (unreadable GFM images left out and recorded)",
    "4172d39": "Addendum 1 to amendment 3, first supplement (images without assets)",
    "487a085": "Addendum 1 to amendment 3, second supplement (record of missing layers)",
    "ec8f82e": "Addendum 2 to amendment 3 (check of GFM masking, after the GFM result had been seen)",
    "ca70160": "Addendum 1 to amendment 5 (element groups)",
    "0311334": "Addendum 2 to amendment 5 (etymology procedure and coding rules)",
    "811ca3c": "Addendum 3 to amendment 5 (fold assignments rerun with the corrected code)",
    "0c81e13": "Addendum 4 to amendment 5 (etymology sample enlarged before coding)",
    "b9f7eac": "Addendum 1 to amendment 6", "04e10f5": "Addendum 2 to amendment 6 (implementation record)",
    "743eb19": "Addendum 3 to amendment 6 (correction of the implementation record)",
    "2631b63": "Addendum 1 to amendment 7",
    "b5ecdcd": "Addendum 2 to amendment 7 (implementation record)",
    "5f32acd": "Addendum 3 to amendment 7 (joint test of the name-feature coefficients)",
    "bdc129b": "Addendum 1 to amendment 8 (implementation record)",
    "53b1fb5": "Addendum 2 to amendment 8 (implementation record, before any result of items B5 to B7 and B13)",
    "8f9d542": "Addendum 3 to amendment 8 (cluster storage failure; computations moved to a local machine)",
    "a80bc9c": "Addendum 4 to amendment 8 (changes to the plan, before any result of items C1 to C5)",
    "cc84ba4": "Addendum 5 to amendment 8 (storage restored; items C1 and C3 back on the cluster, before their results)",
    "51c4f37": "Addendum 6 to amendment 8 (the Pakistan stratum and estimates without South Asia, before computing)",
    "201359f": "Addendum 7 to amendment 8 (precision of the strict reading by language, before computing)",
    "1886f6b": "Addendum 8 to amendment 8 (one job of item C1 rerun on the cluster, before its results)",
        "b2ec94a": "Addendum 9 to amendment 8 (refitted-model odds ratios withdrawn, after they had been seen)",
    "c0bb8c1": "Addendum 1 to amendment 9 (implementation record, before any result of items C1 to C4 and C7)",
}


def _amend_updates() -> list[tuple[str, str]]:
    """改訂のファイルを変えたコミットのうち AMEND_UPDATES にあるもの: (記述, UTC)。"""
    import datetime as dt
    out = []
    for f in AMEND.values():
        try:
            log = subprocess.run(["git", "log", "--format=%h %cI", "--", f"docs/protocol/amendments/{f}"], cwd=REPO,
                                 capture_output=True, text=True, check=True).stdout.split("\n")
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
        for line in log:
            if not line.strip():
                continue
            h, t = line.split()
            if h[:7] in AMEND_UPDATES:
                out.append((AMEND_UPDATES[h[:7]],
                            dt.datetime.fromisoformat(t).astimezone(dt.UTC).strftime("%Y-%m-%dT%H:%MZ")))
    return out
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
    m["PowerSesoiOR"] = f"{np.exp(SESOI_DELTA):.1f}"  # 編集者の指摘（第 5 回）: H4 は情報を与えない
    fj = g._json(DATA / "runs" / "world_v2" / "predictions" / "primary" / "fold0.json") or {}
    ex = ((fj.get("stacks") or {}).get("M1c") or {}).get("extra")  # 地域構成 = 名前の特徴の AOI 平均（同じ数）
    m["NNameFeat"] = str(len(ex)) if ex else missing("name features")
    pw = g._json(DATA / "runs" / "power_v3_null" / "power_summary.json") or {}
    conf = g._json(DATA / "runs" / "world_v2" / "summary_primary.json") or {}
    gs = g._json(DATA / "global" / "gfm_summary.json") or {}
    blind1 = _first_commit_time("docs/research/reviews/2026-10-26_paper2_blind1_editor.md")
    blind2 = _first_commit_time("docs/research/reviews/2026-10-27_paper2_blind2_editor.md")
    blind3 = _first_commit_time("docs/research/reviews/2026-10-27_paper2_blind3_editor.md")
    blind4 = _first_commit_time("docs/research/reviews/2026-10-27_paper2_blind4_editor.md")
    blind5 = _first_commit_time("docs/research/reviews/2026-10-29_paper2_blind5_editor.md")
    interp = _first_commit_time("docs/research/2026-10-26_v2_interpretation_plan.md")
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
            ("Interpretation plan written (after the freeze, before any confirmatory result)", interp or ""),
            ("Third LLM-agent critiques of the manuscript recorded", blind3 or ""),
            ("Amendment 7 (post hoc analyses after the third critiques)", out["AmendTimeBSeven"]),
            ("Fourth LLM-agent critiques of the manuscript recorded", blind4 or ""),
            ("Amendment 8 (exploratory analyses after the fourth critiques)", out["AmendTimeBEight"]),
            ("Fifth LLM-agent critiques of the manuscript recorded", blind5 or ""),
            ("Amendment 9 (exploratory analyses after the fifth critiques)", out["AmendTimeBNine"])]
    rows += _amend_updates()
    del addendum, addenda23  # 追記は _amend_updates から（全ての改訂について同じ規則で）
    rows.sort(key=lambda r: r[1])
    body = "\n".join(f"{a} & {b} \\\\" for a, b in rows)
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table}[tp]\n"
           "\\caption{Timeline of the study (UTC), in order. Times are those of the version-controlled development "
           "repository, of the run records and of the container image (whose tag is its build time). Amendment $k$ is "
           "the file B-$k$ of the protocol's amendments in the code repository; an addendum is a later commit to that "
           "file.}\\label{tab:timeline}\n\\scriptsize\n"
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
    return ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[tp]\n"
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


def _t0_note() -> str:
    """v65: T0 の段の名前の増分の区間が、凍結したコード（表 2）より直したコードで広い理由（編集者の指摘 7、記述）。"""
    v = _json("v65_t0_rung.json")
    if not v:
        return " " + missing("v65")
    fr, co = v["frozen"]["ci95_R_T0"], v["corrected"]["ci95_R_T0"]
    top = v["top_units"][:2]
    names = " and ".join(f"{u['aoi'].replace('_', ' ')} ({u['country']})" for u in top)
    ratio = v["corrected"]["spread"] / v["frozen"]["spread"]
    rest = v["spread_without_top2"]["corrected"] / v["spread_without_top2"]["frozen"]
    return (f" At the T0 rung, the 95\\,\\% interval of $R$ is {rv(co[0])} to {rv(co[1])} with the corrected code "
            f"against {rv(fr[0])} to {rv(fr[1])} with the frozen code (Table~2 of the paper), so a 90\\,\\% interval "
            "here can be wider than the 95\\,\\% interval there. The spread of the location groups' contributions "
            f"grew by a factor of {ratio:.2f}, almost all of it from two large AOIs, {names}, whose contributions changed "
            f"when the models were refitted; without them the factor is {rest:.2f}.")


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
    ems_here = _orc(v41["odds_ems"]["wetside_adjusted_M1"]["terms"].get("wetside:le500m"))[0]
    ems_c2 = (_orc(v41["odds_ems_c2_M1"]["wetside_adjusted_M1"]["terms"].get("wetside:le500m"))[0]
              if v41.get("odds_ems_c2_M1") else missing("v41 EMS c2"))
    cap = ("Independent labels from GFM (post hoc, amendment B-3) on the \\GFMNAoi{} AOIs assessable on both label "
           "sets. Odds ratios come from logistic regressions with an intercept for each AOI and the logit of the "
           "held-out terrain prediction of the run with the frozen code, with standard errors clustered by place; "
           "pct.\\ points are within-AOI contrasts. The EMS odds ratio within 500\\,m of water-side names on these "
           f"AOIs differs in the second decimal with the terrain prediction it is conditioned on: {ems_here} here, "
           "\\SameEmsOR{} with the prediction of the run with the corrected reading (the paper and the table of the "
           f"pooled odds ratio near settlements) and {ems_c2} with the terrain model refitted on GFM labels. "
           "$^{a}$Added after the result above had been seen (B-3 addendum~2). $^{b}$90\\,\\% interval.")
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
                    f"{g._ci(o['ci90_R'])}$^{{a}}$ & {o['decision']} \\\\")  # 表 3 と同じ 90%（査読者 3、第 5 回）
    for key, what in (("notrace", "flood traces counted as not flooded"), ("nonwater", "flooded share over land")):
        lv = (v47.get("label_variants") or {}).get(key)
        if isinstance(lv, dict):
            rows.append(f"$R$, labels with {what} & {g._r(lv['R'])} & {g._ci(lv['ci90_R'])}$^{{a}}$ & "
                        f"{lv['decision']} \\\\")
    cap = ("Post hoc analyses of amendment B-5 on the primary split. The largest odds ratios are upper bounds by test "
           "inversion under the simulation model (the largest effect whose simulated increments place the observed "
           "increment above their 5th percentile). $^{a}$90\\,\\% interval. $^{b}$2.5th to 97.5th percentile of the "
           "reassignments. $^{c}$90\\,\\% interval combining the bootstrap error of the first assignment with the "
           "spread over assignments.")
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
           "500\\,m of a dike or polder name was flooded, so that odds ratio cannot be estimated. The expectation for "
           "island names concerns the site of the settlement, which on a river island stands on its highest ground, not "
           "the ground that the name refers to, which lies within the floodplain. $^{a}$Includes German "
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
         "GFMCSISOne": f"{ag['sensor_s1']['csi']:.2f}", "GFMCSISar": f"{ag['sensor_sar']['csi']:.2f}",
         "GFMHitPct": pct1(ag["all"]["p_gfm1_given_ems1"])}  # GFM が拾った EMS の浸水セル（査読者 1、第 5 回）
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
    act_keys = ("ActPermPct", "ActPermPctOrd", "ActPermN", "ActPermPlacPct", "ActPermRReal", "ActPermRMean", "ActPermRLo",
                "ActPermRHi")
    if isinstance(ap, dict) and ap.get("yardstick"):
        y = float(ap["yardstick"])
        m |= {"ActPermPct": f"{ap['percentile_of_real']:.0f}", "ActPermPctOrd": _ordinal(ap["percentile_of_real"]),
              "ActPermN": str(ap["n_perms"]),
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
    m |= {"HThreeP": f"{pv['H3']['p_delta0_gt_delta4']:.2f}", "HThreeSpearman": f"\\ensuremath{{{pv['H3']['spearman']:.1f}}}",
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
                "logit of the held-out terrain prediction, standard errors clustered by location group (95\\,\\% "
                "interval); pairs: pairs cluster bootstrap over location groups (1,000 replicates, amendment 6; the "
                "2,000 replicates of amendment 7 give \\PairsWetNearLo--\\PairsWetNearHi{} for the first row, "
                "Sect.~\\ref{sec:near}); LOO: range when leaving out one location group. $^{a}$90\\,\\% interval.",
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
            gains = " & ".join(f"{1e4 * x:.2f}" for x in r["gain"])
            rows.append(f"{lab} & {gains} & {_or_tex(r['delta_M'])} & {_or_tex(r['delta_A'])} & "
                        f"{_or_tex(r['delta_B'])} & {_or_tex(r['mde80'])} \\\\")
        head = ("Target & " + " & ".join(f"{np.exp(x):.1f}" for x in pl["deltas"])
                + " & Margin & Bound & Bound$_0$ & MDE")
        out["b6_plasmode"] = _tbl("tab:b6-plasmode", "Effects injected into the observed labels (plasmode, amendment 6): "
                                  "mean gain in $\\Delta$BS ($\\times10^{-4}$) over the uninjected run for each injected odds "
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
                              "the fixed-margin 90\\,\\% interval and the verdict of the decision table." + _t0_note(),
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
APPENDIX_TABLES = ("comparisons", "splits", "sensitivity")
# 補足の表: 主題ごとの節と、各表が何を示すかの一文（編集者の求め、第 4 回）
SUPP_GROUPS = (
    ("Design, data and timeline", (
        ("b6_flow", "Which labelled AOIs entered the analysis, by inclusion criterion, flood type and continent."),
        ("seven_excluded", "The flooded area of the excluded AOIs, by criterion."),
        ("timeline", "The time of every step of the study, including each amendment and addendum."),
        ("seven_status", "The planning status of each group of reported estimates."))),
    ("Pre-specified analyses", (
        ("splits_secondary", "The secondary comparisons of the protocol under each split."),
        ("registered", "Quantities planned in amendment 2 and not shown elsewhere, including the protocol's two-way "
                       "bootstrap."),
        ("yardsticks", "The within-area comparison judged against five yardsticks and two bootstraps."),
        ("ratio", "Ratio-bootstrap intervals of $R$ for every comparison and split."),
        ("power_splits", "Operating characteristics of the decision rule on synthetic labels for each split."),
        ("strata", "The primary comparison within descriptive strata (H5)."),
        ("metrics", "Secondary metrics of the held-out predictions."),
        ("lexfix", "The confirmatory comparisons with the frozen code and with the corrected reading of Norwegian names."),
        ("coefs", "The coefficients of the second stage."),
        ("b6_stats", "Interval coverage in the null simulations, the $p$ values of the protocol and the spatial "
                     "correlation of residuals."))),
    ("Exploratory checks of the verdict", (
        ("posthoc_models", "The models of the exploratory analyses."),
        ("checks", "The within-area comparison under exploratory variants of sample, labels and reading."),
        ("posthoc", "The primary comparison in subsets of cells and under other fold assignments."),
        ("eight_verdict", "The achieved equivalence bound, the location groups needed to settle the inconclusive "
                          "comparisons, the permutation test, the interval with the second stage refitted and the "
                          "operating characteristics at the margin."),
        ("nine_skill", "The Murphy decomposition, the comparison within AOIs and the skill of flood maps."),
        ("nine_runs", "The primary comparison with one definition of permanent water, other subsampling seeds and a "
                      "bare-earth elevation model."),
        ("composition", "What the regional composition of names adds."),
        ("b6_runs", "Exploratory runs and controls of amendment 6."),
        ("b5", "Test inversion and the other analyses of amendment 5."))),
    ("What the pooled test can detect", (
        ("b6_plasmode", "Gains from effects injected into the observed labels, base models fixed."),
        ("seven_refit", "The same gains with every model refitted on the injected labels."),
        ("eight_bounds", "Bounds, detectable effects and margin translations with different no-effect references, with "
                         "Monte Carlo intervals."),
        ("nine_bounds", "Bounds counted three ways, in both directions and for effects that differ between settlements."),
        ("nine_power", "Simulated power of the decision rule with every model refitted."),
        ("seven_translations", "Translations of the margin and bounds under the synthetic labels of the protocol."),
        ("b6_gsw", "Recurrent surface water around settlements, by name class (no flood labels)."))),
    ("Evidence near settlements", (
        ("positive_control", "Where settlements with names of each class stand relative to the drainage, by distance "
                             "(no flood labels)."),
        ("eight_heights", "Heights of settlements above the main rivers, by name class (no flood labels)."),
        ("odds", "Within-AOI odds ratios of flooding by class and distance."),
        ("seven_pairs", "Odds ratios near settlements by class and distance band, with pairs-bootstrap intervals."),
        ("class", "The names increment by class of the nearest name."),
        ("band_contrast", "Flood contrasts and increments by distance from the settlement."),
        ("groups", "Element groups of the lexicon."),
        ("eight_near", "The pooled odds ratio near settlements under every variant, by river size and by class, and "
                       "the two label sets on the same cells."),
        ("eight_where", "Where the flooded cells near water-side names lie, by country and by activation."),
        ("nine_activations", "The pooled odds ratio near settlements with each activation left out."),
        ("nine_classq", "Adjusted $p$ and $q$ values of the class coefficients."),
        ("gfm", "Every quantity of the analyses on GFM labels (amendment 3)."),
        ("b6_gfm", "Diagnostics of the GFM labels."),
        ("eight_calib", "Calibration and skill of the held-out predictions, by split."))),
    ("How the names were read", (
        ("b6_lexicon", "The lexicon by language block."),
        ("b6_forms", "The water-side forms that matched the most cells."),
        ("eight_audit", "Forms excluded by the lexicon's own rules or prone to transferred names, and the most frequent "
                        "forms by country."),
        ("b6_examples", "Examples of names read as water-side, with their recorded derivations."),
        ("b6_example", "The settlements of an example AOI and their readings."),
        ("b6_readings", "Precision and recall of the readings in more detail."),
        ("etym", "Recorded etymologies by lexicon class and language."),
        ("eight_coder", "Agreement of a second coder and the precision and recall from its codes."),
        ("eight_prompt", "The fixed prompt of the language-model reading."),
        ("seven_wikidata", "Settlements whose namesake is recorded on Wikidata."),
        ("seven_languages", "Precision and recall by language."))),
    ("Further quantities", (
        ("nine_text", "Quantities stated in the paper and not given in another table."),)),
)
SUPP_TABLES = tuple(name for _, items in SUPP_GROUPS for name, _ in items)
PAPER_TABLE_WORDS = {"tab:models": "Table~1 of the paper",
                     "tab:si-comparisons": "the table of pre-specified comparisons in Sect.~3 of the paper",
                     "tab:hazard": "the table of Sect.~4 of the paper",
                     "tab:etym": "the table of etymologies in the paper",
                     "tab:amendments": "Table~A1 of the paper", "tab:si-splits": "Table~C1 of the paper",
                     "tab:si-sensitivity": "Table~C2 of the paper"}
SUPP_NOTES = {name: note for _, items in SUPP_GROUPS for name, note in items}
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


# 旧版の表の「place」（場所グループ）を本文の用語 location group にそろえる（populated places・in place of などは触らない）
_GROUP_TERMS = (
    (r"by 10\\,km places but", r"by the location groups formed with a 10\\,km buffer but"),
    (r"by the 5\\,km places of the primary analysis", r"by the location groups of the primary analysis"),
    (r"\\quad as first run, folds by 5\\,km places", r"\\quad as first run, folds by the 5\\,km groups"),
    (r"Places merged within 10\\,km instead of 5\\,km", r"AOIs buffered by 10\\,km instead of 5\\,km to form location groups"),
    (r"Places merged within", r"Location groups merged within"),
    (r"cells per AOI instead of per place\b", r"cells per AOI instead of per location group"),
    (r"cells per place\b", r"cells per location group"),
    (r"& AOIs & Places &", r"& AOIs & Location groups &"),  # この表の「両群」は名前の群なので略さない
    (r"& Places &", r"& Groups &"),
    (r"applying the same decision table\. ", r"applying the same decision table; groups: location groups. "),
    (r"clustered by place\b", r"clustered by location group"),
    (r"over places\b", r"over location groups"),
    (r"one place at a time", r"one location group at a time"),
    (r"leaving out one place\b", r"leaving out one location group"),
    (r"(\d) places\b", r"\1 location groups"),
    (r"^places held out", r"location groups held out"),
)


def _groups(sx: str) -> str:
    for a, b in _GROUP_TERMS:
        sx = re.sub(a, b, sx, flags=re.MULTILINE)
    return sx


def supplement() -> dict[str, str]:
    """補足の表を supp_table_*.tex に書き、付録の表の改訂の名前を直し、S 番号のマクロを返す。"""
    m: dict[str, str] = {}
    k = 0
    lines, contents, current = [], [], None
    for name in SUPP_TABLES:
        f = HERE / f"app_table_{name}.tex"
        if not f.is_file():
            continue
        sx = _groups(_rename_amendments(f.read_text(encoding="utf-8")))
        sx = (sx.replace("Appendix~\\ref{app:prereg}", "Appendix~A of the paper")
              .replace("Appendix~\\ref{app:results}", "Appendix~C of the paper")
              .replace("Sect.~\\ref{sec:methods}", "Sect.~2 of the paper"))
        for lab, words in PAPER_TABLE_WORDS.items():  # 本文・付録の表への参照は言葉で（補足は別の文書）
            sx = sx.replace(f"Table~\\ref{{{lab}}}", words).replace(f"Tables~\\ref{{{lab}}}", words)
        if "\\resizebox" not in sx and "\\begin{longtable}" not in sx:  # 補足の版面より広い表だけを縮める
            lim = "\\textheight" if "sidewaystable" in sx else "\\textwidth"
            sx = re.sub(r"(\\begin\{tabular\}\{.*?\}\n)", lambda mm, lim=lim: f"\\resizebox{{\\ifdim\\width>{lim}{lim}"
                        f"\\else\\width\\fi}}{{!}}{{%\n{mm.group(1)}", sx, count=1)
            sx = sx.replace("\\end{tabular}", "\\end{tabular}}", 1)
        sx = sx.replace("\\emph{Fold assignment, grouping and seeds (verdict of the decision table)}",
                        "\\emph{Fold assignment, grouping and seeds (frozen code; verdict of the decision table)}")
        sx = sx.replace("alternative yardsticks (post hoc, amendment 2):", "alternative yardsticks (post hoc, amendment 2; "
                        "frozen code, as in the confirmatory analysis):")
        sx = sx.replace("without the regional composition (post hoc, amendment 2).",
                        "without the regional composition (post hoc, amendment 2; frozen code).")
        sx = sx.replace("}\\label{tab:si-posthoc}", " Values of $R$ for subsets with a small terrain increment "
                        "(for example single classes of names) are unstable; the increments by class on the Brier scale "
                        "are in the table of increments by class.}\\label{tab:si-posthoc}")
        sx = sx.replace("(the interval of $R$ is then reversed; Appendix~A of the paper)",
                        "(dividing by a negative increment reverses the interval of $R$)")
        sx = sx.replace("}\\label{tab:b6-forms}", " Overlapping forms are each counted: every name in \\emph{-bach} also "
                        "matches \\emph{-ach}.}\\label{tab:b6-forms}")
        sx = sx.replace("}\\label{tab:si-strata}", " Language strata are labelled by the first language of the country "
                        "in GeoNames (for example ur for Pakistan), whereas the lexicon reads the names of a country with "
                        "its blocks (fa and hi for Pakistan).}\\label{tab:si-strata}")
        # 補足は一段組の article: 表をその場に置き（番号の順に並べる）、横置きの表は自分の頁に
        sx = (re.sub(r"\\begin\{table\*\}(\[[a-zA-Z!]*\])?", lambda _: "\\begin{table}[H]", sx)
              .replace("\\end{table*}", "\\end{table}")
              .replace("\\begin{sidewaystable*}", "\\begin{sidewaystable}")
              .replace("\\end{sidewaystable*}", "\\end{sidewaystable}"))
        if "\\begin{sidewaystable}" in sx:  # 横置きの表は横向きの頁に（PDF の閲覧で回して読まなくてよいように）
            sx = (sx.replace("\\begin{sidewaystable}", "\\begin{landscape}\n\\begin{table}[H]")
                  .replace("\\end{sidewaystable}", "\\end{table}\n\\end{landscape}")
                  .replace("\\ifdim\\width>\\textheight\\textheight", "\\ifdim\\width>\\linewidth\\linewidth"))
        nums = []
        for lab in re.findall(r"\\label\{(tab:[^}]*)\}", sx):
            k += 1
            m[_label_macro(lab)] = f"S{k}"
            nums.append(f"S{k}")
        (HERE / f"supp_table_{name}.tex").write_text(sx, encoding="utf-8")
        group = next(gt for gt, items in SUPP_GROUPS if name in dict(items))
        if group != current:  # 主題の節の見出し（表の番号の順に）
            current = group
            lines.append(f"\\section*{{{group}}}")
            contents.append(f"\\multicolumn{{2}}{{l}}{{\\emph{{{group}}}}} \\\\")
        lines.append(f"\\input{{figures/supp_table_{name}}}")
        contents.append(f"{'--'.join(dict.fromkeys([nums[0], nums[-1]]))} & {SUPP_NOTES[name]} \\\\")
    (HERE / "supp_tables.tex").write_text("% generated by figures/make_figures.py — do not edit\n" + "\n".join(lines)
                                          + "\n", encoding="utf-8")
    (HERE / "supp_contents.tex").write_text(
        "% generated by figures/make_figures.py — do not edit\n\\begin{longtable}{lp{14.5cm}}\n\\toprule\n"
        "Table & What it shows \\\\\n\\midrule\n\\endhead\n" + "\n".join(contents)
        + "\n\\bottomrule\n\\end{longtable}\n", encoding="utf-8")
    supp_labels = {lab for name in SUPP_TABLES if (HERE / f"supp_table_{name}.tex").is_file()
                   for lab in re.findall(r"\\label\{(tab:[^}]*)\}", (HERE / f"supp_table_{name}.tex").read_text())}
    for name in APPENDIX_TABLES:  # 付録に残す表から補足の表への参照は S 番号のマクロに
        f = HERE / f"app_table_{name}.tex"
        if f.is_file():
            sx = _groups(_rename_amendments(f.read_text(encoding="utf-8")))
            sx = re.sub(r"Tables?~\\ref\{(tab:[^}]*)\}",
                        lambda mm: (f"Table~\\{_label_macro(mm.group(1))}{{}} of the Supplement"
                                    if mm.group(1) in supp_labels else mm.group(0)), sx)
            f.write_text(sx, encoding="utf-8")
    return m


def table_hazard() -> None:
    """本文の表 2: 事前に固定した R と、名前の地面が全体の R に寄与しうる量・オッズ比への翻訳・対照・直接の推定。"""
    rows = [
        r"\multicolumn{3}{l}{\emph{The pre-specified test}} \\",
        r"Observed $R$, all cells (90\,\% interval) & \multicolumn{2}{l}{\RHat{} (\RLoNinety{} to \RHiNinety)} \\",
        r"\multicolumn{3}{l}{\emph{What this ground contributes to $R$, on the same pooled scale}} \\",
        (r"Observed contribution (90\,\% interval)$^{a}$ & \ObsRAllWet{} (\ObsRAllWetLo{} to \ObsRAllWetHi) & "
         r"\ObsRFiveHundred{} (\ObsRFiveHundredLo{} to \ObsRFiveHundredHi) \\"),
        r"Largest contribution if names predicted this ground perfectly & \CeilRAll & \CeilRFiveHundred \\",
        r"\multicolumn{3}{l}{\emph{Odds ratios of flooding from effects injected into labels}} \\",
        r"Odds ratio at the margin, observed labels, base models fixed & \PlasMarginA & \PlasMarginB \\",
        r"\quad same, synthetic labels of the power analysis & \ORAtMargin & \ensuremath{>}\,20 \\",
        (r"Gain with every model refitted, as a share of the gain with base models fixed$^{h}$ & "
         r"\PlasRefRatioMin--\PlasRefRatioMax & --- \\"),
        (r"Upper bound from the observed increment, fitted models held fixed$^{f}$ & \BoundNewA{} "
         r"(\BoundNewALo--\BoundNewAHi) & \BoundNewB{} (\BoundNewBLo--\BoundNewBHi) \\"),
        r"\quad with the spread over fold assignments & \BoundFoldA & \BoundFoldB \\",
        r"\quad with the second stage refitted within the bootstrap & \BoundRefA & \BoundRefB \\",
        r"\quad same as the first, for truly water-related names$^{b}$ & \BoundNewCA & \BoundNewCB \\",
        r"Lower bound for effects that remove flooding$^{f}$ & \BoundNegA & \BoundNegB \\",
        (r"Detectable with 80\,\% power, normal approximation (fixed; refitted)$^{f}$ & \MdeNewA; \MdeRefA & "
         r"\MdeNewB; \MdeRefB \\"),
        (r"\multicolumn{3}{l}{\emph{Controls: other information in place of the names, assigned to every cell "
         r"($R$, 90\,\% interval)}} \\"),
        (r"Flooded share around the settlement in the same event$^{c}$ & \multicolumn{2}{l}{\LabOracleR{} "
         r"(\LabOracleRLoNinety{} to \LabOracleRHiNinety)} \\"),
        (r"Recurrent surface water around the settlement, 1984--2021 & \multicolumn{2}{l}{\GswOracleR{} "
         r"(\GswOracleRLoNinety{} to \GswOracleRHiNinety)} \\"),
        r"HAND of the settlement's site & \multicolumn{2}{l}{\OracleR{} (\OracleRLoNinety{} to \OracleRHiNinety)} \\",
        (r"Flooding of the cell itself in other EMS events$^{d}$ & \multicolumn{2}{l}{\RecurCellR{} "
         r"(\RecurCellRLoNinety{} to \RecurCellRHiNinety)} \\"),
        (r"\quad the same, averaged around the settlement and assigned through it & \multicolumn{2}{l}{\RecurSetR{} "
         r"(\RecurSetRLoNinety{} to \RecurSetRHiNinety)} \\"),
        (r"Recurrent surface water of the cell itself, 1984--2021 & \multicolumn{2}{l}{\GswCellR{} "
         r"(\GswCellRLoNinety{} to \GswCellRHiNinety)} \\"),
        (r"\multicolumn{3}{l}{\emph{Direct estimates: odds ratio of flooding, water-side against other names "
         r"(95\,\% interval)}} \\"),
        r"EMS labels, all AOIs, unconditional & \ORWetAllRaw{} (\ORWetAllRawLo--\ORWetAllRawHi) & \ORWetNearRaw{} (\ORWetNearRawLo--\ORWetNearRawHi) \\",
        r"EMS labels, all AOIs, given terrain$^{g}$ & \ORWetAllAdj{} (\ORWetAllAdjLo--\ORWetAllAdjHi) & \ORWetNearAdj{} (\PairsWetNearLo--\PairsWetNearHi) \\",
        r"GFM labels, given terrain$^{e,g}$ & \GFMORWetAll{} (\GFMORWetAllLo--\GFMORWetAllHi) & \SameGfmOR{} (\SameGfmORLo--\SameGfmORHi) \\",
        r"\quad EMS labels on the same cells, given terrain$^{g}$ & --- & \SameEmsOR{} (\SameEmsORLo--\SameEmsORHi) \\",
        r"GFM labels with seasonal water assessable (V1), given terrain & \VOneORAll{} (\VOneORAllLo--\VOneORAllHi) & \VOneOR{} (\VOneORLo--\VOneORHi) \\",
    ]
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[tp]\n"
           "\\caption{What the test shows in hazard terms and what it could have missed, for all ground nearest a "
           "settlement with a water-side name (\\MatchWetPct{} of the cells) and for the ground within 500\\,m of such a "
           "settlement. Only the first row is pre-specified; the rest are exploratory (amendments 2, 3 and 5 to 9, or not "
           "planned, footnote a). Contributions are on the scale of the pooled $R$, that is, the change in the Brier "
           "score of all cells divided by the terrain increment. $^{a}$~The increments on this ground (amendment 2) "
           "re-expressed on the pooled scale (not planned). $^{b}$~Effect injected on a latent indicator of truly "
           "water-related names, drawn with the misreading rates measured against recorded etymologies. $^{c}$~Leaks "
           "the outcome by construction (the cell itself is left out): an upper bound, not a possible reading of names. "
           "$^{d}$~Known for \\RecurCovPct{} of the cells; zero with an indicator elsewhere. $^{e}$~AOIs and cells that "
           "radar can assess; built-up land and seasonal water are excluded. $^{f}$~Base models fixed; no-effect reference: the "
           "mean increment of reassigned names with the regional composition rebuilt (amendment 8); in parentheses, Monte "
           "Carlo 95\\,\\% intervals over the simulated replicates (\\PlasRepsA{} for all ground, \\PlasRepsB{} within "
           "500\\,m). "
           "Bounds with other references are in Tables~\\STabEightBounds{} and \\STabNineBounds{} of the Supplement. "
           "$^{g}$~Within 500\\,m, intervals from a pairs cluster bootstrap over location groups (2,000 replicates); "
           "other odds ratios have cluster-robust intervals. $^{h}$~Three injected odds ratios "
           "(\\PlasRefOrMin{} to \\PlasRefOrMax), three replicates each (amendment 7); this simulation has no effect "
           "small enough to read bounds from (amendment 8, addendum 9).}\\label{tab:hazard}\n\\scriptsize\n"
           "\\begin{tabular}{lll}\n\\tophline\n"
           " & All ground nearest the name & Within 500\\,m of the settlement \\\\\n\\middlehline\n"
           + "\n".join(rows) + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    (HERE / "tab_hazard.tex").write_text(tex, encoding="utf-8")


# %% B-7（初見の査読 第 3 回への事後の解析: v55 手元・v56 クラスタ）
def _or3(prefix: str, d: dict | None) -> dict[str, str]:
    """v55 のオッズ比（odds_ratio・ci95）を {prefix}・Lo・Hi に。"""
    if not isinstance(d, dict) or d.get("odds_ratio") is None:
        return {f"{prefix}{s}": missing(prefix) for s in ("", "Lo", "Hi")}
    return {prefix: orv(d["odds_ratio"]), f"{prefix}Lo": orv(d["ci95"][0]), f"{prefix}Hi": orv(d["ci95"][1])}


def _m1(x: float) -> str:
    """HAND などの差（m、小数 1 桁、数式の負号）。"""
    return f"\\ensuremath{{{x:.1f}}}"


def _corr(x: float | None) -> str:
    if x is None:
        return missing("correlation")
    return f"\\ensuremath{{{x:.3f}}}" if abs(x) < 0.01 else f"\\ensuremath{{{x:.2f}}}"


def _ordinal(p: float) -> str:
    k = int(np.floor(p + 0.5))  # 四捨五入（round は偶数への丸めで 80.5 を 80 にする）
    suf = "th" if 10 <= k % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")
    return f"{k}{suf}"


def _status_word(st: str) -> str:
    return st


def plas_margin_range(v37: dict | None, v52: dict | None, v56: dict | None,
                      v59: dict | None = None) -> tuple[str, list]:
    """余白の翻訳の範囲（B-7 の 2 の規則: B-2 の合成ラベル・基礎を固定した plasmode。B-8 の B1 の 50 回の曲線が
    あればそれを、無ければ B-6 の 10 回を。当て直しの読みは追記 9 により入れない）。"""
    vals, above = [], False
    est = ((v37 or {}).get("effect_size_translation") or {}).get("odds_ratio", {}).get("margin")
    if est is not None:
        vals.append(("synthetic labels (amendment 2)", float(est)))
    b1 = (v59 or {}).get("B1") if isinstance((v59 or {}).get("B1"), dict) else None
    pa = ((v52 or {}).get("plasmode") or {})
    a = ((b1.get("margin_or") or {}).get("A") if b1 else
         (pa.get("forms") or {}).get("A", {}).get("delta_M") if isinstance(pa, dict) else None)
    if isinstance(a, dict) and a.get("odds_ratio") is not None:
        vals.append(("plasmode, base models fixed (amendments 6 and 8)", float(a["odds_ratio"])))
    # B-8 追記 9: 当て直したモデルの換算（B-7 の 2）は粗い格子の直線補間で低く偏るので、範囲に入れない
    if len(vals) < 2:
        return missing("margin range (plasmode)"), vals
    lo, hi = min(v for _, v in vals), max(v for _, v in vals)
    if above:
        return f"{lo:.1f} to more than {np.exp(1.5):.1f}", vals
    return f"{lo:.1f}--{hi:.1f}", vals


def b7_macros(v55: dict | None, v56: dict | None, v37: dict | None, v52: dict | None) -> dict[str, str]:
    m: dict[str, str] = {}
    met = g._json(DATA / "runs" / "world_v2" / "metrics_primary.json") or {}
    br = met.get("brier") or {}
    m |= {"BSMOne": f"{br['M1']:.5f}" if "M1" in br else missing("brier M1"),
          "BSMOnes": f"{br['M1s']:.5f}" if "M1s" in br else missing("brier M1s")}
    cut = g._json(DATA / "runs" / "world_v2" / "metrics_llmcut.json") or {}
    m |= _cmp_macros("FlexLlmRX", cut.get("H4_flexible"))  # 言語モデルの読みの柔軟な比較（M1 → M3r）
    m |= {"RFlexLlmR": m.pop("FlexLlmRXR"), "RFlexLlmRLo": m.pop("FlexLlmRXRLo"), "RFlexLlmRHi": m.pop("FlexLlmRXRHi"),
          "VerdictFlexLlmR": m.pop("FlexLlmRXVerdict")}
    for k in ("FlexLlmRXRLoNinety", "FlexLlmRXRHiNinety", "FlexLlmRXRp"):
        m.pop(k, None)
    # H5: 意味・地形の層の R の範囲（記述、凍結したコード）
    h5 = [v["R_within"] for k, v in (met.get("H5_strata") or {}).items()
          if isinstance(v, dict) and v.get("R_within") is not None and "built" not in k]
    m |= ({"HFiveRMin": rv(min(h5)), "HFiveRMax": rv(max(h5))} if h5 else
          {"HFiveRMin": missing("H5 strata"), "HFiveRMax": missing("H5 strata")})
    # 500 m 以内の観測の増分を全体の尺度に（ΔBS_sub × n_sub / N / 物差し; B-2 の帯の増分の言い換え、計画なし）
    b = ((v37 or {}).get("primary_by_distance_wetside") or {}).get("le500m")
    if isinstance(b, dict):
        n_all = float(met.get("n_cells") or 0)
        d01 = float((met.get("primary") or {}).get("delta_bs_M0_M1") or 0)
        f = b["n_cells"] / n_all / d01 if n_all and d01 else None
        m |= ({"ObsRFiveHundred": rv(b["delta_bs_names"] * f), "ObsRFiveHundredLo": rv(b["ci90_delta"][0] * f),
               "ObsRFiveHundredHi": rv(b["ci90_delta"][1] * f)} if f else {})
    cw = ((v37 or {}).get("primary_by_class") or {}).get("wetside")  # 水系側の名の最寄りの全セル（B-2 の区分別の増分）
    if isinstance(cw, dict) and met.get("n_cells") and (met.get("primary") or {}).get("delta_bs_M0_M1"):
        f2 = cw["n_cells"] / float(met["n_cells"]) / float(met["primary"]["delta_bs_M0_M1"])
        m |= {"ObsRAllWet": rv(cw["delta_bs_names"] * f2), "ObsRAllWetLo": rv(cw["ci90_delta"][0] * f2),
              "ObsRAllWetHi": rv(cw["ci90_delta"][1] * f2)}
    for k, what in (("ObsRAllWet", "v37 class"), ("ObsRAllWetLo", "v37 class"), ("ObsRAllWetHi", "v37 class"),
                    ("ObsRFiveHundred", "v37 band"), ("ObsRFiveHundredLo", "v37 band"),
                    ("ObsRFiveHundredHi", "v37 band")):
        if k not in m:
            m[k] = missing(what)
    h5s = met.get("H5_strata") or {}
    if h5s.get("nearest_name_in_lexicon") and met.get("n_cells"):
        m |= {"MatchAnyPctFrozen": pct(h5s["nearest_name_in_lexicon"]["n"] / met["n_cells"]),
              "MatchWetPctFrozen": pct(h5s["nearest_name_wetside"]["n"] / met["n_cells"])}
    else:
        m |= {"MatchAnyPctFrozen": missing("H5 strata"), "MatchWetPctFrozen": missing("H5 strata")}
    rng, _ = plas_margin_range(v37, v52, v56)
    m["PlasMarginRange"] = rng
    # v55（手元）
    if not v55:
        for k in ("FpOREms", "FpORGfm", "CciOREms", "CciORGfm", "PartialOR", "FracOR", "SevLowOR", "PairsWetNear",
                  "AdjHandWet", "FpHandAdj"):
            m |= {f"{k}{s}": missing("v55") for s in ("", "Lo", "Hi")}
        return m | v56_macros(v56)
    oe, og = v55["odds_ems"], v55["odds_gfm"]
    m |= _or3("FpOREms", oe["floodplain"]) | _or3("FpORGfm", og["floodplain"])
    m |= {"FpEmsFlooded": _n(int(oe["floodplain"]["n_group_flooded"])),
          "FpGfmFlooded": _n(int(og["floodplain"]["n_group_flooded"]))}
    m |= _or3("CciOREms", oe["wetside_adjusted_cci"]) | _or3("CciORGfm", og["wetside_adjusted_cci"])
    st = v55["settlement"]
    a = st["wetside_fe"]["adjusted_cci"]
    m |= {"AdjHandWet": _m1(a["estimate"]), "AdjHandWetLo": _m1(a["ci95"][0]), "AdjHandWetHi": _m1(a["ci95"][1])}
    a = st["floodplain_fe"]["adjusted_cci"]
    m |= {"FpHandAdj": _m1(a["estimate"]), "FpHandAdjLo": _m1(a["ci95"][0]), "FpHandAdjHi": _m1(a["ci95"][1])}
    a = st["floodplain_v37_estimator"]["hand_mean"]
    m |= {"FpHand": _m1(a["estimate"]), "FpHandLo": _m1(a["ci95"][0]), "FpHandHi": _m1(a["ci95"][1])}
    p = v55["partial"]
    m |= {"PartialPct": pct(p["share_partial_all_cells"])}
    m |= _or3("PartialOR", p["odds_without_partial"]) | _or3("FracOR", p["odds_fractional_riverine_cover"])
    sv = v55["severity"]
    m |= _or3("SevLowOR", sv["low"]) | {"SevLowFlooded": _n(int(sv["low"]["n_group_flooded"])),
                                         "SevMidOR": orv(sv["middle"]["odds_ratio"]),
                                         "SevHighOR": orv(sv["high"]["odds_ratio"])}
    gd = v55["gfm_discrimination"]["lexfix"]["M3n_minus_M1c"]
    m |= {"GfmAucDiffMed": f"{gd['auc_median']:.4f}".replace("-0.0000", "0.0000"),
          "GfmAucDiffQOne": f"\\ensuremath{{{gd['auc_iqr'][0]:.4f}}}",
          "GfmAucDiffQThree": f"\\ensuremath{{{gd['auc_iqr'][1]:.4f}}}"}
    pb = (v55.get("pairs_bootstrap") or {}).get("wetside", {}).get("wetside:le500m")
    m |= ({"PairsWetNearLo": orv(pb["pairs_ci95"][0]), "PairsWetNearHi": orv(pb["pairs_ci95"][1])} if pb else
          {"PairsWetNearLo": missing("pairs bootstrap"), "PairsWetNearHi": missing("pairs bootstrap")})
    wt = v55.get("wald_tests", {})
    m |= {"WaldWetP": f"{wt['wetside']['p_F']:.2f}" if wt else missing("Wald"),
          "WaldClassP": f"{wt['classes']['p_F']:.2f}" if wt else missing("Wald")}
    hc = v55.get("heldout_coef_test") or {}
    m |= {"HeldoutWaldP": f"{hc['p_F']:.2f}" if hc else missing("held-out coefficient test"),
          "HeldoutNFeat": str(hc["df"]) if hc else missing("held-out coefficient test")}
    nu = v55["null_outliers"]
    m |= {"NullNOut": str(nu["n_outlying_3rsd"]), "NullSdNoOut": f"{nu['sd_R_without_outlying']:.4f}"}
    ex = v55["excluded_area"]
    m |= {"ExclWetMedKm": f"{ex['wet<20']['median_riverine_km2']:.1f}",
          "ExclWetMaxKm": f"{ex['wet<20']['max_riverine_km2']:.1f}",
          "ExclWetNinetyKm": f"{ex['wet<20']['q10_q90_riverine_km2'][1]:.1f}",
          "InclMedKm": f"{ex['included']['median_riverine_km2']:.0f}"}
    dg = v55.get("digraph_matches") or {}
    dd = v55["digraphs"]
    m |= {"DigraphPct": pct(dd["share_settlements"]),
          "DigraphNSettlements": _n(dg["n_settlements"]) if dg else missing("digraphs"),
          "DigraphGainWet": str(dg["gain_wet"]["n_settlements"]) if dg else missing("digraphs"),
          "DigraphLoseWet": str(dg["lose_wet"]["n_settlements"]) if dg else missing("digraphs")}
    rb = v55["residual_between_aois"]
    m |= {"ResidSameAct": _corr(rb["same_activation"]["all"]["corr"]),
          "ResidSameActN": str(rb["same_activation"]["all"]["n_pairs"]),
          "ResidDiffAct": _corr(rb["different_activation"]["all"]["corr"])}
    eu = v55["europe"]
    m |= _cmp_macros("Eu", eu["europe_lexfix"])
    shares = eu["flooded_share_by_continent"]
    non = [(v["share"], v["n_cells"]) for k, v in shares.items() if k != "EU"]
    m |= {"EuFloodPct": pct(shares["EU"]["share"]),
          "NonEuFloodPct": pct(sum(s * n for s, n in non) / sum(n for _, n in non))}
    h3 = v55["h3_slope"]
    m |= {"HThreeSlope": rv(h3["slope_R_per_rung"]), "HThreeSlopeLo": rv(h3["ci95"][0]),
          "HThreeSlopeHi": rv(h3["ci95"][1])}
    cs = v55["fixed_threshold_csi"]
    m |= {"CsiHalf": f"{cs['p_ge_0.5']['csi_pooled']:.2f}", "CsiTrain": f"{cs['p_ge_train_rate']['csi_pooled']:.2f}"}
    lp = (v55.get("block_precision") or {}).get("by_language") or {}
    lr = ((v55.get("language_recall") or {}).get("by_language")) or {}
    precs = [v["precision_lex"]["rate"] for v in lp.values()]
    recs = [v["recall"] for v in lr.values() if v.get("recall") == v.get("recall")]
    m |= {"LangPrecMin": _pct0(min(precs)) if precs else missing("languages"),
          "LangPrecMax": _pct0(max(precs)) if precs else missing("languages"),
          "LangRecallMin": _pct0(min(recs)) if recs else missing("language recall"),
          "LangRecallMax": _pct0(max(recs)) if recs else missing("language recall"),
          "LangNWord": _word(len(lp)) if lp else missing("languages")}
    return m | v56_macros(v56)


def v56_macros(v56: dict | None) -> dict[str, str]:
    m: dict[str, str] = {}
    v56 = v56 or {}

    def blk(key: str, name: str) -> dict | None:
        b = v56.get(key)
        return b.get(name) if isinstance(b, dict) else None

    rc = v56.get("recur") if isinstance(v56.get("recur"), dict) else {}
    rec = ((rc.get("record") or {}).get("recurrence")) or {}
    m |= {"RecurCovPct": pct(rec["share_cells_cell_cov"]) if rec else missing("recurrence record"),
          "RecurNAoi": str(rec["n_aoi_with_other_events"]) if rec else missing("recurrence record")}
    m |= _cmp_macros("RecurCell", blk("recur", "recur_cell")) | _cmp_macros("RecurSet", blk("recur", "recur_settlement"))
    sub = rc.get("subset_aois_with_other_events") or {}
    m |= _cmp_macros("RecurSubCell", sub.get("recur_cell")) | _cmp_macros("RecurSubSet", sub.get("recur_settlement"))
    m |= _cmp_macros("Fp", blk("floodplain", "floodplain")) | _cmp_macros("Partial", blk("partial", "primary"))
    m |= _cmp_macros("ActFix", blk("actfix", "primary"))
    af = v56.get("actfix") if isinstance(v56.get("actfix"), dict) else None
    fz = af.get("frozen_code_run") if af else None
    if af and isinstance(fz, dict):  # B-7 の 3: 凍結・直したコードの M1c・M3n の Brier（5 桁）
        m |= {"BrierAfFrozenMOnec": f"{fz['brier']['M1c']:.5f}", "BrierAfFrozenMThreen": f"{fz['brier']['M3n']:.5f}",
              "BrierAfFixMOnec": f"{af['brier']['M1c']:.5f}", "BrierAfFixMThreen": f"{af['brier']['M3n']:.5f}",
              "BrierAfMOnes": f"{af['brier']['M1s']:.4f}", "BrierAfMOnecFour": f"{af['brier']['M1c']:.4f}"}
    else:
        m |= {k: missing("actfix") for k in ("BrierAfFrozenMOnec", "BrierAfFrozenMThreen", "BrierAfFixMOnec",
                                             "BrierAfFixMThreen", "BrierAfMOnes", "BrierAfMOnecFour")}
    pm = v56.get("permutations")
    if isinstance(pm, dict):
        rc_, ct = pm["recomputed_composition"], pm["centred"]
        m |= {"PermRcPct": _ordinal(rc_["percentile_of_real"]), "PermRcMean": rv(rc_["R_perm_mean"]),
              "PermRcLo": rv(rc_["R_perm_q025"]), "PermRcHi": rv(rc_["R_perm_q975"]),
              "PermRcPctNum": f"{rc_['percentile_of_real']:.0f}",
              "PermCtPct": _ordinal(ct["percentile_of_real"]), "PermCtMean": rv(ct["R_perm_mean"]),
              "PermCtLo": rv(ct["R_perm_q025"]), "PermCtHi": rv(ct["R_perm_q975"]), "PermCtReal": rv(ct["R_real"])}
        b5 = pm.get("b5_single_placebo")
        m["PlacPctBfive"] = _ordinal(b5["percentile_of_single_placebo"]) if isinstance(b5, dict) else missing("B-5")
        m["PlacRBfive"] = rv(b5["R_placebo"]) if isinstance(b5, dict) else missing("B-5")
    else:
        m |= {k: missing("permutations (B-7)") for k in ("PermRcPct", "PermRcMean", "PermRcLo", "PermRcHi",
                                                         "PermRcPctNum", "PermCtPct", "PermCtMean", "PermCtLo",
                                                         "PermCtHi", "PermCtReal", "PlacPctBfive", "PlacRBfive")}
    pr = v56.get("plasmode_refit")
    if isinstance(pr, dict):
        ds = list(pr["deltas"])
        k1 = ds.index(1.0) if 1.0 in ds else 0
        m |= {"PlasRefRatioMin": f"{min(pr['ratio_refit_to_fixed']):.2f}",
              "PlasRefRatioMax": f"{max(pr['ratio_refit_to_fixed']):.2f}",
              "PlasRefOrMin": f"{min(pr['odds_ratios']):.2f}", "PlasRefOrMax": f"{max(pr['odds_ratios']):.2f}"}
        m |= {"PlasRefRatio": f"{pr['ratio_refit_to_fixed'][k1]:.2f}",
              "PlasRefMargin": pr["delta_M"]["odds_ratio_text"], "PlasRefBound": pr["delta_A"]["odds_ratio_text"],
              "PlasRefReps": str(pr["reps"])}
    else:
        m |= {k: missing("plasmode refit") for k in ("PlasRefRatio", "PlasRefMargin", "PlasRefBound", "PlasRefReps",
                                                     "PlasRefRatioMin", "PlasRefRatioMax", "PlasRefOrMin",
                                                     "PlasRefOrMax")}
    return m


# %% B-8（初見の査読 第 4 回への事後の解析）のマクロ: v57（手元の解析 B2〜B12）・v57b（二人目の符号化者 B13）・
# v59（C1〜C5 と B1 の集計）
AUDIT_SHORT = ("de:au:suffix", "sv:as:suffix", "no:as:suffix", "fa:ab:suffix")  # 辞書の規則が除いた 2 文字の接辞
AUDIT_EN = ("en:burn:suffix", "en:bourne:suffix", "en:borne:suffix", "en:pool:suffix", "en:well:suffix",
            "en:wells:suffix")  # 移された名・人名も読む英語の接尾辞
CLASS_TEX = {"riverine": "Riverine", "wet": "Wet", "water": "Water", "crossing": "Crossing", "low": "Low",
             "high": "High"}


def _f2(x: float) -> str:
    return f"{x:.2f}"


def _thousands(n: float) -> str:
    return f"{int(round(n)):,}".replace(",", "{,}")


def _orp(prefix: str, d: dict | None, key: str = "pairs_ci95") -> dict[str, str]:
    """オッズ比と、key の区間（既定は pairs bootstrap）を {prefix}・Lo・Hi に。"""
    if not isinstance(d, dict) or d.get("odds_ratio") is None or d.get(key) is None:
        return {f"{prefix}{s}": missing(prefix) for s in ("", "Lo", "Hi")}
    return {prefix: orv(d["odds_ratio"]), f"{prefix}Lo": orv(d[key][0]), f"{prefix}Hi": orv(d[key][1])}


def _est3(prefix: str, d: dict | None, fmt=None) -> dict[str, str]:
    """推定値と 95% 区間（estimate・ci95）を {prefix}・Lo・Hi に（既定は m の小数 1 桁）。"""
    f = fmt or _m1
    if not isinstance(d, dict) or d.get("estimate") is None:
        return {f"{prefix}{s}": missing(prefix) for s in ("", "Lo", "Hi")}
    return {prefix: f(d["estimate"]), f"{prefix}Lo": f(d["ci95"][0]), f"{prefix}Hi": f(d["ci95"][1])}


def _near_rule(rule: dict | None) -> str:
    """B-8 の B6 で固定した規則（上限の無い標本の区間が 1 を外れ、型のある AOI・掃除した読みと同じ向きなら変える）。"""
    if not isinstance(rule, dict):
        return missing("B6 rule")
    if rule["conclusion_changes"]:
        return "it did, so the conclusion changes"
    if rule["uncapped_excludes_one"]:
        return "its interval excluded one, but not in the direction of both, so the conclusion stands"
    return "its interval included one, so the conclusion stands"


def lexicon_form_macros() -> dict[str, str]:
    """凍結した辞書の項目の数え方（初見の査読 第 5 回: 2,538 は綴り・位置・区分・ブロックの組）。"""
    sys.path.insert(0, str(REPO / "backend" / "packages" / "topsem"))
    from topsem.world.lexicon import Lexicon

    lx = Lexicon(REPO / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")
    ents = [e for v in lx.by_lang.values() for e in v]
    return {"NLexForms": f"{len({e.form for e in ents}):,}".replace(",", "{,}")}


def _rtext(c: dict | None, what: str) -> str:
    """評価の実行の主比較を本文の句に（$R$、90% 区間、判定）。"""
    p = (c or {}).get("primary") if isinstance(c, dict) else None
    if not isinstance(p, dict):
        return missing(what)
    return (f"$R={rv(p['R'])}$, 90\\,\\% interval {rv(p['ci90_R'][0])} to {rv(p['ci90_R'][1])}, "
            f"\\emph{{{p['decision']}}}")


def c9_macros(v64: dict | None) -> dict[str, str]:
    """B-9 の B1・B2・C2〜C7（v64）。B1 の固定した区間の上限は、低い δ の格子を足した値で B-8 の値を置き換える。"""
    m: dict[str, str] = {}
    v = v64 or {}
    b = v.get("B1_B2") if isinstance(v.get("B1_B2"), dict) else {}

    def orx(e: dict | None, what: str) -> dict[str, str]:
        if not isinstance(e, dict):
            return {"": missing(what), "Lo": missing(what), "Hi": missing(what)}
        return {"": _gt(e["odds_ratio_text"]), "Lo": _gt(e["mc95_text"][0]), "Hi": _gt(e["mc95_text"][1])}

    bd = b.get("bounds") or {}
    for kind, key in (("refit", "BoundRef"), ("fold", "BoundFold")):
        for form in ("A", "B", "CA", "CB"):
            m |= {f"{key}{form}{s}": x for s, x in orx(((bd.get(kind) or {}).get("b7_recomputed") or {}).get(form),
                                                        f"v64 {kind} {form}").items()}
    if b.get("low_ranges"):  # 低い δ が揃えば、固定した区間の上限（主の基準）を合わせた格子の値に
        for form in ("A", "CA"):
            e = ((bd.get("fixed") or {}).get("b7_recomputed") or {}).get(form)
            pre = "BoundNewA" if form == "A" else "BoundNewCA"
            if isinstance(e, dict):
                m |= {f"{pre}{s}": x for s, x in orx(e, "v64 fixed").items()}
    md = b.get("mde80") or {}
    m |= {f"MdeRef{f}": orx((md.get("refit") or {}).get(f), "v64 mde")[""] for f in ("A", "B")}
    my = b.get("margin_by_yardstick") or {}
    m |= {"MarginYardLo": orx(my.get("q025"), "v64 margin")[""], "MarginYardHi": orx(my.get("q975"), "v64 margin")[""]}
    c2 = v.get("C2") if isinstance(v.get("C2"), dict) else None
    for f in ("A", "B"):
        e = ((c2 or {}).get("lower_bounds") or {}).get("b7_recomputed", {}).get(f) if c2 else None
        m |= {f"BoundNeg{f}": (e["odds_ratio_text"] if isinstance(e, dict) else missing("v64 C2"))}
    c3 = v.get("C3") if isinstance(v.get("C3"), dict) else None
    m |= {"BoundHetA": orx((c3 or {}).get("bound", {}).get("b7_recomputed") if c3 else None, "v64 C3")[""],
          "MarginHetA": orx((c3 or {}).get("margin_or") if c3 else None, "v64 C3")[""]}
    c4 = v.get("C4") if isinstance(v.get("C4"), dict) else None
    if c4:
        parts = [f"{x['odds_ratio']:.2f} in {x['n_positive']} of {x['reps']}" for x in c4["by_delta"]]
        m["PowerSimList"] = ", ".join(parts[:-1]) + f" and {parts[-1]}"
        m["PowerSimText"] = (f"the smallest injected odds ratio detected in at least 80\\,\\% of replicates was "
                             f"{c4['simulated_detectable_or']:.2f}" if c4.get("simulated_detectable_or")
                             else c4["rule"])
    else:
        m |= {"PowerSimList": missing("v64 C4"), "PowerSimText": missing("v64 C4")}
    c5 = v.get("C5")
    gp = (c5 or {}).get("record", {}).get("b6", {}).get("gsw_perm") if isinstance(c5, dict) else None
    m["GswPermText"] = (f"with one definition of permanent water for all AOIs ({_n(gp['n_dropped'])} cells removed; "
                        f"{_rtext(c5, 'v64 C5')})" if gp else missing("v64 C5"))
    c6 = v.get("C6") or {}
    rs = [x["primary"]["R"] for x in c6.values() if isinstance(x, dict) and isinstance(x.get("primary"), dict)]
    m["SeedText"] = (f"with {_word(len(rs))} other seeds for the subsampling under the cap ($R$ from {rv(min(rs))} to "
                     f"{rv(max(rs))})" if rs else missing("v64 C6"))
    c7 = v.get("C7")
    runs = [c5, *c6.values(), c7]
    if all(isinstance(x, dict) and isinstance(x.get("primary"), dict) for x in runs):
        decs = {x["primary"]["decision"] for x in runs}
        m["RobustBNineText"] = ("left the verdict unchanged" if decs == {"equivalent-null"}
                                else "gave the verdicts " + ", ".join(sorted(decs)))
    else:
        m["RobustBNineText"] = missing("v64 C5-C7")
    m["FabText"] = (f"with a bare-earth elevation model and heights above the main rivers in place of HAND "
                    f"({_rtext(c7, 'v64 C7')})" if isinstance(c7, dict) else missing("v64 C7"))
    # 集落の近くのオッズ比（B-8 の B6 の主の式、pairs bootstrap）: 永久水域の共通の定義（C5）と FABDEM（C7）
    m |= _orp("GswPermNearOR", v.get("C5_near")) | _orp("FabNearOR", v.get("C7_near"))
    sk = v.get("C7_skill") if isinstance(v.get("C7_skill"), dict) else None
    m["FabCsi"] = f"{sk['M1']['csi']['median']:.2f}" if sk else missing("v64 C7 skill")
    m["FabDeltaCsi"] = (f"{sk['delta_csi_M3n_minus_M1c']['median']:.3f}" if sk else missing("v64 C7 skill"))
    # B-9 の C7 の規則: 物差し（地形の利得）の変化を書く
    m["FabYard"] = sci(c7["yardstick"]) if isinstance(c7, dict) else missing("v64 C7")
    m["FabYardLo"], m["FabYardHi"] = ((sci(c7["ci95_yardstick"][0]), sci(c7["ci95_yardstick"][1]))
                                      if isinstance(c7, dict) else (missing("v64 C7"), missing("v64 C7")))
    return m


def b9_macros(v63: dict | None, v42: dict | None) -> dict[str, str]:
    """B-9 の B3〜B8・A7（v63、手元）: Murphy の分解、AOI の中だけの評価、浸水の地図の技能、activation を一つずつ
    除く推定、区分の BH、詳細図の数、方法の細部の数。"""
    m: dict[str, str] = {}
    v = v63 or {}
    b3 = (v.get("B3") or {}).get("diff_M3n_minus_M1c")
    m |= ({"MurphyRelDiff": sci(b3["reliability"]), "MurphyRelDiffLo": sci(b3["reliability_ci95"][0]),
           "MurphyRelDiffHi": sci(b3["reliability_ci95"][1]), "MurphyResDiff": sci(b3["resolution"]),
           "MurphyResDiffLo": sci(b3["resolution_ci95"][0]), "MurphyResDiffHi": sci(b3["resolution_ci95"][1])}
          if b3 else {k: missing("v63 B3") for k in ("MurphyRelDiff", "MurphyRelDiffLo", "MurphyRelDiffHi",
                                                     "MurphyResDiff", "MurphyResDiffLo", "MurphyResDiffHi")})
    b4 = v.get("B4")
    m |= ({"WithinR": rv(b4["R_within"]), "WithinRLoNinety": rv(b4["ci90_R"][0]),
           "WithinRHiNinety": rv(b4["ci90_R"][1]), "WithinRLo": rv(b4["ci95_R"][0]), "WithinRHi": rv(b4["ci95_R"][1]),
           "WithinYard": sci(b4["yardstick_within"]), "WithinNAoi": str(b4["n_aoi_inner"])}
          if b4 else {k: missing("v63 B4") for k in ("WithinR", "WithinRLoNinety", "WithinRHiNinety", "WithinRLo",
                                                     "WithinRHi", "WithinYard", "WithinNAoi")})
    b5 = v.get("B5")
    if b5:
        for key, tex in (("M1", "MOne"), ("M3n", "MThreen"), ("HAND", "Hand")):
            s = b5[key]
            m |= {f"Map{tex}Csi": _f2(s["csi"]["median"]), f"Map{tex}CsiQOne": _f2(s["csi"]["q25_q75"][0]),
                  f"Map{tex}CsiQThree": _f2(s["csi"]["q25_q75"][1]), f"Map{tex}Hit": _f2(s["hit_rate"]["median"]),
                  f"Map{tex}Far": _f2(s["false_alarm_ratio"]["median"])}
        d = b5["delta_csi_M3n_minus_M1c"]
        m |= {"DeltaCsiMed": f"{d['median']:.3f}", "DeltaCsiLo": f"{d['ci95_median'][0]:.3f}",
              "DeltaCsiHi": f"{d['ci95_median'][1]:.3f}", "DeltaCsiMean": f"\\ensuremath{{{d['mean']:.4f}}}",
              "DeltaCsiPosPct": pct1(d["share_aoi_positive"]), "DeltaCsiNegPct": pct1(d["share_aoi_negative"]),
              "MapNAoi": str(b5["M1"]["n_aoi"])}
    else:
        m |= {k: missing("v63 B5") for k in ("MapMOneCsi", "MapMOneHit", "MapMOneFar", "MapHandCsi", "DeltaCsiMed",
                                             "DeltaCsiLo", "DeltaCsiHi", "DeltaCsiPosPct", "DeltaCsiNegPct",
                                             "MapNAoi", "MapMOneCsiQOne", "MapMOneCsiQThree", "MapMThreenCsi",
                                             "DeltaCsiMean")}
    b6 = v.get("B6")
    m |= ({"LooActMin": orv(b6["min"]["odds_ratio"]), "LooActMinAct": b6["min"]["activation"],
           "LooActMax": orv(b6["max"]["odds_ratio"]), "LooActMaxAct": b6["max"]["activation"],
           "LooActN": str(b6["n_activations"])}
          if b6 else {k: missing("v63 B6") for k in ("LooActMin", "LooActMinAct", "LooActMax", "LooActMaxAct",
                                                     "LooActN")})
    b7 = v.get("B7")
    m["ClassMinBH"] = _f2(min(x["q_bh"] for x in b7.values())) if b7 else missing("v63 B7")
    b8 = v.get("B8")
    m |= ({"DetailAoi": str(b8["with_detail"]), "DetailBoth": str(b8["with_detail_and_overview"])}
          if b8 else {"DetailAoi": missing("v63 B8"), "DetailBoth": missing("v63 B8")})
    a7 = v.get("A7")
    m |= ({"NAoiFiftyM": str(a7["n_aoi_50m_subgrid"]), "NCellsOtherFlood": str(a7["n_cells_excluded_other_flood"]),
           "CoderTwoModel": ", ".join(a7["coder2"]["models"]), "CoderTwoCli": ", ".join(a7["coder2"]["cli_versions"]),
           "CoderTwoDate": a7["coder2"]["date"], "CoderTwoSessions": str(a7["coder2"]["sessions"])}
          if a7 else {k: missing("v63 A7") for k in ("NAoiFiftyM", "NCellsOtherFlood", "CoderTwoModel",
                                                     "CoderTwoCli", "CoderTwoDate", "CoderTwoSessions")})
    # B-5 の 8 で登録した予想に照らす: 中州の群の集落の HAND（全セル、AOI 内の差）
    isl = (((v42 or {}).get("groups") or {}).get("by_group") or {}).get("island") or {}
    h = isl.get("pc_hand_mean_all")
    m |= ({"EGIslandHand": _m1(h["estimate"]), "EGIslandHandLo": _m1(h["ci95"][0]),
           "EGIslandHandHi": _m1(h["ci95"][1])}
          if h else {k: missing("v42 island HAND") for k in ("EGIslandHand", "EGIslandHandLo", "EGIslandHandHi")})
    return m


def b8_macros(v57: dict | None, v57b: dict | None, v43: dict | None, v41: dict | None) -> dict[str, str]:
    m: dict[str, str] = {}
    v57 = v57 or {}
    b2 = (v57.get("B2") or {}).get("world_v2")  # 凍結したコードの実行（判定の実行）
    m |= ({"AchievedBoundR": rv(b2["achieved_bound_R"]), "AchievedBoundBS": sci(b2["achieved_bound_brier"]),
           "ObsHalfWidth": f"{(b2['ci90_R'][1] - b2['ci90_R'][0]) / 2:.4f}"} if b2
          else {k: missing("B2") for k in ("AchievedBoundR", "AchievedBoundBS", "ObsHalfWidth")})
    b3 = v57.get("B3") or {}
    tot, near = b3.get("total_M1s_M3s"), b3.get("near_500m_wetside_not_planned")
    m |= ({"NeedEquivTotalX": f"{tot['ratio_equivalence']:.0f}", "NeedPosTotalX": f"{tot['ratio_positive']:.0f}"}
          if tot else {k: missing("B3") for k in ("NeedEquivTotalX", "NeedPosTotalX")})
    m |= ({"NeedPosNearX": f"{near['ratio_positive']:.1f}", "NearR": rv(near["R"]), "NearRLo": rv(near["ci90_R"][0]),
           "NearRHi": rv(near["ci90_R"][1])} if near
          else {k: missing("B3 near") for k in ("NeedPosNearX", "NearR", "NearRLo", "NearRHi")})
    b4 = v57.get("B4")
    m["PermRcP"] = f"{b4['p_one_sided']:.2f}" if b4 else missing("B4")
    b9 = v57.get("B9") or {}
    m |= ({"SkillMOneAuc": _f2(b9["M1"]["auc_within"]["median"]),
           "SkillHandAuc": _f2(b9["HAND"]["auc_within"]["median"])} if b9.get("M1") and b9.get("HAND")
          else {k: missing("B9") for k in ("SkillMOneAuc", "SkillHandAuc")})
    b8 = v57.get("B8") or {}
    for key, split in (("CalibPrimMOnec", "primary"), ("CalibActMOnec", "activation_corrected"),
                       ("CalibGfmMOnec", "gfm_refit")):
        c = (b8.get(split) or {}).get("M1c")
        m[key] = _f2(c["slope"]) if c else missing(f"B8 {split}")
    # B5: 同じセルでの EMS と GFM（pairs bootstrap 2,000 回・t の jackknife）、氾濫原の読み、Sentinel-1 以外の AOI
    b5 = v57.get("B5") or {}
    sc = b5.get("same_cells") or {}
    m |= _orp("SameEmsOR", sc.get("ems"))
    gf, df_ = sc.get("gfm") or {}, sc.get("log_or_difference_gfm_minus_ems") or {}
    m |= ({"SameGfmOR": orv(gf["odds_ratio"]), "SameGfmORLo": orv(gf["pairs_ci95"][0]),
           "SameGfmORHi": orv(gf["pairs_ci95"][1]),
           "SameGfmJkLo": orv(gf["jackknife_t_ci95"][0]), "SameGfmJkHi": orv(gf["jackknife_t_ci95"][1])} if gf
          else {k: missing("B5 gfm") for k in ("SameGfmOR", "SameGfmORLo", "SameGfmORHi", "SameGfmJkLo",
                                               "SameGfmJkHi")})
    m |= ({"SameRatio": orv(df_["ratio_of_odds_ratios"]), "SameRatioLo": orv(df_["pairs_ci95_ratio"][0]),
           "SameRatioHi": orv(df_["pairs_ci95_ratio"][1])} if df_
          else {k: missing("B5 ratio") for k in ("SameRatio", "SameRatioLo", "SameRatioHi")})
    fe = (b5.get("floodplain_same_cells") or {}).get("ems")
    m |= _orp("FpSameEms", fe, "ci95")
    m["FpSameEmsFlooded"] = _thousands(fe["n_group_flooded"]) if fe else missing("B5 floodplain")
    rd = b5.get("refit_diagnosis") or {}
    m |= ({"GfmRefitAuc": _f2(rd["p_M1"]["auc_within_median"]), "GfmEmsAuc": _f2(rd["p_EMS_M1"]["auc_within_median"])}
          if rd else {k: missing("B5 refit") for k in ("GfmRefitAuc", "GfmEmsAuc")})
    ind = ((v41 or {}).get("sensitivity_independent_sensor") or {}).get("odds_ems") or {}
    t = ((ind.get("wetside_adjusted_M1") or {}).get("terms") or {}).get("wetside:le500m")
    m |= _orp("NonSOneEms", t, "ci95_odds_ratio")
    sens = (v41 or {}).get("sensor_counts_included") or {}
    m["GFMNSOne"] = str(sens["s1"]) if "s1" in sens else missing("v41 sensors")
    # B6: 区分・変種・誤分類の補正（pairs bootstrap 2,000 回、区分は Bonferroni の区間と Holm の p）
    b6 = v57.get("B6") or {}
    ca = b6.get("class_adjusted") or {}
    for c, tex in CLASS_TEX.items():
        m |= _orp(f"ClassOR{tex}", ca.get(c), "bonferroni_ci95")
    m["ClassMinHolm"] = _f2(min(v["p_holm"] for v in ca.values())) if ca else missing("B6 class")
    for key, name in (("VarTypedOR", "typed_aois"), ("VarNoHydroOR", "without_no_hydro"),
                      ("VarNoSettlerOR", "without_oceania_americas"), ("VarCleanOR", "cleaned_reading"),
                      ("VarUncappedOR", "uncapped")):
        m |= _orp(key, b6.get(name))
    un, cl, mc = b6.get("uncapped"), b6.get("cleaned_reading"), b6.get("misclassification_corrected")
    m["VarUncappedFlooded"] = _thousands(un["n_group_flooded"]) if un else missing("B6 uncapped")
    m["CleanChangedPct"] = pct1(cl["share_cells_changed"]) if cl else missing("B6 cleaned")
    m |= ({"MisclassOR": orv(mc["median"]), "MisclassORLo": orv(mc["ci95"][0]), "MisclassORHi": orv(mc["ci95"][1])}
          if mc else {k: missing("B6 misclassification") for k in ("MisclassOR", "MisclassORLo", "MisclassORHi")})
    m["NearRuleOutcome"] = _near_rule(b6.get("rule"))
    # B7: 本川からの高さ（集水 100・1000 km² 以上）と、川の大きさごとのオッズ比
    b7 = v57.get("B7") or {}
    h100, h1000 = b7.get("h_above_100") or {}, b7.get("h_above_1000") or {}
    m |= _est3("MainHundredWet", h100.get("wetside"))
    m |= _est3("MainHundredCrossing", (h100.get("by_class") or {}).get("crossing"))
    m["MainThousandWet"] = _m1(h1000["wetside"]["estimate"]) if h1000.get("wetside") else missing("B7 1000")
    rs = b7.get("odds_by_river_size") or {}
    m |= _or3("RiverLtOR", rs.get("lt100")) | _or3("RiverMidOR", rs.get("100_1000")) | _or3("RiverGeOR", rs.get("ge1000"))
    ha = b7.get("odds_adjusted_h_main")
    m["HMainAdjOR"] = orv(ha["odds_ratio"]) if ha else missing("B7 h_main")
    # B11: 規則が除いた形・移された名を読む形の照合の数
    wf = ((v57.get("B11") or {}).get("watched_forms")) or {}
    if wf:
        n = {k: int(wf[k]["settlements"]) if k in wf else 0 for k in (*AUDIT_SHORT, *AUDIT_EN, "hi:char:word")}
        m |= {"AuditShortN": str(sum(n[k] for k in AUDIT_SHORT)), "AuditAuN": str(n["de:au:suffix"]),
              "AuditSvN": str(n["sv:as:suffix"]), "AuditAbN": str(n["fa:ab:suffix"]),
              "AuditEnN": str(sum(n[k] for k in AUDIT_EN)), "AuditCharN": str(n["hi:char:word"])}
    else:
        m |= {k: missing("B11") for k in ("AuditShortN", "AuditAuN", "AuditSvN", "AuditAbN", "AuditEnN", "AuditCharN")}
    # B12: 近くの川・湖の名を含む集落の名（辞書に依らない読み）
    b12 = v57.get("B12") or {}
    if b12:
        ov = b12["overlap"]
        m |= {"HydSharePct": pct1(b12["share_hydronym"]),
              "HydOverlapPct": pct1(ov["hydronym_and_wetside"] / (ov["hydronym_and_wetside"] + ov["hydronym_only"])),
              "HydRecallK": str(b12["recall_wikidata_water_namesakes"]["k"])}
    else:
        m |= {k: missing("B12") for k in ("HydSharePct", "HydOverlapPct", "HydRecallK")}
    m |= _est3("HydHand", b12.get("positive_control_hand")) | _or3("HydOR", b12.get("near_500m"))
    rl = (v43 or {}).get("recall_lexicon")
    m["WDRecallLexK"] = str(rl["k"]) if rl else missing("v43 recall")
    # B13: 二人目の符号化者（別の開発元のモデル）
    if v57b:
        pr = v57b["coder2_precision_recall"]
        tr = v57b["transferred_among_lexicon_wetside_true"]
        m |= {"CodeTwoAgreeWr": _pct0(v57b["water_related_agreement"]),
              "CodeTwoKappaWr": _f2(v57b["water_related_kappa"]), "CodeTwoKappaStatus": _f2(v57b["status_kappa"]),
              "CodeTwoPrec": _pct0(pr["broad"]["lex_wetside"]["precision"]),
              "CodeTwoRecall": _pct0(pr["broad"]["lex_wetside"]["recall"]),
              "CodeTwoPrecStrict": _pct0(pr["strict"]["strict_flag"]["precision"]),
              "CodeTwoRecallStrict": _pct0(pr["strict"]["strict_flag"]["recall"]),
              "CodeTwoPrecNread": _pct0(pr["broad"]["nread_wetside"]["precision"]),
              "CodeTwoRecallNread": _pct0(pr["broad"]["nread_wetside"]["recall"]),
              "CodeTwoTransferredK": str(tr["yes"]), "CodeTwoTransferredN": str(tr["n"])}
    else:
        m |= {k: missing("B13") for k in (
            "CodeTwoAgreeWr", "CodeTwoKappaWr", "CodeTwoKappaStatus", "CodeTwoPrec", "CodeTwoRecall",
            "CodeTwoPrecStrict", "CodeTwoRecallStrict", "CodeTwoPrecNread", "CodeTwoRecallNread",
            "CodeTwoTransferredK", "CodeTwoTransferredN")}
    return m


def c8_macros(v59: dict | None) -> dict[str, str]:
    """B-8 の C1〜C5 と B1（v59 の集計）。"""
    m: dict[str, str] = {}
    v59 = v59 or {}
    d = lambda k: v59[k] if isinstance(v59.get(k), dict) else {}
    c2 = d("C2").get("corrected") or {}
    for lo, hi, name in (("RefitRLoNinety", "RefitRHiNinety", "refit"),
                         ("RefitFixedLoNinety", "RefitFixedHiNinety", "fixed")):
        c = c2.get(name)
        m |= {lo: rv(c["ci90_R"][0]), hi: rv(c["ci90_R"][1])} if c else {lo: missing("C2"), hi: missing("C2")}
    rf = c2.get("refit")  # 当て直した区間の広がり（R の SD の比）と、当て直しでの達成した境界
    m["RefitSdRatio"] = f"{c2['sd_ratio_refit_to_fixed']:.1f}" if rf else missing("C2")
    m["AchievedBoundRefit"] = rv(max(abs(rf["ci90_R"][0]), abs(rf["ci90_R"][1]))) if rf else missing("C2")
    rule = c2.get("rule")  # B-8 の C2 の規則の文
    c4p = d("C4").get("primary")  # 手元の評価が確証の比較を再現したか（追記 3）
    lf = d("C5").get("lexfix") or {}
    m |= {"CFourPrimRFive": f"\\ensuremath{{{c4p['R']:.5f}}}" if c4p else missing("C4 primary"),
          "LexfixRFive": f"\\ensuremath{{{lf['R']:.5f}}}" if lf.get("R") is not None else missing("lexfix R")}
    cl = v59.get("C1_check_local") if isinstance(v59.get("C1_check_local"), dict) else None
    m["ConeLocalSe"] = f"{cl['max_mean_diff_over_se']:.1f}" if cl else missing("C1 local check")
    cc = v59.get("C1_check") if isinstance(v59.get("C1_check"), dict) else None  # クラスタでの同じ fold の再現（追記 5）
    m["ConeClusterText"] = ((f"where the same fold reproduced the earlier cluster values (mean differences of at most "
                             f"{sci(cc['max_mean_diff_over_se'])} standard errors)" if cc["reproduced"] else
                             "where the same fold did not reproduce the earlier cluster values either") if cc
                            else missing("C1 cluster check"))
    m["RefitRuleText"] = ({"equivalent-null even with the variability of refitting":
                           "so the verdict is \\emph{equivalent-null} also when the variability of refitting is included",
                           "the fixed-model interval understated the uncertainty":
                           "so the fixed-model interval understated the uncertainty"}.get(rule) or missing("C2 rule"))
    m["RefitDiscussText"] = ({"equivalent-null even with the variability of refitting":
                              "and the verdict holds when the variability of refitting is carried into the interval, which "
                              "widens that margin to \\AchievedBoundRefit",
                              "the fixed-model interval understated the uncertainty":
                              "although carrying the variability of refitting into the interval widens it beyond that"}
                             .get(rule) or missing("C2 rule"))
    c3 = d("C3")
    if c3:
        dec = sorted(c3["decisions"].items(), key=lambda kv: -kv[1])
        m["CThreeLabels"] = (f"every replicate was labelled \\emph{{{dec[0][0]}}}" if len(dec) == 1 else
                             "the labels were " + ", ".join(f"\\emph{{{k}}} {v}" for k, v in dec))
        m["CThreeRMean"] = rv(c3["R_mean"])
        m["CThreeOR"] = f"{c3['odds_ratio']:.2f}"  # 注入したオッズ比（B-6 の基礎を固定した余白の換算、10 回）
        m |= {"CThreeK": str(c3["n_equivalent"]), "CThreeN": str(c3["reps"]), "CThreeLo": pct1(c3["wilson95"][0]),
              "CThreeHi": pct1(c3["wilson95"][1]),
              "CThreeRule": ("so the fixed-model interval is anti-conservative at the margin"
                             if c3["wilson95"][0] > 0.05 else "so there is no evidence that the fixed-model interval is "
                             "anti-conservative at the margin")}
    else:
        m |= {k: missing("C3") for k in ("CThreeK", "CThreeN", "CThreeLo", "CThreeHi", "CThreeRule", "CThreeLabels",
                                         "CThreeRMean", "CThreeOR")}
    m |= _cmp_macros("GswCell", d("C4").get("gsw_cell"))
    m["GswCellRule"] = ({"assignment is the bottleneck for GSW too": "assignment is the bottleneck for surface water too",
                         "recurrent water carries little about mapped extent once terrain is known, however assigned":
                         "recurrent surface water carries little about mapped extent once terrain is known, however it "
                         "is assigned",
                         "neither reading of the rule applies": "neither reading fixed in advance applies"}
                        .get(d("C4").get("rule")) or missing("C4 rule"))
    m |= _cmp_macros("LexClean", d("C5").get("lexclean"))
    b1 = v59.get("B1") if isinstance(v59.get("B1"), dict) else {}
    bd = b1.get("bounds") or {}

    def orx(e: dict | None, what: str) -> dict[str, str]:  # オッズ比の文字列と、モンテカルロの 95% 区間
        if not isinstance(e, dict):
            return {"": missing(what), "Lo": missing(what), "Hi": missing(what)}
        return {"": _gt(e["odds_ratio_text"]), "Lo": _gt(e["mc95_text"][0]), "Hi": _gt(e["mc95_text"][1])}

    for ref, key in (("b7_recomputed", "BoundNew"), ("zero", "BoundZero"), ("b7_centred", "BoundCent"),
                     ("b5", "BoundBfive")):
        for form, tex in (("A", "A"), ("B", "B"), ("CA", "CA"), ("CB", "CB")):
            m |= {f"{key}{tex}{s}": v for s, v in orx((bd.get(ref) or {}).get(form), f"B1 {ref} {form}").items()}
    for form in ("A", "B"):
        m |= {f"MdeNew{form}{s}": v for s, v in orx((b1.get("mde80") or {}).get(form), f"B1 MDE {form}").items()}
    reps = b1.get("reps") or {}
    m["PlasRepsA"] = str(reps["A"]) if "A" in reps else missing("B1 reps")
    m["PlasRepsB"] = str(reps["B"]) if "B" in reps else missing("B1 reps")
    return m


def _gt(sx: str) -> str:
    """「> 20」のような格子の外の値を本文で組める形に。"""
    return sx.replace("> ", "\\ensuremath{>}\\,") if sx.startswith("> ") else sx


# %% B-8 の補足の表（量は全て補足に、B-8 の報告の規則）
REF_LABELS = (("b7_recomputed", "reassigned names, composition rebuilt (amendment 7; primary)"),
              ("zero", "zero"), ("b7_centred", "reassigned names centred within AOIs (amendment 7)"),
              ("b5", "reassignments of amendment 5"))
CLASS_LABELS = (("water", "Water bodies"), ("wet", "Wetland"), ("riverine", "Riverine ground"),
                ("crossing", "Crossings"), ("low", "Lowland"), ("high", "Highland"))


def _orci(d: dict | None, key: str = "pairs_ci95", what: str = "") -> str:
    if not isinstance(d, dict) or d.get("odds_ratio") is None or d.get(key) is None:
        return f"\\multicolumn{{2}}{{l}}{{{missing(what or 'odds ratio')}}}"
    return f"{orv(d['odds_ratio'])} & {orv(d[key][0])}--{orv(d[key][1])}"


def _cnt(d: dict | None) -> str:
    if not isinstance(d, dict) or d.get("n_group_flooded") is None:
        return "& & "
    return f"{_thousands(d['n_group_flooded'])} & {d.get('n_group_units', '')} & {d.get('n_aoi', '')}"


def _orm(e: dict | None) -> str:
    """B1 のオッズ比とモンテカルロの 95% 区間。"""
    if not isinstance(e, dict):
        return "---"
    lo, hi = (_gt(x) for x in e["mc95_text"])
    return f"{_gt(e['odds_ratio_text'])} ({lo}--{hi})"


def table_eight_verdict(v57: dict, v59: dict) -> str:
    c3or = (f"{v59['C3']['odds_ratio']:.2f}" if isinstance(v59.get("C3"), dict) else missing("C3"))
    b2, b3, b4 = v57.get("B2") or {}, v57.get("B3") or {}, v57.get("B4") or {}
    rows = []
    for key, lab in (("world_v2", "frozen code (verdict)"), ("world_v2_lexfix", "corrected reading")):
        e = b2.get(key)
        rows.append(f"Achieved equivalence bound, {lab}: $R$; Brier units & "
                    + (f"{rv(e['achieved_bound_R'])}; {sci(e['achieved_bound_brier'])}" if e else missing("B2")) + r" \\")
    mx = b2.get("max_over_combos")
    rows.append(r"\quad largest over five yardsticks and two bootstraps (frozen code) & "
                + (f"{rv(mx['achieved_bound_R'])} ({mx['which'].replace('_', '-').replace('/', ', ')})" if mx
                   else missing("B2 max")) + r" \\")
    for key, lab in (("total_M1s_M3s", "total increment"), ("flexible_M1_M3", "flexible models refitted"),
                     ("near_500m_wetside_not_planned", "within 500\\,m of water-side names (not planned)")):
        e = b3.get(key)
        if not e:
            rows.append(f"Location groups needed, {lab} & {missing('B3')} \\\\")
            continue
        eq = "---" if not np.isfinite(e["ratio_equivalence"]) else f"{e['ratio_equivalence']:.1f}"
        rows.append(f"Location groups needed, {lab}: equivalence; positive (multiple of {e['n_units']}) & "
                    f"{eq}; {e['ratio_positive']:.1f} \\\\")
    if b4:
        rows.append(r"Permutations with rebuilt composition: one-sided $p$; mean and 95th percentile of $R$ & "
                    f"{b4['p_one_sided']:.2f}; {rv(b4['R_perm_mean'])}, {rv(b4['R_perm_q95'])} \\\\")
        rows.append(r"\quad label with the permutation distribution as the placebo & unchanged (equivalent-null) \\")
    c2 = ((v59.get("C2") or {}) if isinstance(v59.get("C2"), dict) else {}).get("corrected")
    if isinstance(c2, dict):
        for key, lab in (("refit", "second stage refitted in each replicate"), ("fixed", "same replicates, models fixed")):
            c = c2[key]
            rows.append(f"Corrected run, {lab}: 90\\,\\% (95\\,\\%) interval of $R$ & "
                        f"{rv(c['ci90_R'][0])} to {rv(c['ci90_R'][1])} ({rv(c['ci95_R'][0])} to {rv(c['ci95_R'][1])}) \\\\")
        rows.append(f"\\quad standard deviation of $R$, refitted over fixed & {c2['sd_ratio_refit_to_fixed']:.2f} \\\\")
    else:
        rows.append(f"Second stage refitted in the bootstrap & {missing('C2')} \\\\")
    c3 = v59.get("C3")
    if isinstance(c3, dict):
        dec = ", ".join(f"{k} {v}" for k, v in sorted(c3["decisions"].items(), key=lambda kv: -kv[1]))
        rows.append(f"Effect at the margin injected, every model refitted: replicates judged equivalent "
                    f"(Wilson 95\\,\\%) & {c3['n_equivalent']} of {c3['reps']} ({pct1(c3['wilson95'][0])}--"
                    f"{pct1(c3['wilson95'][1])}) \\\\")
        rows.append(f"\\quad labels; mean $R$ & {dec}; {rv(c3['R_mean'])} \\\\")
    else:
        rows.append(f"Effect at the margin injected, every model refitted & {missing('C3')} \\\\")
    cap = ("The verdict in more detail (exploratory, amendment 8). Achieved bound: the smallest symmetric margin within "
           "which the 90\\,\\% interval of $R$ lies. Location groups needed: the multiple of the present number of groups "
           "at which the 90\\,\\% interval would lie within $\\pm0.05$ (equivalence) or the lower 95\\,\\% limit would "
           "exceed zero (positive), if the half-width shrinks with the square root of the number of groups and the point "
           "estimate stays; --- where the point estimate lies outside the margins. These multiples are power "
           "calculations at the observed point estimate and add nothing to the intervals; the paper no longer uses "
           "them. Permutations: 200 reassignments of "
           "names among the settlements of each AOI with the regional composition rebuilt. Refitted interval: 500 "
           "bootstrap replicates of location groups in which the second stage of both models is refitted with the "
           "replicate's weights, against the same replicates with the fitted models held fixed. Effect at the margin: "
           f"an odds ratio of {c3or} for all ground nearest a water-side name injected into the labels, every model "
           "refitted, and each replicate judged by the protocol's decision table without the placebo condition.")
    return _tbl("tab:eight-verdict", cap, "p{11.5cm}l", "Quantity & Value", rows)


def table_eight_bounds(v59: dict) -> str:
    b1 = v59.get("B1") if isinstance(v59.get("B1"), dict) else None
    if not b1:
        return _tbl("tab:eight-bounds", "Bounds (amendment 8).", "l", "", [missing("B1") + r" \\"])
    forms = ("A", "B", "CA", "CB")
    rows = []
    for ref, lab in REF_LABELS:
        rows.append(f"Upper bound; reference: {lab} & " + " & ".join(_orm(b1["bounds"][ref].get(f)) for f in forms)
                    + r" \\")
    rows.append(r"Detectable with 80\,\% power; reference: rebuilt composition & "
                + " & ".join(_orm(b1["mde80"].get(f)) for f in forms) + r" \\")
    rows.append(r"Odds ratio at the margin & " + " & ".join(_orm(b1["margin_or"].get(f)) for f in forms) + r" \\")
    rows.append(r"Upper bound from the first 10 replicates; reference: rebuilt composition & --- & "
                + _orm(b1["bounds_10reps"]["b7_recomputed"].get("B")) + " & --- & "
                + _orm(b1["bounds_10reps"]["b7_recomputed"].get("CB")) + r" \\")
    reps = b1.get("reps") or {}
    cap = ("Odds ratios of flooding from effects injected into the observed labels, read with different no-effect "
           "references (exploratory, amendment 8). A: all ground nearest a water-side name; B: the part within "
           "500\\,m of the settlement; CA, CB: the same for a latent indicator of truly water-related names, drawn with "
           "the measured misreading rates. Base models fixed; replicates: "
           + ", ".join(f"{f} {reps.get(f, '?')}" for f in forms)
           + ". Upper bound: the odds ratio at which the reference plus the mean injected gain reaches the upper 90\\,\\% "
           "limit of the observed increment; in parentheses, the 2.5th and 97.5th percentiles over 2,000 resamples of "
           "the replicates (Monte Carlo interval); values above the simulated range (odds ratio "
           f"{b1['or_grid_max']:.0f}) are marked $>$. The detectable odds ratio uses the normal approximation of "
           "Appendix~C. Bounds with every model refitted are not given: that simulation (amendment 7) has no effect "
           "below an odds ratio of 1.65, and reading it on a straight line from zero across the convex gain curve "
           "would understate them; its gains were 0.91--0.95 of those with the base models fixed, so refitted bounds "
           "would be slightly larger (amendment 8, addendum 9).")
    return _tbl("tab:eight-bounds", cap, "p{6.4cm}llll", "Quantity & A & B & CA & CB", rows)


def table_eight_near(v57: dict, v41: dict | None, v61: dict | None = None) -> str:
    b5, b6, b7, b12 = (v57.get(k) or {} for k in ("B5", "B6", "B7", "B12"))
    v61 = v61 or {}
    rows = [r"\multicolumn{6}{l}{\emph{Pooled water-side class, EMS labels (pairs cluster bootstrap)}} \\"]
    for key, lab in (("main", "All AOIs"), ("typed_aois", "AOIs with a recorded flood type"),
                     ("without_no_hydro", "Without AOIs lacking a hydrography layer"),
                     ("without_oceania_americas", "Without Oceania and the Americas"),
                     ("without_australia", "Without Australia"),
                     ("cleaned_reading", "Cleaned reading"), ("uncapped", "Every evaluable cell (no cap)")):
        d = b6.get(key)
        rows.append(f"{lab} & {_orci(d, what=f'B6 {key}')} & {_cnt(d)} \\\\")
    for d, lab in ((v61.get("without_pakistan"), "Without AOIs in Pakistan"),
                   (v61.get("without_pakistan_bangladesh"), "Without AOIs in Pakistan and Bangladesh")):
        rows.append(f"{lab} & {_orci(d, what=lab)} & {_cnt(d)} \\\\")
    mc = b6.get("misclassification_corrected")
    rows.append(r"Truly water-related names (misclassification corrected) & "
                + (f"{orv(mc['median'])} & {orv(mc['ci95'][0])}--{orv(mc['ci95'][1])}" if mc
                   else f"\\multicolumn{{2}}{{l}}{{{missing('B6 misclassification')}}}") + r" & & & \\")
    rs = b7.get("odds_by_river_size") or {}
    rows.append(r"\multicolumn{6}{l}{\emph{By the river nearest the settlement, and adjusted for its height "
                r"(cluster-robust)}} \\")
    for key, lab in (("lt100", "Draining less than 100\\,km$^2$"), ("100_1000", "Draining 100--1000\\,km$^2$"),
                     ("ge1000", "Draining at least 1000\\,km$^2$")):
        rows.append(f"{lab} & {_orci(rs.get(key), 'ci95', f'B7 {key}')} & {_cnt(rs.get(key))} \\\\")
    for key, lab in (("odds_same_cells_unadjusted", "Cells with a height above the main river, unadjusted"),
                     ("odds_adjusted_h_main", "\\quad adjusted for the height above the river")):
        rows.append(f"{lab} & {_orci(b7.get(key), 'ci95', f'B7 {key}')} & {_cnt(b7.get(key))} \\\\")
    rows.append(f"Names containing a nearby river or lake name & {_orci(b12.get('near_500m'), 'ci95', 'B12')} & "
                f"{_cnt(b12.get('near_500m'))} \\\\")
    rows.append(r"\multicolumn{6}{l}{\emph{By class, all classes in one model (Bonferroni 95\,\% interval; "
                r"Holm-adjusted $p$ in the last column)}} \\")
    ca = b6.get("class_adjusted") or {}
    for c, lab in CLASS_LABELS:
        d = ca.get(c)
        rows.append(f"{lab} & {_orci(d, 'bonferroni_ci95', f'B6 class {c}')} & & & "
                    + (f"{d['p_holm']:.2f}" if d else "") + r" \\")
    fr = b6.get("fractional_all_bands") or {}
    rows.append(r"\multicolumn{6}{l}{\emph{Riverine flood cover as a fractional outcome, water-side class by distance "
                r"band}} \\")
    for key, lab in (("le500m", "Within 500\\,m"), ("500m_1km", "500\\,m--1\\,km"), ("gt1km", "Beyond 1\\,km")):
        rows.append(f"{lab} & {_orci(fr.get(key), 'ci95_odds_ratio', f'B6 fractional {key}')} & & & \\\\")
    sc = b5.get("same_cells") or {}
    rows.append(r"\multicolumn{6}{l}{\emph{The two label sets on the same cells (pairs bootstrap; jackknife "
                r"with $t$ quantiles)}} \\")
    for key, lab in (("ems", "EMS"), ("gfm", "GFM")):
        d = sc.get(key)
        jk = (f"; jackknife {orv(d['jackknife_t_ci95'][0])}--{orv(d['jackknife_t_ci95'][1])}" if d else "")
        rows.append(f"{lab} & {_orci(d, what=f'B5 {key}')} & \\multicolumn{{3}}{{l}}{{{sc.get('n_units', '')} groups, "
                    f"{sc.get('n_aoi', '')} AOIs{jk}}} \\\\")
    df_ = sc.get("log_or_difference_gfm_minus_ems")
    rows.append(r"Ratio of the odds ratios, GFM over EMS & "
                + (f"{orv(df_['ratio_of_odds_ratios'])} & {orv(df_['pairs_ci95_ratio'][0])}--"
                   f"{orv(df_['pairs_ci95_ratio'][1])}" if df_ else f"\\multicolumn{{2}}{{l}}{{{missing('B5')}}}")
                + r" & & & \\")
    fp = b5.get("floodplain_same_cells") or {}
    for key, lab in (("ems", "Floodplain reading, EMS"), ("gfm", "Floodplain reading, GFM")):
        rows.append(f"{lab} & {_orci(fp.get(key), 'ci95', f'B5 floodplain {key}')} & {_cnt(fp.get(key))} \\\\")
    ind = (v41 or {}).get("sensitivity_independent_sensor") or {}
    for key, lab in (("odds_gfm", "EMS map not from Sentinel-1: GFM"), ("odds_ems", "EMS map not from Sentinel-1: EMS")):
        t = ((((ind.get(key) or {}).get("wetside_adjusted_M1") or {}).get("terms")) or {}).get("wetside:le500m")
        n_in, n_all = (ind.get(key) or {}).get("n_aoi", ""), ind.get("n_aoi", "")
        rows.append(f"{lab} & {_orci(t, 'ci95_odds_ratio', key)} & & & {n_in} of {n_all}$^{{s}}$ \\\\")
    cl = b6.get("cleaned_reading")
    cap = ("Odds ratios of flooding within 500\\,m of settlements, given the held-out terrain prediction (AOI fixed "
           "effects, distance bands, logit of the terrain model's prediction; exploratory, amendment 8, all not planned "
           "in the protocol). Pairs: cluster bootstrap over location groups with the regression refitted, 2,000 "
           "replicates. Cleaned reading: without two-letter affixes and the English suffixes that also read transferred "
           "names"
           + (f" ({pct1(cl['share_cells_changed'])} of the cells change their reading)" if cl else "")
           + ". Uncapped: every evaluable cell of the analysed AOIs, with the predictions of the uncapped sensitivity "
           "run. Misclassification: probabilistic bias analysis with the rates measured against recorded etymologies "
           "by the model that compiled the lexicon. $^{s}$AOIs that enter the model, which needs flooded and dry cells "
           "in an AOI, of the AOIs whose EMS map was not drawn from Sentinel-1. "
           "Columns 4--6: flooded cells within 500\\,m of water-side names, "
           "location groups with such cells, and AOIs.")
    return _tbl("tab:eight-near", cap, "p{6.6cm}llrrr",
                r"Estimate & Odds ratio & 95\,\% interval & Flooded & Groups & AOIs", rows)


def table_eight_heights(v57: dict) -> str:
    b7, b12 = v57.get("B7") or {}, v57.get("B12") or {}
    h1, h2, lu = (b7.get(k) or {} for k in ("h_above_100", "h_above_1000", "log_upland"))

    def cell(d: dict | None, f=_m1) -> str:
        if not isinstance(d, dict):
            return "---"
        return f"{f(d['estimate'])} ({f(d['ci95'][0])} to {f(d['ci95'][1])})"

    lg = lambda x: f"\\ensuremath{{{x:.2f}}}"
    rows = []
    for key, lab in (("wetside", "Water-side (pooled)"), ("floodplain", "Floodplain reading (riverine, wetland)")):
        d = h1.get(key) or {}
        rows.append(f"{lab} & {cell(h1.get(key))} & {cell(h2.get(key))} & {cell(lu.get(key), lg)} & "
                    f"{d.get('n_settlements', '')} \\\\")
    for c, lab in CLASS_LABELS[:5]:
        d = (h1.get("by_class") or {}).get(c) or {}
        rows.append(f"\\quad {lab} & {cell((h1.get('by_class') or {}).get(c))} & "
                    f"{cell((h2.get('by_class') or {}).get(c))} & {cell((lu.get('by_class') or {}).get(c), lg)} & "
                    f"{d.get('n_settlements', '')} \\\\")
    pc = b12.get("positive_control_hand") or {}
    rows.append(f"Names containing a nearby river or lake name (HAND of the site) & {cell(pc)} & --- & --- & "
                f"{pc.get('n_settlements', '')} \\\\")
    cap = ("Where settlements with water-side names stand (exploratory, amendment 8; no flood labels). Within-AOI "
           "difference between settlements whose name is in the class and other settlements of the same AOI in the "
           "height of the settlement's cell above the nearest point of a river draining at least 100 or 1000\\,km$^2$ "
           "(HydroRIVERS; mean elevation of the settlement's cell minus the lowest elevation of the river's cell, m), and "
           "in the logarithm of the drainage area of the nearest river; 95\\,\\% cluster-bootstrap intervals over "
           "location groups. Last row: the HAND of the settlement's site, for names that contain the name of a river, "
           "stream or lake within 5\\,km in the gazetteer.")
    return _tbl("tab:eight-heights", cap, "p{5.2cm}llll",
                r"Names & Above rivers $\geq$100\,km$^2$ & $\geq$1000\,km$^2$ & log drainage area & Settlements", rows)


SPLIT_LABELS = (("primary", "Primary"), ("country", "Countries held out"), ("temporal", "Floods from 2024"),
                ("llmcut", "After the LM cut-off"), ("activation_frozen", "Events held out whole, frozen code"),
                ("activation_corrected", "Events held out whole, corrected code"), ("gfm_refit", "GFM labels, refitted"))


def table_eight_calib(v57: dict) -> str:
    b8, b9 = v57.get("B8") or {}, v57.get("B9") or {}
    rows = []
    for key, lab in SPLIT_LABELS:
        s = b8.get(key)
        if not s:
            rows.append(f"{lab} & \\multicolumn{{7}}{{l}}{{{missing(f'B8 {key}')}}} \\\\")
            continue
        for j, mdl in enumerate(("M1", "M1s", "M1c", "M3n")):
            e = s[mdl]
            rows.append(f"{lab if j == 0 else ''} & {mdl} & {e['intercept']:.2f} & {e['slope']:.2f} & {e['brier']:.4f} & "
                        f"{sci(e['reliability'])} & {sci(e['resolution'])} & {e['auc_within_median']:.2f} \\\\")
    rows.append(r"\middlehline")
    rows.append(r"\multicolumn{8}{l}{\emph{Skill in flood-mapping terms, primary split, median over AOIs "
                r"(interquartile range)}} \\")
    for key, lab in (("M1", "Terrain model (M1)"), ("HAND", "HAND alone")):
        e = b9.get(key)
        if not e:
            continue
        q = lambda d: f"{d['median']:.2f} ({d['q25_q75'][0]:.2f}--{d['q25_q75'][1]:.2f})"
        rows.append(f"\\multicolumn{{2}}{{l}}{{{lab}}} & \\multicolumn{{6}}{{l}}{{CSI {q(e['csi'])}; hit rate "
                    f"{q(e['hit_rate'])}; false alarms {q(e['false_alarm_ratio'])}; AUC {q(e['auc_within'])}; "
                    f"{e['n_aoi']} AOIs}} \\\\")
    cap = ("Calibration of the held-out predictions (exploratory, amendment 8): intercept and slope of a logistic "
           "regression of the labels on the logit of the prediction, Brier score with the reliability and resolution "
           "terms of Murphy's decomposition (ten equal-count bins), and median within-AOI area under the ROC curve. "
           "M1: terrain model; M1s: recalibrated; M1c: with the regional composition of names (reference of the primary "
           "comparison); M3n: with the names. Skill: each AOI predicted flooded on as many cells as the labels flood "
           "(CSI: critical success index).")
    return _tbl("tab:eight-calib", cap, "llrrrrrr",
                r"Split & Model & Intercept & Slope & Brier & Reliability & Resolution & AUC", rows)


def table_eight_audit(v57: dict) -> str:
    b11 = v57.get("B11") or {}
    wf = b11.get("watched_forms") or {}
    rows = []
    for key in (*AUDIT_SHORT, *AUDIT_EN, "hi:char:word", "fa:bela:word", "fa:kacha:word", "fa:pattan:word",
                "fa:pattan:prefix", "fa:tibba:word"):
        e = wf.get(key)
        if not e:
            continue
        blk, form, pos = key.split(":")
        cc = ", ".join(f"{c} {n}" for c, n in sorted(e["by_country"].items(), key=lambda kv: -kv[1]))
        ex = ", ".join(k.title() for k in list(e["examples"])[:4])
        rows.append(f"{blk} & \\emph{{{form}}} ({pos}) & {e['settlements']} & {_thousands(e['cells'])} & {cc} & {ex} \\\\")
    rows.append(r"\middlehline")
    rows.append(r"\multicolumn{6}{l}{\emph{Most frequent matched forms by country (class; settlements / cells)}} \\")
    for cc, lst in (b11.get("top_forms_by_country") or {}).items():
        sx = "; ".join(f"\\emph{{{e['form']}}} ({e['cat']}; {e['settlements']}/{_thousands(e['cells'])})" for e in lst[:5])
        rows.append(f"{cc} & \\multicolumn{{5}}{{p{{14cm}}}}{{{sx}}} \\\\")
    cap = ("Audit of matched forms (exploratory, amendment 8): forms that the lexicon's own rules excluded (two-letter "
           "affixes), English suffixes that also read transferred and personal names, and regional words, with the "
           f"number of analysed settlements ({_thousands(b11.get('n_settlements', 0))} in all) whose name they match, the "
           "cells nearest those settlements, the countries and examples (most cells first).")
    return _tbl("tab:eight-audit", cap, "llrrp{3.2cm}p{5.6cm}",
                r"Block & Form & Settlements & Cells & Countries & Examples", rows)


def table_eight_coder(v57b: dict | None) -> str:
    if not v57b:
        return _tbl("tab:eight-coder", "Second coder (amendment 8).", "l", "", [missing("B13") + r" \\"])
    pr, c1 = v57b["coder2_precision_recall"], v57b["coder1_precision_recall_same_rows"]
    rows = [f"Excerpts coded by both; stating a derivation, agreement and $\\kappa$ & {v57b['n_items']}; "
            f"{_pct0(v57b['status_agreement'])}, {v57b['status_kappa']:.2f} \\\\",
            f"Water-related derivation, agreement and $\\kappa$ ({v57b['n_both_coded']} excerpts) & "
            f"{_pct0(v57b['water_related_agreement'])}, {v57b['water_related_kappa']:.2f} \\\\",
            "$\\kappa$ by class & " + ", ".join(f"{c} {v:.2f}" for c, v in v57b["category_kappa"].items()) + r" \\",
            f"Lexicon, broad: precision; recall (coder 1 on the same excerpts; coder 2) & "
            f"{_pct0(c1['lex_wetside']['precision'])}; {_pct0(c1['lex_wetside']['recall'])} --- "
            f"{_pct0(pr['broad']['lex_wetside']['precision'])}; {_pct0(pr['broad']['lex_wetside']['recall'])} \\\\",
            f"Language model, broad: precision; recall (coder 1; coder 2) & "
            f"{_pct0(c1['nread_wetside']['precision'])}; {_pct0(c1['nread_wetside']['recall'])} --- "
            f"{_pct0(pr['broad']['nread_wetside']['precision'])}; {_pct0(pr['broad']['nread_wetside']['recall'])} \\\\",
            f"Lexicon, strict: precision; recall (coder 2) & {_pct0(pr['strict']['strict_flag']['precision'])}; "
            f"{_pct0(pr['strict']['strict_flag']['recall'])} \\\\"]
    tr, st = v57b["transferred_among_lexicon_wetside_true"], v57b["at_site_no_among_lexicon_wetside_true"]
    rows.append(f"Names read as water-side with a water-related derivation: transferred or commemorative; feature not at "
                f"the site & {tr['yes']} of {tr['n']}; {st['no']} of {st['n']} \\\\")
    cap = ("A second coder of the Wikipedia excerpts (exploratory, amendment 8): a language model of another developer "
           "that took no part in the lexicon or the first coding, given the same guide and blind to the readings. "
           "Precision and recall weight the two sampling strata by their size. The rule of amendment 8: if either rate "
           "differed by 10 percentage points or more between the coders, the paper would say that the measured rates "
           "depend on the coder; they did not. Agreement between two models is not a validation by toponymists.")
    return _tbl("tab:eight-coder", cap, "p{9.6cm}l", "Quantity & Value", rows)


def table_eight_prompt() -> str:
    """言語モデルの読み（N-read）の固定した prompt（system と user の雛形）をそのまま。"""
    from topsem.world.nread import _USER, SYSTEM, prompt_sha256

    def tt(sx: str) -> str:
        sx = (sx.replace("\\", "\\textbackslash{}").replace("{", "\\{").replace("}", "\\}")
              .replace("_", "\\_").replace("<", "\\textless{}").replace(">", "\\textgreater{}"))
        return "\\newline ".join(sx.split("\n"))

    rows = [f"System & \\texttt{{{tt(SYSTEM)}}} \\\\", f"User & \\texttt{{{tt(_USER)}}} \\\\"]
    cap = ("The fixed prompt of the language-model reading (N-read): the system message and the template of the user "
           "message, in which the name, its local spelling where different and the country are filled in for each "
           "name (doubled braces stand for single braces in the text sent). Temperature 0; sha256 of the two strings "
           f"\\texttt{{{prompt_sha256()[:16]}}}\\ldots")
    return _tbl("tab:eight-prompt", cap, "lp{14cm}", "Message & Text", rows)


def table_eight_where(v60: dict | None, v61: dict | None = None) -> str:
    """水系側の名の集落から 500 m 以内の浸水セル（EMS）の国・事象ごとの数（編集者の求め、記述）。"""
    if not v60:
        return _tbl("tab:eight-where", "Where the flooded cells near water-side names lie.", "l", "",
                    [missing("v60") + r" \\"])
    rows = [r"\multicolumn{5}{l}{\emph{By country of the nearest settlement}} \\"]
    for r in v60["by_country"][:15]:
        rows.append(f"{r['cc']} & {r['flooded_near']} & {_thousands(r['cells_near'])} & {pct1(r['share_of_flooded_near'])} & "
                    f"{r['aois']} \\\\")
    rest = v60["by_country"][15:]
    if rest:
        rows.append(f"Other {len(rest)} countries & {sum(r['flooded_near'] for r in rest)} & "
                    f"{_thousands(sum(r['cells_near'] for r in rest))} & "
                    f"{pct1(sum(r['share_of_flooded_near'] for r in rest))} & {sum(r['aois'] for r in rest)} \\\\")
    rows.append(r"\multicolumn{5}{l}{\emph{The activations with the most such cells}} \\")
    for r in v60["by_activation"][:10]:
        lab = f"{r['activation']}: {_tex(r.get('name', ''))[:48]} ({r.get('year', '')})"
        rows.append(f"{lab} & {r['flooded_near']} & {_thousands(r['cells_near'])} & {pct1(r['share_of_flooded_near'])} & "
                    f"{r['aois']} \\\\")
    pk = (v61 or {}).get("pakistan_by_activation") or []
    if pk:  # 編集者の細かい指摘 17: パキスタンの層を支配する activation（追記 6 の c）
        rows.append(r"\multicolumn{5}{l}{\emph{AOIs in Pakistan, by activation (evaluable cells; flooded near / near)}} \\")
        for r in pk:
            rows.append(f"{r['activation']}: {_tex(r['name'])[:48]} ({r['year']}) & {r['flooded_near']} & "
                        f"{_thousands(r['cells_near'])} & {_thousands(r['flooded'])} of {_thousands(r['cells'])} & "
                        f"{r['aois']} \\\\")
    cap = ("Where the flooded cells within 500\\,m of settlements with water-side names lie (EMS labels; descriptive): "
           f"{v60['flooded_near_total']} flooded cells among {_thousands(v60['cells_near_total'])} such cells, in the AOIs "
           "of the odds ratios near settlements (AOIs with both flooded and dry cells). Share: of all flooded cells near "
           "water-side names. For the AOIs in Pakistan the fourth column gives the flooded cells of all evaluable cells "
           "(amendment 8, addendum 6).")
    return _tbl("tab:eight-where", cap, "p{8.4cm}rrrr", r"Country or activation & Flooded & Cells & Share & AOIs", rows)


def supp_tables_b8(v57: dict | None, v57b: dict | None, v59: dict | None, v41: dict | None,
                   v60: dict | None = None, v61: dict | None = None) -> None:
    v57, v59 = v57 or {}, v59 or {}
    (HERE / "app_table_eight_prompt.tex").write_text(table_eight_prompt(), encoding="utf-8")
    (HERE / "app_table_eight_where.tex").write_text(table_eight_where(v60, v61), encoding="utf-8")
    for name, tex in (("eight_verdict", table_eight_verdict(v57, v59)), ("eight_bounds", table_eight_bounds(v59)),
                      ("eight_near", table_eight_near(v57, v41, v61)), ("eight_heights", table_eight_heights(v57)),
                      ("eight_calib", table_eight_calib(v57)), ("eight_audit", table_eight_audit(v57)),
                      ("eight_coder", table_eight_coder(v57b))):
        (HERE / f"app_table_{name}.tex").write_text(tex, encoding="utf-8")


# %% B-9 の補足の表（初見の査読 第 5 回への事後の解析、量は全て補足に）
def _lim(e: dict | None, what: str) -> str:
    """v64 の上限・下限（オッズ比とモンテカルロの 95% 区間）。無ければ未計算の印。"""
    return _orm(e) if isinstance(e, dict) else missing(what)


def table_nine_bounds(v64: dict | None) -> str:
    """B-9 の B1・B2・C2・C3: 区間の数え方ごとの上限、下限、集落ごとに違う効果、物差しの区間での余白の換算。"""
    v = v64 or {}
    b = v.get("B1_B2") if isinstance(v.get("B1_B2"), dict) else None
    if not b:
        return _tbl("tab:nine-bounds", "Bounds (amendment 9).", "l", "", [missing("v64 B1") + r" \\"])
    forms = ("A", "B", "CA", "CB")
    bd = b["bounds"]
    rows = [r"\multicolumn{5}{l}{\emph{Upper bounds}} \\"]
    for kind, lab in (("fixed", "Base models fixed"), ("fold", "Spread over fold assignments added"),
                      ("refit", "Second stage refitted within the bootstrap")):
        rows.append(f"\\quad {lab} & " + " & ".join(_lim(bd[kind]["b7_recomputed"].get(f), "v64 B1")
                                                     for f in forms) + r" \\")
    for ref, lab in REF_LABELS[1:]:
        rows.append(f"\\quad Base models fixed; reference: {lab} & "
                    + " & ".join(_lim(bd["fixed"][ref].get(f), "v64 B1") for f in forms) + r" \\")
    rows.append(r"\multicolumn{5}{l}{\emph{Detectable with 80\,\% power}} \\")
    for kind, lab in (("fixed", "Interval with the base models fixed"), ("refit", "Interval with the second stage refitted")):
        rows.append(f"\\quad {lab} & " + " & ".join(_lim(b["mde80"][kind].get(f), "v64 B1") for f in forms) + r" \\")
    my = b["margin_by_yardstick"]
    rows.append(r"\multicolumn{5}{l}{\emph{Odds ratio at the margin, terrain gain at}} \\")
    for q, lab in (("point", "its estimate"), ("q025", "the 2.5th percentile of its 95\\,\\% interval"),
                   ("q975", "the 97.5th percentile")):
        rows.append(f"\\quad {lab} & {_lim(my.get(q), 'v64 B2')} & --- & --- & --- \\\\")
    c2 = v.get("C2") if isinstance(v.get("C2"), dict) else None
    rows.append(r"\multicolumn{5}{l}{\emph{Effects that remove flooding: lower bound}} \\")
    rows.append("\\quad Base models fixed & " + (" & ".join(_lim(c2["lower_bounds"]["b7_recomputed"].get(f), "v64 C2")
                                                         for f in ("A", "B")) if c2 else
                                                f"{missing('v64 C2')} & {missing('v64 C2')}") + r" & --- & --- \\")
    c3 = v.get("C3") if isinstance(v.get("C3"), dict) else None
    rows.append(r"\multicolumn{5}{l}{\emph{Effects that differ between settlements (mean odds ratio)}} \\")
    rows.append("\\quad Upper bound, base models fixed & "
                + (_lim(c3["bound"]["b7_recomputed"], "v64 C3") if c3 else missing("v64 C3")) + r" & --- & --- & --- \\")
    rows.append("\\quad Odds ratio at the margin & " + (_lim(c3["margin_or"], "v64 C3") if c3 else missing("v64 C3"))
                + r" & --- & --- & --- \\")
    grid_a = b["grid"]["A"]
    low = [x for x in grid_a if x < 0.25 - 1e-9]  # noqa: PLR2004
    low_txt = (f"; for A and CA also {', '.join(f'{np.exp(x):.2f}' for x in low)} (amendment 9)" if low else "")
    neg_txt = (f" Lower bounds: effects that set flooded labels to dry with probability $1-p'/p$, at odds ratios "
               f"{', '.join(f'{np.exp(x):.2f}' for x in sorted(c2['deltas'], reverse=True))}; the bound is the odds "
               "ratio below which the gain would exceed the same limit." if c2 else "")
    het_txt = (" Effects that differ between settlements: odds ratio $e^{\\delta w_s}$ around settlement $s$ with "
               f"weights $w_s$ drawn from an exponential distribution with mean one ({c3['reps']} replicates); "
               "bound and margin are read on the mean." if c3 else "")
    cap = ("Odds ratios of flooding that the data exclude, from effects injected into the observed labels "
           "(exploratory; amendments 8 and 9). A: all ground nearest a water-side name; B: the part within 500\\,m "
           "of the settlement; CA, CB: the same targets for truly water-related names, given the measured misreadings. "
           "Upper bound: the odds ratio at which the no-effect reference plus the mean injected gain reaches the upper "
           "90\\,\\% limit of the observed increment, with that limit taken from the bootstrap of the primary analysis "
           "(base models fixed), with the spread over 20 assignments of location groups to folds added, or with the "
           "second stage refitted within each bootstrap replicate. Reference: reassigned names with rebuilt "
           "composition unless stated. Detectable: the odds ratio whose mean gain exceeds the reference by "
           "2.49 standard errors of the increment. Simulated odds ratios: "
           f"{', '.join(f'{np.exp(x):.2f}' for x in grid_a if x >= 0.25 - 1e-9)}{low_txt}; "  # noqa: PLR2004
           f"{b['reps']['A']} replicates per effect and target; values beyond them are marked $>$ or $<$. In "
           "parentheses: 2.5th and 97.5th percentiles over 2,000 resamples of the replicates (Monte Carlo interval)."
           + neg_txt + het_txt)
    return _tbl("tab:nine-bounds", cap, "p{6.0cm}llll", "Quantity & A & B & CA & CB", rows)


def table_nine_power(v64: dict | None) -> str:
    """B-9 の C4: 全てのモデルを当て直した模擬の検出力。"""
    c4 = (v64 or {}).get("C4")
    if not isinstance(c4, dict):
        return _tbl("tab:nine-power", "Simulated power (amendment 9).", "l", "", [missing("v64 C4") + r" \\"])
    rows = []
    for x in c4["by_delta"]:
        lo, hi = x["wilson95_positive"]
        rows.append(f"{x['odds_ratio']:.2f} & {x['reps']} & {x['n_positive']} ({lo:.2f}--{hi:.2f}) & "
                    f"{x['n_equivalent']} & {rv(x['R_mean'])} \\\\")
    cap = ("Verdicts of the decision rule when effects of the stated odds ratio are injected on all ground nearest a "
           "water-side name and every model, including the terrain model, is refitted on the injected labels in each "
           "fold (exploratory; amendment 9). Positive: lower 95\\,\\% limit of $R$ above zero (the decision table "
           "without its placebo condition); in parentheses, the Wilson 95\\,\\% interval of the share. Equivalent: "
           "90\\,\\% interval of $R$ within $\\pm 0.05$. Replicates "
           f"{min(c4['rep_ids'])}--{max(c4['rep_ids'])}. {c4['rule'][:1].upper() + c4['rule'][1:]}"
           + (f": {c4['simulated_detectable_or']:.2f}." if c4.get("simulated_detectable_or") else "."))
    return _tbl("tab:nine-power", cap, "lllll",
                "Odds ratio & Replicates & Positive & Equivalent & Mean $R$", rows)


def _run_row(lab: str, c: dict | None, near: dict | None, what: str) -> str:
    p = (c or {}).get("primary") if isinstance(c, dict) else None
    if not isinstance(p, dict):
        return f"{lab} & \\multicolumn{{6}}{{l}}{{{missing(what)}}} \\\\"
    nr = (f"{orv(near['odds_ratio'])} ({orv(near['pairs_ci95'][0])}--{orv(near['pairs_ci95'][1])})"
          if isinstance(near, dict) and near.get("pairs_ci95") else "---")
    return (f"{lab} & {_n(int(c['n_cells']))} & {1e3 * c['yardstick']:.2f} & {g._r(p['R'])} & "
            f"{g._ci(p['ci90_R'])} & {p['decision']} & {nr} \\\\")


def table_nine_runs(v64: dict | None, v57: dict | None, v63: dict | None) -> str:
    """B-9 の C5〜C7: 永久水域の共通の定義、上限の抽出の種、FABDEM の地形。主の行は直したコードの実行。"""
    v = v64 or {}
    lex = g._json(DATA / "runs" / "world_v2_lexfix" / "metrics_primary.json") or {}
    prim = None
    if lex.get("primary"):
        prim = {"primary": lex["primary"] | {"decision": lex["primary"].get("decision")},
                "n_cells": lex.get("n_cells") or lex["primary"].get("n_cells"),
                "yardstick": lex["primary"]["delta_bs_M0_M1"]}
    rows = [_run_row("Primary (corrected code)",
                     prim if prim and prim.get("n_cells") else None,
                     ((v57 or {}).get("B6") or {}).get("main"), "lexfix metrics")]
    rows.append(_run_row("GSW permanent water removed", v.get("C5"), v.get("C5_near"), "v64 C5"))
    for k, run in enumerate(sorted((v.get("C6") or {}).keys()), start=1):
        rows.append(_run_row(f"Subsampling seed {k}", v["C6"][run], None, f"v64 C6 {k}"))
    rows.append(_run_row("FABDEM terrain (bare earth)", v.get("C7"), v.get("C7_near"), "v64 C7"))
    gp = ((v.get("C5") or {}).get("record") or {}).get("b6", {}).get("gsw_perm") if isinstance(v.get("C5"), dict) else None
    gsw_txt = (f" Permanent water: cells whose mean JRC GSW occurrence at 25 points is at least 50\\,\\% were removed in "
               f"every AOI ({_n(gp['n_dropped'])} cells, {pct1(gp['share_dropped'])}, in {gp['n_aoi_with_drop']} AOIs; "
               f"{_n(gp['n_dropped_with_hydro'])} in AOIs with a hydrography layer and {_n(gp['n_dropped_without_hydro'])} "
               "in AOIs without one), in addition to the water of each product's hydrography layer." if gp else "")
    sk = v.get("C7_skill") if isinstance(v.get("C7_skill"), dict) else None
    fab_txt = (" Bare-earth elevation: the elevation features rebuilt from FABDEM v1-2 (Hawker et al., 2022), and HAND "
               "replaced by the heights of each cell above the nearest vertex of a HydroRIVERS river draining at least "
               "100 and 1000\\,km$^2$; the location model M0 is unchanged.")
    b5 = (v63 or {}).get("B5") if isinstance((v63 or {}).get("B5"), dict) else None
    if sk and b5:
        fab_txt += (f" Its terrain model flagged the observed number of flooded cells per AOI with a median critical "
                    f"success index of {sk['M1']['csi']['median']:.2f} (with GLO-30 and HAND: "
                    f"{b5['M1']['csi']['median']:.2f}); with the names, the median difference was "
                    f"{sk['delta_csi_M3n_minus_M1c']['median']:.3f}.")
    cap = ("The primary comparison under the robustness runs of amendment 9 (exploratory; corrected code; the "
           "primary row uses the frozen subsampling seed 20261025). Cells: "
           "evaluated cells; terrain gain: $\\mathrm{BS}(\\mathrm{M0})-\\mathrm{BS}(\\mathrm{M1})$ in units of "
           "$10^{-3}$; odds ratio: flooding within 500\\,m of water-side names relative to other names, given the "
           "held-out terrain prediction of the run, with a pairs cluster bootstrap 95\\,\\% interval (2,000 "
           "replicates); not computed for the seeds." + gsw_txt
           + " Seeds: the seed that draws cells within location groups under the cap (frozen: 20261025)." + fab_txt
           + " By the rule fixed in amendment 9, these runs do not replace the confirmatory verdict, and an odds ratio "
           "whose interval excludes one would remain an uncorroborated description.")
    return _tbl("tab:nine-runs", cap, "lllllll",
                "Run & Cells & Terrain gain & $R$ & 90\\,\\% interval & Verdict & Odds ratio (95\\,\\%)", rows)


def table_nine_skill(v63: dict | None) -> str:
    """B-9 の B3〜B5: Murphy の分解、AOI の中だけの評価、浸水の地図の技能。"""
    v = v63 or {}
    b3, b4, b5 = v.get("B3"), v.get("B4"), v.get("B5")
    if not (isinstance(b3, dict) and isinstance(b4, dict) and isinstance(b5, dict)):
        return _tbl("tab:nine-skill", "Skill (amendment 9).", "l", "", [missing("v63 B3-B5") + r" \\"])
    d = b3["diff_M3n_minus_M1c"]
    rows = [r"\multicolumn{4}{l}{\emph{Murphy decomposition of the Brier score}} \\",
            "\\quad Brier score & " + f"{b3['M1c']['brier']:.5f} & {b3['M3n']['brier']:.5f} & {sci(d['brier'])} \\\\",
            "\\quad Reliability (lower is better) & " + f"{sci(b3['M1c']['reliability'])} & {sci(b3['M3n']['reliability'])}"
            f" & {sci(d['reliability'])} ({sci(d['reliability_ci95'][0])} to {sci(d['reliability_ci95'][1])}) \\\\",
            "\\quad Resolution (higher is better) & " + f"{sci(b3['M1c']['resolution'])} & {sci(b3['M3n']['resolution'])}"
            f" & {sci(d['resolution'])} ({sci(d['resolution_ci95'][0])} to {sci(d['resolution_ci95'][1])}) \\\\",
            f"\\quad Uncertainty & \\multicolumn{{2}}{{l}}{{{sci(b3['M1c']['uncertainty'])}}} & \\\\",
            r"\multicolumn{4}{l}{\emph{Within AOIs (each AOI's flooded fraction given)}} \\",
            f"\\quad Terrain gain within AOIs & \\multicolumn{{2}}{{l}}{{{sci(b4['yardstick_within'])}}} & \\\\",
            f"\\quad Names increment $\\mathrm{{BS}}(\\mathrm{{M1c}})-\\mathrm{{BS}}(\\mathrm{{M3n}})$ & & & "
            f"{sci(b4['delta_bs_within'])} "
            f"({sci(b4['ci95_delta_bs'][0])} to {sci(b4['ci95_delta_bs'][1])}) \\\\",
            f"\\quad $R$ (90\\,\\% interval) & & & {g._r(b4['R_within'])} {g._ci(b4['ci90_R'])} \\\\",
            r"\multicolumn{4}{l}{\emph{Maps flagging the observed number of flooded cells per AOI (median; quartiles)}} \\"]
    for key, lab in (("csi", "Critical success index"), ("hit_rate", "Hit rate"), ("false_alarm_ratio", "False-alarm ratio")):
        cells = []
        for mm in ("M1c", "M3n"):
            x = b5[mm][key]
            cells.append(f"{x['median']:.3f} ({x['q25_q75'][0]:.2f}--{x['q25_q75'][1]:.2f})")
        h = b5["HAND"][key]
        rows.append(f"\\quad {lab} & {cells[0]} & {cells[1]} & HAND: {h['median']:.3f} "
                    f"({h['q25_q75'][0]:.2f}--{h['q25_q75'][1]:.2f}) \\\\")
    dc = b5["delta_csi_M3n_minus_M1c"]
    rows.append(f"\\quad Difference in the critical success index & & & median {dc['median']:.3f} "
                f"({dc['ci95_median'][0]:.3f} to {dc['ci95_median'][1]:.3f}) \\\\")
    cap = ("The skill of the held-out predictions of the corrected run with and without the names, decomposed "
           "(exploratory; amendment 9). Murphy decomposition (Murphy, 1973) over ten bins of predicted "
           "probability; within AOIs: each model's predictions shifted on the logit scale so that its mean in each AOI "
           f"equals the observed flooded fraction ({b4['n_aoi_inner']} AOIs with both flooded and dry cells), and the "
           "comparison and yardstick recomputed; maps: in each AOI, the cells with the highest predicted probability "
           "(for HAND, the lowest HAND, ties broken by mean HAND and then at random) flagged as flooded, as many as "
           f"were observed flooded ({b5['M1c']['n_aoi']} AOIs). Differences of the decomposition and of the critical "
           "success index: M3n minus M1c; the names raised the index in "
           f"{pct1(dc['share_aoi_positive'])} of the AOIs and lowered it in {pct1(dc['share_aoi_negative'])}. M1c: "
           "without names; M3n: with names; M1 ranks cells within AOIs as M1c does. Intervals: "
           f"95\\,\\% from a cluster bootstrap over location groups ({b3['replicates']:,} replicates) unless stated.")
    return _tbl("tab:nine-skill", cap, "llll", "Quantity & M1c & M3n & Difference or reference", rows)


def table_nine_activations(v63: dict | None) -> str:
    """B-9 の B6: 浸水セルのある activation を一つずつ除いた、集落の近くのオッズ比。"""
    b6 = (v63 or {}).get("B6")
    if not isinstance(b6, dict):
        return _tbl("tab:nine-activations", "Leave one activation out (amendment 9).", "l", "",
                    [missing("v63 B6") + r" \\"])
    rs = b6["rows"]
    half = -(-len(rs) // 2)

    def cell(r: dict | None) -> str:
        if r is None:
            return "& & &"
        return (f"{r['activation']} & {r['n_flooded_near']} & {orv(r['odds_ratio'])} & "
                f"{orv(r['ci95'][0])}--{orv(r['ci95'][1])}")
    rows = [f"{cell(rs[i])} & {cell(rs[i + half] if i + half < len(rs) else None)} \\\\" for i in range(half)]
    cap = ("The pooled odds ratio of flooding within 500\\,m of water-side names (EMS labels, the main model of "
           "amendment 8) with each activation that has flooded cells near water-side names left out in turn "
           f"(exploratory; amendment 9; {b6['n_activations']} activations). Flooded: flooded cells near water-side "
           "names in the activation left out. Intervals: 95\\,\\%, cluster-robust over location groups.")
    return _tbl("tab:nine-activations", cap, "llll|llll",
                "Left out & Flooded & OR & 95\\,\\% & Left out & Flooded & OR & 95\\,\\%", rows)


def table_nine_classq(v63: dict | None) -> str:
    """B-9 の B7: 6 区分の係数の Holm の p と Benjamini–Hochberg の q。"""
    b7 = (v63 or {}).get("B7")
    if not isinstance(b7, dict):
        return _tbl("tab:nine-classq", "Class multiplicity (amendment 9).", "l", "", [missing("v63 B7") + r" \\"])
    rows = [f"{lab} & {b7[k]['p']:.3f} & {b7[k]['p_holm']:.3f} & {b7[k]['q_bh']:.3f} \\\\"
            for k, lab in CLASS_LABELS if k in b7]
    cap = ("The six class coefficients of the model with all lexicon classes within 500\\,m (EMS labels; "
           "the figure of odds ratios near settlements in the paper): unadjusted $p$, Holm-adjusted $p$ "
           "(Holm, 1979) and Benjamini--Hochberg $q$ "
           "(Benjamini and Hochberg, 1995) (exploratory; amendment 9).")
    return _tbl("tab:nine-classq", cap, "llll", "Class & $p$ & Holm $p$ & BH $q$", rows)


def table_nine_text() -> str:
    """本文にだけ出ていた量（編集者の指摘 8、第 5 回）: 本文と同じマクロで並べる（値は numbers.tex から）。"""
    rows = [
        r"\multicolumn{2}{l}{\emph{What the pooled test can detect (Sect.~4)}} \\",
        r"\quad Share of cells within 500\,m of a settlement with a water-side name & \CeilShareFiveHundred \\",
        r"\quad Largest possible $R$ from perfect prediction near water-side names: within 250\,m; 500\,m; 1\,km & "
        r"\CeilRTwoFifty; \CeilRFiveHundred; \CeilROneKm \\",
        r"\multicolumn{2}{l}{\emph{Evidence near settlements (Sect.~5)}} \\",
        r"\quad Joint Wald test, class coefficients in the three distance bands: $p$ & \WaldClassP \\",
        r"\quad Joint Wald test, water-side coefficients in the three distance bands: $p$ & \WaldWetP \\",
        r"\quad Odds ratio near water-side names without partly flooded cells (EMS) & \PartialOR{} "
        r"(\PartialORLo--\PartialORHi) \\",
        r"\quad Odds ratio adjusted for built-up and tree cover: EMS; GFM & \CciOREms{} (\CciOREmsLo--\CciOREmsHi); "
        r"\CciORGfm{} (\CciORGfmLo--\CciORGfmHi) \\",
        r"\quad Odds ratio by flood size, thirds of AOIs: smallest floods; middle; largest & \SevLowOR{} "
        r"(\SevLowORLo--\SevLowORHi; \SevLowFlooded{} flooded cells); \SevMidOR; \SevHighOR \\",
        r"\quad EMS-evaluable cells of the GFM AOIs that GFM could not assess: all; within 500\,m of water-side names; "
        r"near other settlements & \MaskShareAll; \MaskShareWet; \MaskShareOther \\",
        r"\multicolumn{2}{l}{\emph{How the names were read (Sect.~6 and Appendix~B)}} \\",
        r"\quad Settlements whose name contains the name of a river, stream or lake within 5\,km: share of the analysed "
        r"settlements; of them read as water-side by the lexicon & \HydSharePct; \HydOverlapPct \\",
        r"\quad Water namesakes recorded on Wikidata that this reading recovers & \HydRecallK{} of \WDNHydro \\",
        r"\quad Height above the drainage of such settlements minus that of other settlements of their AOI (m) & "
        r"\HydHand{} (\HydHandLo{} to \HydHandHi) \\",
        r"\quad Names of each lexicon class put in that class by the language model (range over classes) & "
        r"\NReadAgreeMin--\NReadAgreeMax \\",
        r"\quad Distinct names read by the language model & \NReadNames \\",
        r"\quad Analysed settlements with an ASCII name containing \emph{ae}, \emph{oe} or \emph{ue} & \DigraphPct \\",
        r"\quad Water-side matches gained; lost if these digraphs were read as plain vowels & \DigraphGainWet; "
        r"\DigraphLoseWet \\",
    ]
    cap = ("Quantities stated in the paper and not given in another table (exploratory unless pre-specified; "
           "amendments 6 to 8). Intervals: 95\\,\\%, cluster-robust over location groups for the odds ratios. The "
           "values are typeset from the same generated file as the paper.")
    return _tbl("tab:nine-text", cap, "p{11.5cm}l", "Quantity & Value", rows)


def supp_tables_b9(v63: dict | None, v64: dict | None, v57: dict | None) -> None:
    for name, tex in (("nine_bounds", table_nine_bounds(v64)), ("nine_power", table_nine_power(v64)),
                      ("nine_runs", table_nine_runs(v64, v57, v63)), ("nine_skill", table_nine_skill(v63)),
                      ("nine_activations", table_nine_activations(v63)), ("nine_classq", table_nine_classq(v63)),
                      ("nine_text", table_nine_text())):
        (HERE / f"app_table_{name}.tex").write_text(tex, encoding="utf-8")


def _cmp_row(label: str, amend: str, code: str, c: dict | None) -> str:
    if not isinstance(c, dict) or c.get("R") is None:
        return f"{label} & {amend} & {code} & \\multicolumn{{3}}{{l}}{{{missing(label)}}} \\\\"
    dec = c.get("decision")
    dec = dec.get("label", "---") if isinstance(dec, dict) else (dec or "---")
    return f"{label} & {amend} & {code} & {g._r(c['R'])} & {g._ci(c.get('ci90_R'))} & {dec} \\\\"


def _decision_label(c: dict) -> str:
    """判定の label（文字列、または {"label": ...} の dict）。"""
    d = c.get("decision")
    return d["label"] if isinstance(d, dict) else str(d)


def posthoc_checks(v41: dict | None, v42: dict | None, v47: dict | None, v48: dict | None, v52: dict | None,
                   v55: dict | None, v56: dict | None, v59: dict | None = None,
                   ) -> list[tuple[str, str, str, dict | None]]:
    """事後の確認（本文の図 4b と付録の表で同じ並び）: (表示, 改訂, コード, verdict の形の dict)。"""
    v47, v52, v56 = v47 or {}, v52 or {}, v56 or {}
    rows: list[tuple[str, str, str, dict | None]] = []
    act = g._json(DATA / "runs" / "world_v2_actgrp" / "metrics_primary.json")
    rows.append(("Flood events held out whole", "2", "frozen", act.get("primary") if act else None))
    rows.append(("Flood events held out whole", "7", "corrected", (v56.get("actfix") or {}).get("primary")
                 if isinstance(v56.get("actfix"), dict) else None))
    for key, lab in (("notrace", "Flood traces counted as not flooded"), ("nonwater", "Flooded share over land only")):
        rows.append((lab, "5", "corrected", (v47.get("label_variants") or {}).get(key)))
    pt = v56.get("partial")
    rows.append(("Partly flooded cells excluded", "7", "corrected", pt.get("primary") if isinstance(pt, dict) else None))
    if v48:
        rows.append(("Every AOI weighted equally", "6", "frozen", v48["aoi_weighted"]["frozen"]))
    if v55:
        rows.append(("European AOIs alone, own yardstick", "7", "corrected", v55["europe"]["europe_lexfix"]))
    if v42:
        st = dict(v42["settler_excluded"])
        rows.append(("Without the Americas and Oceania", "5", "corrected", st))
    for key, cmp_, lab, am in (("strict", "strict", "Strict reading", "6"),
                               ("gn", "primary", "Gazetteer cleaned of city sections and duplicates", "6"),
                               ("nohydro", "primary", "Without AOIs lacking a hydrography layer", "6"),
                               ("typed", "primary", "Without AOIs lacking a flood type", "6"),
                               ("nread", "H4", "Language-model reading, primary split", "6")):
        b = v52.get(key)
        rows.append((lab, am, "corrected", b.get(cmp_) if isinstance(b, dict) else None))
    fp = v56.get("floodplain")
    rows.append(("Floodplain reading (riverine and wetland classes)", "7", "corrected",
                 fp.get("floodplain") if isinstance(fp, dict) else None))
    c5 = (v59 or {}).get("C5")
    rows.append(("Cleaned lexicon (two-letter affixes and English suffixes removed)", "8", "corrected",
                 c5.get("lexclean") if isinstance(c5, dict) else None))
    if v41:
        rows.append(("GFM labels, existing predictions rescored", "3", "frozen", v41["C1_gfm"]))
    return rows


def table_checks(rows: list[tuple[str, str, str, dict | None]]) -> str:
    body = [_cmp_row(lab, am, code, c) for lab, am, code, c in rows]
    return _tbl("tab:si-checks", "Post hoc checks of the within-area comparison (M1c to M3n or its analogue): $R$ "
                "with its 90\\,\\% cluster-bootstrap interval and the verdict of the protocol's decision table. "
                "Amendment: the amendment that planned the check; code: the frozen code of the confirmatory run or "
                "the code with the corrected reading of Norwegian names. The yardstick is that of each run (for the "
                "European AOIs, that of these AOIs).", "lllrll",
                "Check & Amendment & Code & $R$ & 90\\,\\% interval & Verdict", body)


def table_posthoc_models() -> str:
    rows = [r"M3d & logit(M1) + name features centred within each AOI, without composition terms & 2 \\",
            r"M1co, M3o & as M1c and M3n with the HAND of the nearest settlement's site in place of the names & 5 \\",
            (r"M1cL, M3oL & as M1c and M3n with the flooded share within 500\,m of the nearest settlement in the same "
             r"event (the cell itself left out) & 6 \\"),
            (r"M1cG, M3oG & as M1c and M3n with recurrent surface water within 500\,m of the nearest settlement "
             r"(occurrence and seasonality, 1984--2021) & 6 \\"),
            r"M1cS, M3nS, M3nSp & as M1c, M3n and M3np with the strict reading (no springs or wells) & 6 \\",
            (r"T0$^-$ & T0 without river distances and without the water and wetland classes of land cover "
             r"(rung of the ladder) & 6 \\"),
            (r"M1cEc, M3oEc & as M1c and M3n with the share of other EMS events (more than 30 days apart) in "
             r"which the cell itself flooded, and whether there is such an event & 7 \\"),
            (r"M1cE, M3oE & the same share averaged over the cells within 500\,m of each settlement and assigned "
             r"through the nearest settlement & 7 \\"),
            r"M1cF, M3nF, M3nFp & as M1c, M3n and M3np with the floodplain reading (riverine and wetland classes) & 7 \\",
            (r"M1 + c(q), M1 + c(q) + q & reference and names models of a reassignment $q$, with the regional "
             r"composition recomputed from the reassigned names & 7 \\"),
            (r"M1cW, M3oW & as M1c and M3n with recurrent surface water of the cell itself (occurrence, seasonality "
             r"and maximum extent, 1984--2021, averaged over the cell) & 8 \\")]
    return _tbl("tab:posthoc-models", "Models of the post hoc analyses. All use the second stage of the primary "
                "comparison (Table~\\ref{tab:models}); the amendment is the one that planned the model.", "lp{11.4cm}l",
                "Code & Information & Amendment", rows)


def supp_tables_b7(v37: dict | None, v43: dict | None, v44: dict | None, v55: dict | None, v56: dict | None,
                   v52: dict | None) -> None:
    """補足の表（B-7）: 推定値の区分・B-2 の翻訳・Wikidata・言語別の照合・除外 AOI の面積・B-7 の結果。"""
    # 推定値の区分（B-7 の 11(k)）
    P, PH, CB, NP, CO = ("pre-specified", "post hoc, planned before computed", "post hoc, computed before planned",
                         "not planned", "correction")
    status = [
        ("H2 verdict, $R$, readings, placebo, capture; protocol's two-way bootstrap", P, "protocol", "3.2"),
        ("H1, H3, H4, H5; total and flexible increments; secondary splits; seven sensitivity analyses", P,
         "protocol", "3.3--3.4, C"),
        ("$p$ values of H2--H4 (computed after the critiques)", P, "6, add.~1", "3.2--3.3, A"),
        ("Rerun with the corrected reading of Norwegian names; corrected sensitivity analysis 6", CO, "4; 1", "3.2, C"),
        (("Positive control on all cells; alternative yardsticks; two-way bootstrap variant; distance bands from "
         "500\\,m; within-AOI odds ratios on all cells; fold assignments; events held out whole"), PH, "2",
         "3.4, 4, 6.1"),
        (("Class split of the positive control; overlap of top deciles; bound translations under other fold "
         "assignments; split-specific power ratios; interval coverage; class intervals; real minus permuted by "
         "class; operating characteristics against the real margin; interval combining fold assignments"), CB, "2",
         "2.6, 3.4, S"),
        (("Positive control at settlement points and by distance band; flood contrasts and increments below "
         "500\\,m; odds ratios by class and band and within 500\\,m of water-side names"), NP, "2", "5.1, 6.1"),
        ("GFM labels, designated replication, rescoring, odds ratios", PH, "3", "5.3"),
        ("Check of GFM masking (planned after the GFM result had been seen)", NP, "3, add.~2", "5.3"),
        (("Permutation distribution; coefficients; thirteen fold assignments; test inversion; terrain oracle; "
         "within-area metrics; label variants; element groups; settler societies excluded; etymology checks"), PH,
         "5", "3.4, 4, 5, 6.3"),
        (("Ceiling; plasmode with fixed base models; flood-information controls; GFM diagnostics and V1; strict "
         "reading; sensitivity runs; AOI weighting; rounded coordinates; ladder rung T0$^-$"), PH, "6", "3, 4, 5, 6"),
        (("Permutations with recomputed composition; plasmode with refitted models; activation split with the "
         "corrected code; recurrence control; floodplain reading; partly flooded cells; magnitude terciles; "
         "built-up and tree adjustment; GFM within-AOI discrimination; pairs bootstrap; descriptions (a)--(k)"), PH,
         "7", "3.4, 4, 5, A, B"),
        ("Joint tests of the near-settlement coefficients; digraph matching; recall by language", PH, "7, add.~1",
         "5.1, 6.3, B"),
        ("Joint test of the lexicon coefficients on the held-out predictions", PH, "7, add.~3", "3.4"),
        ("Observed increment within 500\\,m on the pooled scale (re-expression of an amendment-2 quantity)", NP,
         "--", "4.1"),
        (("Achieved equivalence bound; location groups needed; permutation test; second stage refitted in the "
         "bootstrap; operating characteristics at the margin; within-500\\,m increment on its own scale"), PH, "8",
         "3.2--3.3, 4.1, S"),
        (("No-effect reference and Monte Carlo intervals of the bounds; more plasmode replicates; recurrent surface "
         "water of each cell; cleaned lexicon"), PH, "8", "3.4, 4"),
        (("Estimates near settlements by class and variant (uncapped sample, misclassification, river size, the two "
         "label sets on the same cells); heights above main rivers; calibration and skill; audit of matched forms; "
         "reading by nearby river and lake names; second coder"), PH, "8", "5, 6, S"),
    ]
    rows = [f"{a} & {b} & {c} & {d} \\\\" for a, b, c, d in status]
    (HERE / "app_table_seven_status.tex").write_text(_tbl(
        "tab:seven-status", "Status of the reported estimates, by group. Pre-specified: fixed in the frozen protocol; "
        "post hoc, planned before computed: fixed in an amendment before it was computed; post hoc, computed before "
        "planned: first computed by the agents whose critiques prompted the amendment, which then fixed how it is "
        "computed and reported; not planned: added without an amendment fixing it in advance; correction: rerun after "
        "an implementation error. Section: where the estimates are reported (A--C: appendices; S: Supplement).",
        "p{9.6cm}p{3.2cm}ll", "Estimates & Status & Amendment & Section", rows), encoding="utf-8")
    # B-2 の合成ラベルの翻訳（B-7 の 11(a)）
    e = (v37 or {}).get("effect_size_translation") or {}
    if e:
        def orr(d: dict) -> str:
            if d.get("odds_ratio") is None:  # 模擬した効果の外（格子の上限を書く、編集者の指摘 10）
                dm = d.get("delta_max")
                return f"\\ensuremath{{>}}\\,{np.exp(dm):.0f}" if dm is not None else "beyond the grid"
            lo, hi = d.get("odds_ratio_mc_range") or [None, None]
            return orv(d["odds_ratio"]) + (f" ({orv(lo)}--{orv(hi)})" if lo is not None else "")
        rows = []
        for key, lab in (("margin", "Margin ($R=0.05$)"), ("achieved_upper90", ("Upper 90\\,\\% limit of the observed "
                                                                               "increment")),
                         ("two_way_upper90", "Same, two-way bootstrap")):
            rows.append(f"{lab} & {orr(e['homogeneous'][key])} & {orr(e['local_500m'][key])} \\\\")
        rows.append("\\middlehline")
        names = {"world_v2": "confirmatory run", "world_v2_fs1": "fold assignment 1",
                 "world_v2_fs2": "fold assignment 2", "world_v2_s6_buffer_b1": "location groups formed with a 10\\,km buffer",
                 "world_v2_actgrp": "events held out whole"}
        for run, d in e["by_run"].items():
            h, loc = d["homogeneous"], d["local_500m"]
            rows.append(f"Upper limit, {names.get(run, run)} & {orr(h)} & {orr(loc)} \\\\")
        (HERE / "app_table_seven_translations.tex").write_text(_tbl(
            "tab:seven-translations", "Odds ratios of flooding that the margin and the observed upper limits correspond "
            "to under the synthetic labels of the protocol's power analysis (amendment 2): an effect on all ground "
            "nearest a water-side name, and on the ground within 500\\,m of such a settlement. In parentheses the range "
            "over Monte Carlo resamples of the simulated curve; $>$: beyond the largest simulated effect (odds ratio "
            f"{np.exp(e['homogeneous']['margin']['delta_max']):.2f} for all ground, "
            f"{np.exp(e['local_500m']['margin']['delta_max']):.1f} within 500\\,m). These translations assume the label model of the power "
            "analysis; the simulations on the observed labels (Table~\\ref{tab:hazard}) give lower values.", "lll",
            "Target of the translation & All ground nearest the name & Within 500\\,m", rows), encoding="utf-8")
    # Wikidata P138（B-7 の 11(c)）
    if v43:
        def kr(d: dict) -> str:
            return f"{d['k']} of {d['n']} ({100 * d['rate']:.0f}\\,\\%)"
        rows = [f"Analysed settlements & {_n(v43['n_settlements'])} \\\\",
                f"\\quad with a Wikidata item & {_n(v43['n_with_item'])} \\\\",
                f"\\quad with a recorded namesake (P138) & {_n(v43['n_with_p138'])} \\\\",
                f"\\quad\\quad named after a water feature & {_n(v43['n_hydro'])} \\\\",
                f"\\quad\\quad named after a person, event or other commemoration & {_n(v43['n_commemorative'])} \\\\",
                "\\middlehline",
                f"Water-feature namesakes read as water-side by the lexicon & {kr(v43['recall_lexicon'])} \\\\",
                f"Same, by the language model & {kr(v43['recall_nread'])} \\\\",
                f"Lexicon water-side names with a water-feature namesake & {kr(v43['hydro_among_lexicon_wetside'])} \\\\",
                f"Language-model water-side names with a water-feature namesake & {kr(v43['hydro_among_nread_wetside'])} \\\\",
                (f"Lexicon water-side names with a commemorative namesake & "
                f"{kr(v43['commemorative_among_lexicon_wetside'])} \\\\")]
        (HERE / "app_table_seven_wikidata.tex").write_text(_tbl(
            "tab:seven-wikidata", "Settlements whose namesake is recorded on Wikidata (property P138), classed by the "
            "type of the namesake (amendment 5). The lexicon contains no river names by design, so it reads few "
            "settlements named after a river; most recorded namesakes are in France.", "lr", "Quantity & Count", rows),
            encoding="utf-8")
    # 言語別の精度・再現率（B-7 の 11(f) と 14）
    if v55 and v44:
        lp = v55["block_precision"]["by_language"]
        lr = (v55.get("language_recall") or {}).get("by_language") or {}
        st = (_json("v62_strict_by_language.json") or {}).get("by_language") or {}  # B-8 追記 7
        rows = []
        for lang in sorted(lp, key=lambda k: -lp[k]["n"]):
            p_ = lp[lang]["precision_lex"]
            r_ = lr.get(lang) or {}
            rec = (f"{_pct0(r_['recall'])} ({_pct0(r_['ci95'][0])}--{_pct0(r_['ci95'][1])})"
                   if r_.get("recall") is not None else "---")
            s_ = (st.get(lang) or {}).get("strict")
            sp = (f"{s_['k']} of {s_['n']} ({_pct0(s_['ci95'][0])}--{_pct0(s_['ci95'][1])})" if s_ and s_["n"]
                  else "---")
            rows.append(f"{lang} & {lp[lang]['n']} & {p_['k']} of {p_['n']} ({_pct0(p_['rate'])}; "
                        f"{_pct0(p_['ci95'][0])}--{_pct0(p_['ci95'][1])}) & {sp} & {rec} \\\\")
        (HERE / "app_table_seven_languages.tex").write_text(_tbl(
            "tab:seven-languages", "Precision and recall of the lexicon's water-side reading by language of the "
            "article, for languages with at least 30 coded names (amendment 5; recall by language, amendment 7, "
            "addendum 1). Recall is estimated as a domain with the stratum weights of the whole sample, with "
            "percentile intervals from 2,000 resamples within strata. The derivations were coded by the language model "
            "that compiled the lexicon. Strict: the strict reading (water without springs and wells, wetland, "
            "riverine) against strict derivations (amendment 8, addendum 7).", "lrlll",
            "Language & Coded names & Precision (95\\,\\% interval) & Strict precision & Recall "
            "(95\\,\\% interval)", rows), encoding="utf-8")
    # 除外した AOI の浸水面積（B-7 の 11(d)）
    if v55:
        ex = v55["excluded_area"]
        lab = {"included": "Analysed", "wet<20": "Fewer than 20 flooded cells", "eval<400": "Fewer than 400 evaluable "
               "cells", "wetfrac>0.8": "Flooded fraction above 0.8", "not-riverine": "Not riverine",
               "no_P_name": "No settlement"}
        rows = []
        for k in ("included", "wet<20", "eval<400", "wetfrac>0.8", "not-riverine", "no_P_name"):
            d = ex.get(k)
            if not isinstance(d, dict):
                continue
            rows.append(f"{lab[k]} & {d['n_aoi']} & {d['median_riverine_km2']:.1f} & "
                        f"{d['q10_q90_riverine_km2'][0]:.1f}--{d['q10_q90_riverine_km2'][1]:.1f} & "
                        f"{d['max_riverine_km2']:.1f} & {d['median_any_km2']:.1f} \\\\")
        (HERE / "app_table_seven_excluded.tex").write_text(_tbl(
            "tab:seven-excluded", "Flooded area of the labelled AOIs by inclusion criterion (amendment 7): riverine "
            "flood cover on the 25\\,m grid summed over evaluable ground (km$^2$; untyped flooding counted as riverine "
            "where the label rule does), and flood cover of any type.", "lrrrrr",
            "AOIs & Number & Median & 10th--90th percentile & Maximum & Median, any type", rows), encoding="utf-8")
    # 当て直しの plasmode（B-7 の 2）: 同じ δ・同じ回の基礎固定と当て直しの利得
    pr = (v56 or {}).get("plasmode_refit")
    if isinstance(pr, dict):
        b6g = pr.get("b6_fixed_gain_10reps") or [None] * len(pr["deltas"])
        rows = [(f"{np.exp(d):.2f} & {1e4 * pr['gain_fixed'][k]:.2f} & "
                 + (f"{1e4 * b6g[k]:.2f}" if b6g[k] is not None else "---")
                 + f" & {1e4 * pr['gain_refit'][k]:.2f} & {pr['ratio_refit_to_fixed'][k]:.2f} & "
                 f"{1e4 * pr['yardstick_injected_refit'][k]:.1f} & {100 * pr['label_flip_share'][k]:.2f} \\\\")
                for k, d in enumerate(pr["deltas"])]
    else:
        rows = [f"\\multicolumn{{7}}{{l}}{{{missing('plasmode refit')}}} \\\\"]
    (HERE / "app_table_seven_refit.tex").write_text(_tbl(
        "tab:seven-refit", "Effects injected into the observed labels for all ground nearest a water-side name, with the "
        "base models fixed and with every model refitted on the injected labels (amendment 7; three repetitions with "
        "the injections of amendment 6): mean gain in $\\Delta$BS over the uninjected run ($\\times10^{-4}$), the same "
        "with fixed base models over the ten repetitions of amendment 6, the ratio of the refitted to the fixed gain, "
        "the terrain increment (M0 to M1) of the refitted models in the injected labels ($\\times10^{-4}$) and the "
        "share of held-out labels changed by the injection (\\%). Odds ratios at the margin and bounds are not read "
        "from this simulation: it has no effect below an odds ratio of 1.65, and a straight line from zero across the "
        "convex gain curve would understate them (amendment 8, addendum 9).", "rrrrrrr",
        "Odds ratio & Fixed & Fixed, 10 rep. & Refitted & Ratio & Terrain increment & Labels changed", rows),
        encoding="utf-8")
    # B-7 の結果（近さのオッズ比と、その他の記述）
    if v55:
        rows = []
        pb = (v55.get("pairs_bootstrap") or {})
        for cls_key, block in (("wetside", pb.get("wetside") or {}), ("classes", pb.get("classes") or {})):
            for term, d in block.items():
                if term.startswith("_") or not isinstance(d, dict):
                    continue
                c, band = term.split(":")
                bl = {"le500m": "$\\leq$500\\,m", "500m_1km": "0.5--1\\,km", "gt1km": "$>$1\\,km"}[band]
                rows.append(f"{g.CLASS_LABELS.get(c, c) if c != 'wetside' else 'water-side'}, {bl} & {orv(d['odds_ratio'])} "
                            f"& {orv(d['sandwich_ci95'][0])}--{orv(d['sandwich_ci95'][1])} & "
                            f"{orv(d['pairs_ci95'][0])}--{orv(d['pairs_ci95'][1])} \\\\")
        if rows:
            (HERE / "app_table_seven_pairs.tex").write_text(_tbl(
                "tab:seven-pairs", "Odds ratios of flooding near settlements by class and distance band on the EMS "
                "labels, given the terrain model's prediction (AOI fixed effects), with cluster-robust and pairs "
                "cluster bootstrap 95\\,\\% intervals over location groups (2,000 replicates; amendment 7). The two "
                "blocks are separate models (water-side names; the six classes).", "lrll",
                "Names, distance & Odds ratio & Cluster-robust & Pairs bootstrap", rows), encoding="utf-8")


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
    v55, v56 = _json("v55_b7_local.json"), _json("v56_b7_cluster.json")
    m |= b7_macros(v55, v56, post["v37"], v52)
    v57, v57b, v59 = _json("v57_b8_local.json"), _json("v57b_second_coder.json"), _json("v59_b8_aggregate.json")
    m |= b8_macros(v57, v57b, v43, v41) | c8_macros(v59)
    m |= b9_macros(_json("v63_b9_local.json"), _json("v42_blind_review_analyses.json"))
    m |= lexicon_form_macros()
    m |= c9_macros(_json("v64_b9_aggregate.json"))
    b1 = (v59 or {}).get("B1") if isinstance((v59 or {}).get("B1"), dict) else None
    if b1 and b1.get("margin_or"):  # B-8 の B1: 50 回の曲線を主に（余白の換算も）
        m["PlasMarginRange"] = plas_margin_range(post["v37"], v52, v56, v59)[0]
        for f in ("A", "B"):
            m[f"PlasMargin{f}"] = _gt(b1["margin_or"][f]["odds_ratio_text"])
    # 割り当てで残る情報の割合（集落を通じた R / セルに直接の R）: 繰り返す浸水（B-7 の 4）と水域（B-6・B-8 の C4）
    rc = (v56 or {}).get("recur") if isinstance((v56 or {}).get("recur"), dict) else {}
    sub = (rc or {}).get("subset_aois_with_other_events") or {}
    g6 = ((v52 or {}).get("gsw_oracle") or {}).get("oracle_gsw") or {}
    g8 = ((v59 or {}).get("C4") or {}).get("gsw_cell") if isinstance((v59 or {}).get("C4"), dict) else None
    m["RecurKeptPct"] = (pct1(rc["recur_settlement"]["R"] / rc["recur_cell"]["R"]) if rc.get("recur_cell")
                         else missing("recurrence"))
    m["RecurSubKeptPct"] = (pct1(sub["recur_settlement"]["R"] / sub["recur_cell"]["R"]) if sub.get("recur_cell")
                            else missing("recurrence subset"))
    m["GswKeptPct"] = pct1(g6["R"] / g8["R"]) if g6.get("R") is not None and g8 else missing("GSW kept")
    m["NAoiAmericas"] = (str(int(m["NAoiNorthAmerica"]) + int(m["NAoiSouthAmerica"]))
                         if "NAoiNorthAmerica" in m and "NAoiSouthAmerica" in m else missing("Americas"))
    pw = g._json(DATA / "runs" / "power_v3_null" / "power_summary.json")
    m["SimHalfWidth"] = f"{pw['median_ci90_halfwidth']:.4f}" if pw else missing("power null")
    met = g._json(DATA / "runs" / "world_v2" / "metrics_primary.json") or {}
    br = met.get("brier") or {}
    m["MarginPctBS"] = (pct1(SESOI * (br["M0"] - br["M1"]) / br["M1c"]) if {"M0", "M1", "M1c"} <= br.keys()
                        else missing("brier for margin"))
    checks = posthoc_checks(v41, v42, v47, v48, v52, v55, v56, v59)
    have = [c for *_, c in checks if isinstance(c, dict) and c.get("decision")]  # 図 4b の事後の確認の判定の内訳
    # 判定は文字列のものと {"label": ...} の dict のもの（GFM の採点し直し）がある（初見の査読 第 5 回: dict を
    # 「equivalent-null 以外」と数えていた）
    labels = [_decision_label(c) for c in have]
    c2c = ((v59 or {}).get("C2") or {}).get("corrected") if isinstance((v59 or {}).get("C2"), dict) else None
    n_rows = (sum(1 for *_, c in checks if isinstance(c, dict) and c.get("R") is not None)
              + int(bool(isinstance(c2c, dict) and c2c.get("refit"))))  # 図 4b の行（先頭の当て直しの区間を含む）
    m |= {"NChecksWord": _word(len(have)), "NChecksRows": _word(n_rows),
          "NChecksOtherWord": _word(sum(lab != "equivalent-null" for lab in labels))}
    if primary:
        fig_main(primary, extra, post["v36"], v41, v42, v47, v48, v52, checks, v56, v59)
        if primary.get("ladder"):
            fig_ladder_b6(primary["ladder"], primary["primary"]["delta_bs_M0_M1"], v52)
    if post["v37"]:
        fig_names(post["v37"], post["v39"], v41, v42, v50, v52, v56)
    fig_plasmode(v52, v56, v59, _json("v64_b9_aggregate.json"))
    fig_near(v55, v57, _json("v61_leaveout.json"))
    t = pd.read_parquet(DATA / "verification" / "v41_table.parquet",
                        columns=["aoi", "y", "ix", "iy", "p_M1c", "p_M3n", "P_d1", "lex_wetside_n1"])
    gs = g._json(DATA / "global" / "gfm_summary.json")
    ex = pick_example(t, gs) if gs else None
    ex2 = pick_example_near(t, gs) if gs else None
    m |= (fig_example([ex] + ([ex2] if ex2 and ex2 != ex else []), t) if ex else
          {"ExampleAoiName": missing("example AOI")})
    m |= timeline(m) | gfm_items()
    write_new_tables(v41, v42, v44, v45, v47)
    table_hazard()
    (HERE / "app_table_gfm_compact.tex").write_text(table_gfm_compact(v41, v47, v50), encoding="utf-8")
    b6_supp_tables(v48, v50, v51, v52, v53, v54)
    (HERE / "app_table_posthoc_models.tex").write_text(table_posthoc_models(), encoding="utf-8")
    (HERE / "app_table_checks.tex").write_text(table_checks(checks), encoding="utf-8")
    supp_tables_b7(post["v37"], v43, v44, v55, v56, v52)
    v60, v61 = _json("v60_where_near.json"), _json("v61_leaveout.json")
    supp_tables_b8(v57, v57b, v59, v41, v60, v61)
    supp_tables_b9(_json("v63_b9_local.json"), _json("v64_b9_aggregate.json"), v57)
    m |= _orp("NoPkOR", (v61 or {}).get("without_pakistan")) | _orp("NoSAsiaOR", (v61 or {}).get("without_pakistan_bangladesh"))
    if v60:  # 水系側の名の近くの浸水セルの散らばり（編集者の求め、記述）
        m |= {"WhereFlooded": str(v60["flooded_near_total"]),
              "WhereNCountries": str(v60["n_countries_with_flooded_near"]),
              "WhereNActivations": str(v60["n_activations_with_flooded_near"]),
              "WhereTopFivePct": pct1(v60["share_top5_countries"]),
              "WhereTopActPct": pct1(v60["largest_activation_share"]),
              "WhereTopCountries": ", ".join(r["cc"] for r in v60["by_country"][:5])}
    else:
        m |= {k: missing("v60") for k in ("WhereFlooded", "WhereNCountries", "WhereNActivations", "WhereTopFivePct",
                                          "WhereTopActPct", "WhereTopCountries")}
    tables = copernicus_tables()
    m |= supplement()
    write_numbers(m)
    print(f"macros: {len(m)}; tables: {len(tables)}; missing: {len(g.MISSING)}")
    for w in sorted(set(g.MISSING)):
        print("  missing:", w)


if __name__ == "__main__":
    main()
