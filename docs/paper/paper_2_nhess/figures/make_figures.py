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
        g.DBS_COLUMN = True  # 表 2・C1 に Brier の単位の変化の列（初見の査読 第 7 回、査読者 2 の M2）
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
    ax.text(17.5, 28.5, "baseline gain\n(denominator of R)", fontsize=7.4, ha="center", va="top",
            color=grey)
    ax.text(74.8, 28.5, "names increment\nR = names increment / baseline gain", fontsize=7.4, ha="center", va="top",
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
        # B-12 の A12: R が未定義の分割には点を打たない（文字だけ）
        spec.append(("Time held out (floods from 2024)", (None, None, None) if und else blk(met["primary"]),
                     "undefined" if und else "main"))
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
            ax.text(0.0, y, "R undefined: the baseline did not improve on position", ha="center", va="center",
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
    # 初見の査読 第 6 回（編集者の 5 節）: 本文の図 4 は事前に固定した比較だけ、事後の確認と範囲は補足の図 S1 に
    # B-12 の A12: 点線の下に、活性化をまとめた探索の分割（凍結したコードと直したコード）
    af = (v56 or {}).get("actfix") if isinstance((v56 or {}).get("actfix"), dict) else {}
    for mp, lab in (((g._json(DATA / "runs" / "world_v2_actgrp" / "metrics_primary.json") or {}).get("primary") or {},
                     "Activations kept whole (exploratory), frozen code"),
                    (af.get("primary") or {}, "   Same, corrected code")):
        if mp.get("R") is not None:
            rows_a.append((lab, mp["R"], mp.get("ci90_R"), mp.get("ci95_R"), "post"))
    na = len(rows_a)
    fig, ax0 = plt.subplots(figsize=(6.6, 0.122 * na + 0.95))
    _draw_rows(ax0, rows_a, style)
    for y, (nm, *_) in zip(np.arange(na)[::-1], rows_a, strict=True):
        if nm.startswith("Countries held out"):
            ax0.axhline(y + 0.5, color="0.8", lw=0.6)
        if nm.startswith("Activations kept whole"):
            ax0.axhline(y + 0.5, color="0.6", lw=0.6, ls=":")
    ax0.set_xlabel("R (improvement from names / baseline gain)", fontsize=8)
    ax0.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "fig_main.pdf")
    plt.close(fig)
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 0.122 * (nb + nc) + 1.1), sharex=True,
                            gridspec_kw={"height_ratios": [nb + 0.6, nc + 0.6]})
    _draw_rows(axs[0], rows_b, style)
    ax = axs[1]
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
    for a, (lab, title) in zip(axs, (("a", "Exploratory checks"), ("b", "Ranges, not intervals")), strict=True):
        a.set_title(f"({lab}) {title}", fontsize=7.5, loc="left", fontweight="bold")
    axs[1].set_xlabel("R (improvement from names / baseline gain)", fontsize=8)
    axs[1].tick_params(axis="x", labelsize=7)
    g.SOURCE["Fig4"] = [{"panel": p, "comparison": n.strip(), "R": r, "ci90_low": c9[0] if c9 else None,
                         "ci90_high": c9[1] if c9 else None, "ci95_low": c5[0] if c5 else None,
                         "ci95_high": c5[1] if c5 else None, "group": st}
                        for p, rows in (("a", rows_a), ("S1a", rows_b)) for n, r, c9, c5, st in rows]
    fig.tight_layout(h_pad=0.6)
    fig.savefig(HERE / "fig_checks.pdf")
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
    ax.set_ylabel("names within areas\n(fraction of baseline gain)", fontsize=8)
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
              ("flooded, net of reference model (M1c)",
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
def fig_near_more(v55: dict | None, v57: dict | None, v61: dict | None = None) -> None:
    """補足の図 S3: 集落の近くの他の推定（区分・氾濫原の読み・標本と読みの変種・センサー・川の大きさ・誤読の補正）。

    初見の査読 第 7 回: 群ごとに区間の種類を一つにし、群の見出しに書く。
    """
    v57, v61 = v57 or {}, v61 or {}
    b6, b7, b12 = (v57.get(k) or {} for k in ("B6", "B7", "B12"))
    blue, grey, black = "#1f6fb4", "0.45", "0.15"
    groups: list[tuple[str, list]] = []
    ca = b6.get("class_adjusted") or {}
    groups.append(("Lexicon classes, one model (Bonferroni-adjusted 95 % intervals)",
                   [(lab, ca[c]["odds_ratio"], ca[c]["bonferroni_ci95"], blue, True) for c, lab in CLASS_LABELS
                    if c in ca]))
    v66b2 = ((_json("v66_b10_local.json") or {}).get("B2")) or {}
    rows = []
    for key, lab in (("main", "Floodplain reading (defined after the class estimates)"),
                     ("reread_without_forms", "  without bank, au, dike and polder forms")):
        d = v66b2.get(key) or {}
        if d.get("pairs_ci95"):
            rows.append((lab, d["odds_ratio"], d["pairs_ci95"], blue, key == "main"))
    groups.append(("Floodplain reading (pairs bootstrap)", rows))
    rows = []
    for key, lab in (("main", "All AOIs"), ("typed_aois", "AOIs with a recorded flood type"),
                     ("without_no_hydro", "Without AOIs lacking hydrography"),
                     ("without_oceania_americas", "Without Oceania and the Americas"),
                     ("cleaned_reading", "Cleaned reading")):
        d = b6.get(key)
        if isinstance(d, dict) and d.get("pairs_ci95"):
            rows.append((lab, d["odds_ratio"], d["pairs_ci95"], black, key == "main"))
    d = v61.get("without_pakistan_bangladesh")
    if isinstance(d, dict) and d.get("pairs_ci95"):
        rows.append(("Without Pakistan and Bangladesh", d["odds_ratio"], d["pairs_ci95"], black, False))
    v64 = _json("v64_b9_aggregate.json") or {}  # B-9 の C5・C7（pairs bootstrap）
    for key, lab in (("C5_near", "One definition of permanent water"), ("C7_near", "Bare-earth elevation model")):
        d = v64.get(key)
        if isinstance(d, dict) and d.get("pairs_ci95"):
            rows.append((lab, d["odds_ratio"], d["pairs_ci95"], black, False))
    groups.append(("Capped sample: variants of sample and reading (pairs bootstrap)", rows))
    sens = ((_json("v66_b10_local.json") or {}).get("B7") or {}).get("sensor") or {}
    rows = []
    for k, lab in (("s1", "Sentinel-1"), ("sar", "Other radar"), ("optical", "Optical"), ("none", "Not recorded")):
        o = (sens.get(k) or {}).get("or_le500m") or {}
        if o.get("pairs_ci95") and o.get("pairs_log_or_sd", 0.0) <= PAIRS_UNSTABLE_SD:
            rows.append((lab, o["odds_ratio"], o["pairs_ci95"], black, True))
    groups.append(("Capped sample, by the images behind the EMS maps (pairs bootstrap)", rows))
    rs = b7.get("odds_by_river_size") or {}
    rows = [(lab, rs[k]["odds_ratio"], rs[k]["ci95"], grey, True)
            for k, lab in (("lt100", "Nearest river < 100 km$^2$"), ("100_1000", "Nearest river 100–1000 km$^2$"),
                           ("ge1000", "Nearest river ≥ 1000 km$^2$")) if k in rs]
    if b12.get("near_500m"):
        rows.append(("Name contains a nearby river or lake name", b12["near_500m"]["odds_ratio"],
                     b12["near_500m"]["ci95"], grey, False))
    groups.append(("Capped sample, by river size; reading by river names (cluster-robust)", rows))
    mc = b6.get("misclassification_corrected")
    if isinstance(mc, dict):
        groups.append(("Corrected for misreadings (2.5th–97.5th percentiles of a bias analysis)",
                       [("Truly water-related names, illustrative rates", mc["median"], mc["ci95"], black, False)]))
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
    ax.set_xlabel("odds ratio of flooding within 500 m, given the baseline model's prediction", fontsize=7.5)
    ax.tick_params(axis="x", labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    g.SOURCE["FigNearMore"] = src
    fig.savefig(HERE / "fig_near_more.pdf")
    plt.close(fig)


def fig_plasmode(v52: dict | None, v56: dict | None, v59: dict | None = None, v64: dict | None = None,
                 v67: dict | None = None) -> None:
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
    c4 = (v67 or {}).get("C4_ac") if isinstance((v67 or {}).get("C4_ac"), dict) else None
    if c4:  # B-10 の C4: 反転を集落の地面の中で相関させた注入（形 AC、基礎は固定）
        ors = np.exp(np.concatenate([[0.0], np.asarray(c4["deltas"], dtype=float)]))
        yv = 1e4 * np.concatenate([[0.0], np.asarray(c4["gain_mean"], dtype=float)])
        se = 1e4 * np.concatenate([[0.0], np.asarray(c4["gain_mc_se"], dtype=float)])
        # 初見の査読 第 7 回（査読者 2 の細かい指摘 12）: 余白の水平線の破線と区別するため、点 2 つの鎖線にする
        ax.plot(ors, yv, ls=(0, (4, 1.2, 1, 1.2, 1, 1.2)), color=blue, lw=1.1, marker="h", ms=3.0, mfc="white",
                mew=0.7,
                label=f"all ground, flips correlated within each settlement's ground ({c4['reps']} replicates)")
        ax.fill_between(ors, yv - se, yv + se, color=blue, alpha=0.10, lw=0)
        src["AC"] = {"odds_ratio": ors.tolist(), "gain": (yv / 1e4).tolist(), "rho": c4.get("rho")}
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
        for j, (ax, title) in enumerate(zip(row, ("EMS map", "GFM labels", "HAND", "reference model (M1c)", "change from names"),
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
            s = _add_flex_intervals(s)
            s = s.replace("$R$ is the change in held-out Brier score as a fraction",
                          "$R$ is the reduction in held-out Brier score as a fraction", 1)
            s = s.replace("}\\label{tab:si-comparisons}", " $\\Delta$BS: the reduction in held-out Brier score from the "
                          "comparison (negative when the predictions got worse), in units of $10^{-5}$, with the 95\\,\\% "
                          "cluster-bootstrap interval; the margin is "
                          "\\MarginDBSPrimary.}\\label{tab:si-comparisons}", 1)
        if f.name == "si_table_sensitivity.tex":  # 初見の査読 第 7 回（査読者 1 の 14）: プロトコル §7 の番号を行に
            a_, b_ = s.index("\\middlehline\n") + len("\\middlehline\n"), s.index("\\bottomhline")
            body = s[a_:b_].rstrip("\n").split("\n")
            nums = ("1:", "2a:", "2b:", "3:", "4:", "5:", "6:", None, "7:")
            if len(body) == len(nums):
                body = [ln if n is None else f"{n} {ln}" for ln, n in zip(body, nums, strict=True)]
                body[7] = (body[7].replace("\\quad as first run,", "\\quad 6, first run:")
                           .replace("folds by 5\\,km places", "folds by the 5\\,km location groups"))
                s = s[:a_] + "\n".join(body) + "\n" + s[b_:]
            s = s.replace("}\\label{tab:si-sensitivity}", " Numbers: the sensitivity analyses of protocol section 7; the "
                          "second has two runs (2a, 2b), and the sixth a first run with folds assigned by the wrong "
                          "location groups (amendment 1).}\\label{tab:si-sensitivity}", 1)
        if f.name == "si_table_splits.tex":
            s = s.replace("}\\label{tab:si-splits}", " Baseline $\\Delta$BS: the baseline gain of the split (the "
                          "denominator of $R$); names $\\Delta$BS: the reduction in held-out Brier score from names "
                          "(negative when the predictions got worse), in units of $10^{-5}$, with the 95\\,\\% "
                          "cluster-bootstrap interval.}\\label{tab:si-splits}", 1)
        out = HERE / f.name.replace("si_table_", "app_table_")
        out.write_text(s, encoding="utf-8")
        names.append(out.name)
    return names


def _d5(x: float) -> str:
    """表 2 の ΔBS（10^-5 の単位、小数 1 桁、負は ${-}$）。"""
    v = x * 1e5
    return f"{'${-}$' if v < 0 else ''}{abs(v):.1f}"


def _add_h3_rows(s: str) -> str:
    """表 2 に H3 の事前に決めた量（T0 と T4 の差と片側の p、T0 で取り戻した地形の情報の割合）を足す（編集者の指摘 8）。
    B-12 の A12: T4 の行は H2 と柔軟な当て直しの行の重複なので消し、差の区間を正しい列に置く。"""
    met = g._json(DATA / "runs" / "world_v2" / "metrics_primary.json") or {}
    lad = met.get("ladder")
    key = "H3 ladder, within areas, T4 = M1 (M1c$\\to$M3n)"
    if not lad or key not in s:
        return s
    d01 = met["primary"]["delta_bs_M0_M1"]
    j = s.index("\\\\", s.index(key)) + 2
    rows = (f"\n\\quad T0 minus T4 (T4 is H2), one-sided $p={lad['p_one_sided_delta0_gt_delta4']:.2f}$ & "
            f"{g._r(lad['delta0_minus_delta4'] / d01)} & --- & "
            f"{g._ci([x / d01 for x in lad['ci95_delta0_minus_delta4']])} & --- & --- & --- & --- \\\\"
            f"\n\\quad share of terrain information recovered at T0$^{{\\ddagger}}$ & "
            f"{g._r(lad['terrain_info_recovered_at_T0'])} & --- & --- & --- & {g._ci(lad['ci95_recovered_at_T0'])} & --- & "
            "--- \\\\")
    s = s[:j] + rows + s[j:]
    for k in (key, "H3 ladder, flexible, T4 = M1 (M1$\\to$M3)"):  # 重複の行を消す
        a_ = s.index(k)
        s = s[:a_] + s[s.index("\\\\", a_) + 3:]
    return s.replace("}\\label{tab:si-comparisons}", " $^{\\ddagger}$The names increment at T0 as a share of what "
                     "elevation and HAND add (T0 to T4), with a ratio-bootstrap interval. The top rung of the ladder "
                     "(T4) is the H2 comparison within areas and the flexible refit M1$\\to$M3."
                     "}\\label{tab:si-comparisons}", 1)


def _add_flex_intervals(s: str) -> str:
    """表 2 の柔軟な「辞書だけ」「学習した得点だけ」に、ΔBS の区間と固定した余白の 95 % 区間（B-12 の B8）を足す。"""
    b8 = (_json("v69_b12_local.json") or {}).get("B8") or {}
    for arm, lab in (("M3x", "\\quad lexicon only, M1$\\to$M3x"), ("M3l", "\\quad learned score only, M1$\\to$M3l")):
        x = b8.get(arm) or {}
        if not x or lab not in s:
            continue
        a_ = s.index(lab)
        e_ = s.index("\\\\", a_)
        cells = [c.strip() for c in s[a_:e_].split("&")]
        if len(cells) != 8:
            continue
        lo, hi = x["ci95_delta_bs"]
        cells[2] = f"{_d5(x['delta_bs'])} [{_d5(lo)}, {_d5(hi)}]"
        cells[3] = g._ci(x["ci95_R_fixed"])
        s = s[:a_] + " & ".join(cells) + " " + s[e_:]
    return s


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
         "BNine": "B-9_blind_review5_analyses.md", "BTen": "B-10_blind_review6_analyses.md",
         "BEleven": "B-11_blind_review7_analyses.md", "BTwelve": "B-12_blind_review8_analyses.md"}
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
    "38e8dc9": "Addendum 1 to amendment 10 (implementation record and a dilution table, before any result of part C)",
    "e682643": "Addendum 2 to amendment 10 (placebo comparison of the GFM within-area evaluation, after its result had "
               "been seen)",
    "13de450": "Addendum 1 to amendment 12 (census flags near settlements corrected, an implementation error; "
               "Benjamini--Yekutieli; 50 names added to the audit sample; before the quantities were computed)",
    "dd007b0": "Addendum 2 to amendment 12 (label flags of the specification curve corrected, an implementation error; "
               "numbering of the archived tables; before recomputing)",
    "1881566": "Addendum 3 to amendment 12 (the specification curve takes the estimator's versions flagged by keys; "
               "before recomputing)",
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
    blind6 = _first_commit_time("docs/research/reviews/2026-10-29_paper2_blind6_editor.md")
    blind7 = _first_commit_time("docs/research/reviews/2026-10-29_paper2_blind7_editor.md")
    # 初見の査読 第 7 回（編集者の 4 節）: 記録を見られた批評の回も年表に
    access = [("", "docs/research/reviews/2026-10-26_paper2_opus_referee.md")] + [
        (f" round {k}", f"docs/research/reviews/2026-10-26_paper2_opus_referee_round{k}.md") for k in (2, 3, 4, 5)]
    pub = g._json(VER / "public_repo.json") or {}  # GitHub の API の作成時刻（v66 と同じ日に取った）
    out["PublicRepoCreatedUTC"] = (pub.get("created_at") or "")[:16] + "Z" if pub.get("created_at") \
        else missing("public repository")
    out["ConfirmatoryUTC"] = (conf.get("created_utc") or "")[:16] + "Z" if conf.get("created_utc") \
        else missing("confirmatory record")
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
            ("Amendment 9 (exploratory analyses after the fifth critiques)", out["AmendTimeBNine"]),
            ("Public code repository created (GitHub)", out["PublicRepoCreatedUTC"]),
            ("Sixth LLM-agent critiques of the manuscript recorded", blind6 or ""),
            ("Amendment 10 (corrections and exploratory analyses after the sixth critiques)", out["AmendTimeBTen"]),
            ("Seventh LLM-agent critiques of the manuscript recorded", blind7 or ""),
            ("Amendment 11 (refocusing, corrections and exploratory analyses after the seventh critiques)",
             out["AmendTimeBEleven"])]
    rows += [(f"LLM-agent critique with access to the records{lab or ' round 1'} recorded", _first_commit_time(f) or "")
             for lab, f in access]
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
    cap = ("Post hoc analyses of amendment B-5 on the primary split. Superseded: the upper bounds of this table "
           "understated the uncertainty and are replaced by the simulations on the observed labels "
           "(Table~\\STabNineBounds). The largest odds ratios are upper bounds by test "
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
           "the ground that the name refers to, which lies within the floodplain. The expected directions were written "
           "after the odds ratios of the wetland, riverine and crossing classes near settlements had been computed "
           "(amendment 2, addendum 3), and the wetland group is the wetland class. Contrary to the expectation, "
           "wetland names did not lie lower than other names. The island group does not contain Bengali \\emph{char}, "
           "a river-island term that matched \\AuditCharN{} settlements, and the springs-and-wells group contains English "
           "\\emph{-bourne} and \\emph{-borne}, which continue Old English \\emph{burna}, ``stream''. "
           "$^{a}$Includes German \\emph{Au} and \\emph{Aue}, which mean water meadow rather than island.")
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
           "the regional composition (AOI means). The reading of Norwegian names is the corrected one (B-4). The "
           "indicator that the lexicon covers a language of the country is one for every analysed cell with the real "
           "names, so its standardised coefficient is zero and the second stage has effectively fourteen name "
           "features.")
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
    # B-12 の A12: どの複製も余白の外に出ないときは p < 1/B（B は bootstrap の回数）と書く
    m["HTwoTostP"] = (("< 1/B = " + f"{tost['resolution']:.5f}".rstrip("0")) if tost["p_tost"] == 0
                      else f"= {tost['p_tost']:.4f}")
    if "H4_tost" in pv:
        t4 = pv["H4_tost"]
        m["HFourTostP"] = (("< 1/B = " + f"{t4['resolution']:.5f}".rstrip("0")) if t4["p_tost"] == 0
                           else f"= {t4['p_tost']:.4f}")
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
                                  "class (a name can match at more than one position, so the position rows overlap and "
                                  "their counts sum to more than the names); recall by class is unweighted and therefore "
                                  "not comparable with the weighted "
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
        rows.append("Null simulations: between-replicate SD of $R$ (robust SD; outlying replicates); median bootstrap "
                    "SE of $R$ & \\NullSdBetween{} (\\NullRobustSd; \\NullNOut); \\NullMedianSe \\\\")
        rows.append("Null simulations: verdicts & " + ", ".join(f"{k} {v}" for k, v in dc.items()) + " \\\\")
        pv = v54["p_values"]["world_v2"]
        rows.append(f"H2: TOST $p$ & $<1/B$ ({pv['H2_tost']['resolution']:.5f}) \\\\")  # B-12 の A12
        rows.append(f"H3: one-sided $p$ ($\\Delta_0>\\Delta_4$); Spearman $\\rho$ and exact one-sided $p$ & "
                    f"{pv['H3']['p_delta0_gt_delta4']:.2f}; {pv['H3']['spearman']:.1f}, {pv['H3']['p_spearman_exact']:.2f} \\\\")
        if "H4_tost" in pv:
            rows.append(f"H4: TOST $p$; one-sided $p$ for an improvement & $<1/B$ ({pv['H4_tost']['resolution']:.5f}); "
                        f"{pv['H4_improvement_p']:.2f} \\\\")
        rows.append(f"Holm-adjusted $p$, H3 / H4 & {pv['holm']['H3']:.2f} / {pv['holm']['H4']:.2f} \\\\")
        rc = v54["residual_correlogram"]
        rows.append("Correlation of held-out M1 residuals within AOIs by distance (km) & " + "; ".join(
            f"{b}: {c:.2f}" for b, c in zip(rc["bins_km"], rc["corr"], strict=True)) + " \\\\")
        rp = v54["reference_points"]
        hb = _hand_b15()
        b15 = (_json("v69_b12_local.json") or {}).get("B15") or {}
        if hb and b15.get("M1") and b15.get("M3n"):  # B-12 の B15: HAND だけの地図を一つの定義に
            rows.append(f"CSI pooled over all cells of maps flagging the observed flooded share per AOI: mean HAND "
                        f"(Table~\\STabTenHand) / M1 / M3n & {hb['csi_pooled']:.2f} / {b15['M1']['csi_pooled']:.2f} / "
                        f"{b15['M3n']['csi_pooled']:.2f} \\\\")
        rows.append(f"Brier score: constant (training rate) / M0 / M1 & {rp['brier']['const_train_rate']:.4f} / "
                    f"{rp['brier']['M0']:.4f} / {rp['brier']['M1']:.4f} \\\\")
        out["b6_stats"] = _tbl("tab:b6-stats", "Interval coverage in the pre-specified null simulations (amendment 2), "
                               "the $p$ values required by the protocol (amendment 6, addendum 1), the spatial correlation of "
                               "held-out residuals and reference points for the baseline model (primary split). The CSI "
                               "maps use the observed flooded count of each AOI and are therefore optimistic. TOST $p$: "
                               "no bootstrap replicate fell outside the margin, so $p$ is below $1/B$ for $B$ replicates "
                               "(the value in brackets).",
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
        ("amendments", "The content of each amendment."),
        ("seven_status", "The planning status of each group of reported estimates."))),
    ("Pre-specified analyses", (
        ("splits_secondary", "The secondary comparisons of the protocol under each split."),
        ("registered", "Quantities planned in amendment 2 and not shown elsewhere, including the protocol's two-way "
                       "bootstrap."),
        ("yardsticks", "The within-area comparison judged against five measures of the baseline gain and two bootstraps."),
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
        ("ten_hand", "The skill of flood maps, with one definition of the HAND-only map, per AOI and pooled."),
        ("ten_robust", "The primary comparison by sensor, product generation and without Africa, and with the "
                       "alternative second stages and references of amendment 10."),
        ("ten_splits", "Brier scores of every split."),
        ("eleven_checks", "The activation split in two code versions, further intervals, the GFM domain and the "
                          "denominator of the labels."),
        ("eleven_census", "A census of the exploratory interval estimates."),
        ("twelve_checks", "Operating characteristics of the two decision rules, the reference model, constant forecasts, "
                          "outlying null simulations, inclusion, flood typing, GFM images, built-up land near "
                          "settlements, where the flooded cells near water-side names lie and a penalty on the names "
                          "alone."),
        ("composition", "What the regional composition of names adds."),
        ("ten_between", "Whether areas with more water-side names flood more, within language regions."),
        ("b6_runs", "Exploratory runs and controls of amendment 6."),
        ("b5", "Test inversion and the other analyses of amendment 5."))),
    ("What assignment through settlements can carry", (
        ("b6_plasmode", "Gains from effects injected into the observed labels, base models fixed."),
        ("seven_refit", "The same gains with every model refitted on the injected labels."),
        ("eight_bounds", "Bounds, detectable effects and margin translations with different no-effect references, with "
                         "Monte Carlo intervals."),
        ("nine_bounds", "Bounds counted three ways, in both directions and for effects that differ between settlements."),
        ("nine_power", "Simulated power of the decision rule with every model refitted."),
        ("ten_injected", "The estimator of Sect. 5 on labels with a known injected odds ratio."),
        ("ten_dilution", "How an effect on the ground a name refers to is diluted on all ground nearest the name."),
        ("eleven_scale", "The margin and the bounds on the scale of the estimator of Sect. 5."),
        ("seven_translations", "Translations of the margin and bounds under the synthetic labels of the protocol."),
        ("b6_gsw", "Recurrent surface water around settlements, by name class (no flood labels)."))),
    ("Flooding near settlements", (
        ("positive_control", "Where settlements with names of each class stand relative to the drainage, by distance "
                             "(no flood labels)."),
        ("eight_heights", "Heights of settlements above the main rivers, by name class (no flood labels)."),
        ("odds", "Within-AOI odds ratios of flooding by class and distance."),
        ("seven_pairs", "Odds ratios near settlements by class and distance band, with pairs-bootstrap intervals."),
        ("class", "The names increment by class of the nearest name."),
        ("band_contrast", "Flood contrasts and increments by distance from the settlement."),
        ("groups", "Element groups of the lexicon."),
        ("ten_islands", "The group of island and water-meadow names at the settlement sites and without the meadow "
                        "terms."),
        ("eight_near", "The pooled odds ratio near settlements under every variant, by river size and by class, and "
                       "the two label sets on the same cells."),
        ("eight_where", "Where the flooded cells near water-side names lie, by country and by activation."),
        ("nine_activations", "The pooled odds ratio near settlements with each activation left out."),
        ("nine_classq", "Adjusted $p$ and $q$ values of the class coefficients."),
        ("ten_misclass", "The odds ratio near water-side names corrected over a grid of misreading rates."),
        ("eleven_near", "The odds ratio near water-side names on every evaluable cell, with both label sets."),
        ("ten_floodplain", "The floodplain reading by country, without one activation and without doubtful forms."),
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
        ("ten_nread", "The classes of the language-model reading."),
        ("eleven_gazetteer", "The settlements by GeoNames feature code, and the pilot simulations."),
        ("seven_wikidata", "Settlements whose namesake is recorded on Wikidata."),
        ("seven_languages", "Precision and recall by language."),
        ("ten_valid", "Precision and recall by country and by form, for both coders."),
        ("ten_other_lang", "Cells whose names are read in a language other than their own."),
        ("eleven_headings", "The Wikipedia headings taken as explanations of the name, by language."))),
    ("Further quantities", (
        ("ten_text", "Quantities stated in the paper and not given in another table."),)),
)
SUPP_TABLES = tuple(name for _, items in SUPP_GROUPS for name, _ in items)
PAPER_TABLE_WORDS = {"tab:models": "Table~1 of the paper",
                     "tab:si-comparisons": "the table of pre-specified comparisons in Sect.~3 of the paper",
                     "tab:ceiling": "Table~3 of the paper", "tab:bounds": "Table~4 of the paper",
                     "tab:etym": "Table~5 of the paper", "tab:reading": "Table~5 of the paper",
                     "tab:si-splits": "Table~C1 of the paper",
                     "tab:si-sensitivity": "Table~C2 of the paper"}
SUPP_NOTES = {name: note for _, items in SUPP_GROUPS for name, note in items}
# 本文（main.tex が読む節）と本文の表（B-12 の A3〜A8 で組み直した後）
PAPER_FILES = ("abstract.tex", "intro.tex", "methods.tex", "results.tex", "uncertainty.tex", "detect.tex", "limits.tex",
               "discussion.tex", "conclusions.tex", "appendix.tex", "availability.tex")
PAPER_TABLE_FILES = ("tab_ceiling.tex", "tab_bounds.tex", "tab_reading.tex")


def _uncommented(f: Path) -> str:
    return "\n".join(ln for ln in f.read_text(encoding="utf-8").splitlines() if not ln.lstrip().startswith("%"))
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


TERMS = (  # 初見の査読 第 7 回（A8）: 用語を本文に揃える（地形だけでないので baseline、yardstick は使わない）
    ("bare-earth terrain model", "bare-earth elevation model"), ("Bare-earth terrain model", "Bare-earth elevation model"),
    ("terrain-only", "baseline-only"),
    ("recalibrated terrain model", "recalibrated baseline model"), ("Terrain model", "Baseline model"),
    ("terrain model", "baseline model"), ("terrain gain", "baseline gain"), ("Terrain gain", "Baseline gain"),
    ("terrain increment", "baseline gain"), ("Terrain increment", "Baseline gain"),
    ("Terrain $", "Baseline $"),  # 「Terrain $\\Delta$BS」の見出し（\\Delta は命令として守られるので前だけを）
    ("alternative yardsticks", "alternative measures of the baseline gain"),
    ("five yardsticks", "five measures of the baseline gain"),
    ("yardsticks", "measures of the baseline gain"), ("Yardsticks", "Measures of the baseline gain"),
    ("yardstick", "measure of the baseline gain"), ("Yardstick", "Measure of the baseline gain"),
)


def _terms(sx: str) -> str:
    """表の文の用語を本文に揃える（\\label や \\ref の中の名前は変えない）。"""
    keep = {}

    def _hold(mm: re.Match) -> str:
        keep[f"@@{len(keep)}@@"] = mm.group(0)
        return f"@@{len(keep) - 1}@@"

    sx = re.sub(r"\\(?:label|ref|input)\{[^}]*\}|\\[A-Za-z]+", _hold, sx)  # 命令の名前（マクロ）も守る
    for a, b in TERMS:
        sx = sx.replace(a, b)
    for k, v in keep.items():
        sx = sx.replace(k, v)
    return sx


def _strata_digits(sx: str) -> str:
    """層の表で、R と 90 % 区間の端が同じ文字列になる行を 4 桁で書き直す（B-12 の A12、表 S8 の fa の行）。"""
    st = (((_json("v36_yardstick_robustness.json") or {}).get("splits") or {}).get("primary") or {}).get("strata") or {}

    def f4(x: float) -> str:
        return f"${x:+.4f}$".replace("+", "{+}").replace("-", "{-}")
    for k, v in st.items():
        if not k.startswith("language_") or v.get("R") is None or not v.get("ci90_R"):
            continue
        r, (lo, hi) = g._r(v["R"]), v["ci90_R"]
        if r not in (g._r(lo), g._r(hi)):
            continue
        old = f"{r} & {g._ci(v['ci90_R'])} \\\\"
        new = f"{f4(v['R'])} & [{f4(lo)}, {f4(hi)}] \\\\"
        sx = sx.replace(f"Language: {k[9:]} & {g._n(v['n_cells'])} & {100 * v['wet_rate']:.1f} & {old}",
                        f"Language: {k[9:]} & {g._n(v['n_cells'])} & {100 * v['wet_rate']:.1f} & {new}")
    return sx


def supplement() -> dict[str, str]:
    """補足の表を supp_table_*.tex に書き、付録の表の改訂の名前を直し、S 番号のマクロを返す。"""
    m: dict[str, str] = {}
    k = 0
    lines, contents, current = [], [], None
    # 初見の査読 第 7 回（編集者の M6、A7）: 本文・付録・補足の本文が引く表を先に、残りを「追加の表」に分ける
    cited: set[str] = set()
    paper = HERE.parent
    for f in [*(paper / x for x in (*PAPER_FILES, "supp_text.tex", "supplement.tex")),
              *(HERE / f"app_table_{x}.tex" for x in APPENDIX_TABLES), *(HERE / x for x in PAPER_TABLE_FILES)]:
        if f.is_file():
            txt = _uncommented(f)
            cited |= set(re.findall(r"\\(STab[A-Za-z]+)", txt))
            cited |= {_label_macro(x) for x in re.findall(r"\\ref\{(tab:[^}]*)\}", txt)}

    def _text(name: str) -> str:
        f = HERE / f"app_table_{name}.tex"
        return f.read_text(encoding="utf-8") if f.is_file() else ""

    lab_of = {n: re.findall(r"\\label\{(tab:[^}]*)\}", _text(n)) for n in SUPP_TABLES}
    name_of = {_label_macro(lab): n for n, labs in lab_of.items() for lab in labs}
    # B-12 の A10（項目 11）: 本文と補足の本文が引く表と、それらの表が引く表だけを補足に。残りは別の PDF（保管庫）に
    keep = {n for n, labs in lab_of.items() if any(_label_macro(lab) in cited for lab in labs)}
    grew = True
    while grew:
        grew = False
        for n in list(keep):
            t = _text(n)
            refs = {_label_macro(x) for x in re.findall(r"\\ref\{(tab:[^}]*)\}", t)} | set(re.findall(r"\\(STab[A-Za-z]+)", t))
            for r in refs:
                if r in name_of and name_of[r] not in keep:
                    keep.add(name_of[r])
                    grew = True
    order = [n for n in SUPP_TABLES if n in keep]
    archive = [n for n in SUPP_TABLES if n not in keep and (HERE / f"app_table_{n}.tex").is_file()]
    ka = 0
    for name in archive:  # 保管庫の表の番号（A1, A2, ...）を先に決める（補足の表から引かれることはない）
        for lab in lab_of[name]:
            ka += 1
            m[_label_macro(lab)] = f"A{ka}"
    def _prep(name: str) -> str:
        f = HERE / f"app_table_{name}.tex"
        sx = _terms(_groups(_rename_amendments(f.read_text(encoding="utf-8"))))
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
        sx = sx.replace("}\\label{tab:si-posthoc}", " Values of $R$ for subsets with a small baseline gain "
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
        # B-12 の A12: 座標が分に丸められた集落の割合の分母を書く（この解析の集落の点と、解析したセルの最寄りの集落）
        sx = re.sub(r"(\d+(?:\.\d)?\\%) of settlement points\)",
                    lambda mm: (f"{mm.group(1)} of the settlement points of this analysis, against "
                                "\\RoundedSettlementsPct{} of the settlements nearest the analysed cells)"), sx)
        if name == "strata":  # B-12 の A12: 点と区間の端が 3 桁で同じに見える行（fa）は 4 桁で
            sx = _strata_digits(sx)
        if name == "yardsticks":  # 横置きの頁に収まらず空の頁ができていた（B-12 の A10: p. 13 の空白）
            sx = sx.replace("\\footnotesize\\setlength{\\tabcolsep}{4pt}", "\\scriptsize\\setlength{\\tabcolsep}{3pt}")
        if "\\begin{sidewaystable}" in sx:  # 横置きの表は横向きの頁に（PDF の閲覧で回して読まなくてよいように）
            sx = (sx.replace("\\begin{sidewaystable}", "\\begin{landscape}\n\\begin{table}[H]")
                  .replace("\\end{sidewaystable}", "\\end{table}\n\\end{landscape}")
                  .replace("\\ifdim\\width>\\textheight\\textheight", "\\ifdim\\width>\\linewidth\\linewidth"))
        return sx

    for name in order:
        if not (HERE / f"app_table_{name}.tex").is_file():
            continue
        sx = _prep(name)
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
    supp_labels = {lab for name in order for lab in lab_of[name]}
    # 保管庫（B-12 の A10）: 引かれない表を番号 A で別の文書に。補足の表への参照は S 番号のマクロに
    alines, acontents, current = [], [], None
    for name in archive:
        sx = _prep(name)
        sx = re.sub(r"Tables?~\\ref\{(tab:[^}]*)\}", lambda mm: (f"Table~\\{_label_macro(mm.group(1))}{{}}"
                                                                if mm.group(1) in supp_labels else mm.group(0)), sx)
        (HERE / f"arch_table_{name}.tex").write_text(sx, encoding="utf-8")
        nums = [m[_label_macro(lab)] for lab in lab_of[name]]
        group = next(gt for gt, items in SUPP_GROUPS if name in dict(items))
        if group != current:
            current = group
            alines.append(f"\\section*{{{group}}}")
            acontents.append(f"\\multicolumn{{2}}{{l}}{{\\emph{{{group}}}}} \\\\")
        alines.append(f"\\input{{figures/arch_table_{name}}}")
        acontents.append(f"{'--'.join(dict.fromkeys([nums[0], nums[-1]]))} & {SUPP_NOTES[name]} \\\\")
    if archive:
        contents.append(f"A1--A{ka} & Tables that neither the paper nor the text above cites, kept for completeness in a "
                        "separate file, \\texttt{supplement\\_archive.pdf}, in the public repository, which lists them. "
                        "\\\\")
    (HERE / "supp_tables.tex").write_text("% generated by figures/make_figures.py — do not edit\n" + "\n".join(lines)
                                          + "\n", encoding="utf-8")
    head = "\\begin{longtable}{lp{14.5cm}}\n\\toprule\nTable & What it shows \\\\\n\\midrule\n\\endhead\n"
    tail = "\n\\bottomrule\n\\end{longtable}\n"
    gen = "% generated by figures/make_figures.py — do not edit\n"
    (HERE / "supp_contents.tex").write_text(gen + head + "\n".join(contents) + tail, encoding="utf-8")
    (HERE / "archive_tables.tex").write_text(gen + "\n".join(alines) + "\n", encoding="utf-8")
    (HERE / "archive_contents.tex").write_text(gen + head + "\n".join(acontents) + tail, encoding="utf-8")
    for name in APPENDIX_TABLES:  # 付録に残す表から補足の表への参照は S 番号のマクロに
        f = HERE / f"app_table_{name}.tex"
        if f.is_file():
            sx = _terms(_groups(_rename_amendments(f.read_text(encoding="utf-8"))))
            sx = re.sub(r"Tables?~\\ref\{(tab:[^}]*)\}",
                        lambda mm: (f"Table~\\{_label_macro(mm.group(1))}{{}} of the Supplement"
                                    if mm.group(1) in supp_labels else mm.group(0)), sx)
            f.write_text(sx, encoding="utf-8")
    return m


def table_hazard() -> None:
    """本文の表 2: 事前に固定した R と、名前の地面が全体の R に寄与しうる量・オッズ比への翻訳・対照・直接の推定。"""
    # 初見の査読 第 6 回: 事前に固定した量（R）は本文と表 2 にあるので外し、この表は探索的な量だけにする
    rows = [
        r"\multicolumn{3}{l}{\emph{What this ground contributes to $R$, on the same pooled scale}} \\",
        (r"Observed contribution (90\,\% interval)$^{a}$ & \ObsRAllWet{} (\ObsRAllWetLo{} to \ObsRAllWetHi) & "
         r"\ObsRFiveHundred{} (\ObsRFiveHundredLo{} to \ObsRFiveHundredHi) \\"),
        r"Largest contribution if names predicted this ground perfectly & \CeilRAll & \CeilRFiveHundred \\",
        r"\multicolumn{3}{l}{\emph{Odds ratios of flooding from effects injected into labels}} \\",
        r"Odds ratio at the margin, observed labels, base models fixed & \PlasMarginA & \PlasMarginB \\",
        r"\quad same, synthetic labels of the power analysis & \ORAtMargin & \ensuremath{>}\,20 \\",
        (r"Upper bound from the observed increment, fitted models held fixed$^{f}$ & \BoundNewA{} "
         r"(\BoundNewALo--\BoundNewAHi) & \BoundNewB{} (\BoundNewBLo--\BoundNewBHi) \\"),
        r"\quad with the spread over fold assignments & \BoundFoldA & \BoundFoldB \\",
        r"\quad with the second stage refitted within the bootstrap & \BoundRefA & \BoundRefB \\",
        r"\quad same as the first, for truly water-related names$^{b}$ & \BoundNewCA & \BoundNewCB \\",
        r"Lower bound for effects that remove flooding$^{f}$ & \BoundNegA & \BoundNegB \\",
        r"Upper bound, flips correlated within the ground of each settlement$^{f}$ & \BoundACA & --- \\",
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
           "settlement. All quantities are exploratory (amendments 2, 3 and 5 to 10, or not planned, footnote a); the "
           "pre-specified result is in Table~\\ref{tab:si-comparisons}. Contributions are on the scale of the pooled $R$, "
           "that is, the change in the Brier "
           "score of all cells divided by the terrain gain. $^{a}$~The increments on this ground (amendment 2) "
           "re-expressed on the pooled scale (not planned). $^{b}$~Effect injected on a latent indicator of truly "
           "water-related names, drawn with the misreading rates measured against recorded etymologies. $^{c}$~Leaks "
           "the outcome by construction (the cell itself is left out): an upper bound, not a possible reading of names. "
           "$^{d}$~Known for \\RecurCovPct{} of the cells; zero with an indicator elsewhere. $^{e}$~AOIs and cells that "
           "radar can assess; built-up land and seasonal water are excluded. $^{f}$~Base models fixed; no-effect reference: the "
           "mean increment of reassigned names with the regional composition rebuilt (amendment 8); in parentheses, Monte "
           "Carlo 95\\,\\% intervals over the simulated replicates (\\PlasRepsA{} each). "
           "Bounds with other references are in Tables~\\STabEightBounds{} and \\STabNineBounds{} of the Supplement. "
           "$^{g}$~Within 500\\,m, intervals from a pairs cluster bootstrap over location groups (2,000 replicates); "
           "other odds ratios have cluster-robust intervals.}\\label{tab:hazard}\n\\scriptsize\n"
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
    # 初見の査読 第 6 回（編集者の 4 節の 7）: 市街地の層も含め、プロトコルの全ての層の範囲
    h5 = [v["R_within"] for k, v in (met.get("H5_strata") or {}).items()
          if isinstance(v, dict) and v.get("R_within") is not None]
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
    m["FabText"] = (f"with a bare-earth elevation model \\citep{{hawker2022}} and heights above the main rivers in "
                    f"place of HAND ({_rtext(c7, 'v64 C7')})" if isinstance(c7, dict) else missing("v64 C7"))
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
    cap = ("Superseded in part by Table~\\STabNineBounds{}: rows that use the withdrawn reassignments as reference are "
           "kept for the record. Odds ratios of flooding from effects injected into the observed labels, read with "
           "different no-effect "
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
                     ("cleaned_reading", "Cleaned reading"), ("uncapped", "Every evaluable cell (no cap; frozen-code readings)")):
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
           "Columns 4--6: flooded cells within 500\\,m of water-side names, location groups with any cell within "
           "500\\,m of a water-side name (flooded or not), and the AOIs that enter each model (AOIs with flooded and dry "
           "cells). The uncapped row uses the readings of the frozen code; the corrected readings are in the table of "
           "uncapped estimates of amendment 11.")
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
    for key, lab in (("M1", "Terrain model (M1)"),):
        e = b9.get(key)
        if not e:
            continue
        q = lambda d: f"{d['median']:.2f} ({d['q25_q75'][0]:.2f}--{d['q25_q75'][1]:.2f})"
        rows.append(f"\\multicolumn{{2}}{{l}}{{{lab}}} & \\multicolumn{{6}}{{l}}{{CSI {q(e['csi'])}; hit rate "
                    f"{q(e['hit_rate'])}; false alarms {q(e['false_alarm_ratio'])}; AUC {q(e['auc_within'])}; "
                    f"{e['n_aoi']} AOIs}} \\\\")
    hb = _hand_b15()  # B-12 の B15: HAND だけの地図はセルの平均 HAND（直した実行）
    if hb:
        rows.append(f"\\multicolumn{{2}}{{l}}{{HAND alone (mean HAND of the cell; corrected run)}} & "
                    f"\\multicolumn{{6}}{{l}}{{CSI {hb['csi_median']:.2f} ({hb['csi_q25_q75'][0]:.2f}--"
                    f"{hb['csi_q25_q75'][1]:.2f}); hit rate {hb['hit_median']:.2f}; false alarms {hb['far_median']:.2f}; "
                    f"{hb['n_aoi']} AOIs}} \\\\")
    cap = ("Calibration of the held-out predictions (exploratory, amendment 8): intercept and slope of a logistic "
           "regression of the labels on the logit of the prediction, Brier score with the reliability and resolution "
           "terms of Murphy's decomposition (ten equal-count bins), and median within-AOI area under the ROC curve. "
           "M1: terrain model; M1s: recalibrated; M1c: with the regional composition of names (reference of the primary "
           "comparison); M3n: with the names. Skill: each AOI predicted flooded on as many cells as the labels flood "
           "(CSI: critical success index). The primary, country, time and language-model splits use the confirmatory run "
           "(frozen code); the Brier scores of the corrected run, which differ in the fourth decimal, are in "
           "Table~\\ref{tab:ten-splits}.")
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
    rows.append(f"Names read as water-side with a water-related derivation in coder 2's codes: transferred or "
                f"commemorative; feature not at the site & {tr['yes']} of {tr['n']}; {st['no']} of {st['n']} \\\\")
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
           "(amendment 8, addendum 6). AOIs: the AOIs with such cells whose nearest settlement lies in the country or "
           "activation; an AOI whose cells are nearest settlements in several countries counts in each, so the column "
           f"sums to more than the number of AOIs ({sum(r['aois'] for r in v60.get('by_country') or [])} by country).")
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
           "(base models fixed), with the spread over 13 assignments of location groups to folds added, or with the "
           "second stage refitted within each bootstrap replicate. Reference: reassigned names with rebuilt "
           "composition unless stated. Detectable: the odds ratio whose mean gain exceeds the reference by "
           "2.80 standard errors of the increment (two-sided 5\\,\\%, power 80\\,\\%). Simulated odds ratios: "
           f"{', '.join(f'{np.exp(x):.2f}' for x in grid_a if x >= 0.25 - 1e-9)}{low_txt}; "  # noqa: PLR2004
           f"{b['reps']['A']} replicates per effect and target; values beyond them are marked $>$ or $<$. In "
           "parentheses: 2.5th and 97.5th percentiles over 2,000 resamples of the replicates (Monte Carlo interval); "
           "the resampling is drawn afresh for this table, so a limit can differ in the last digit from the same "
           "interval in Table~\\ref{tab:eight-bounds}."
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
        hb = _hand_b15()  # B-12 の B15: HAND だけの地図はセルの平均 HAND（四分位は CSI だけ）
        hk = {"csi": "csi", "hit_rate": "hit", "false_alarm_ratio": "far"}[key]
        hq = (f" ({hb['csi_q25_q75'][0]:.2f}--{hb['csi_q25_q75'][1]:.2f})" if key == "csi" else "")
        rows.append(f"\\quad {lab} & {cells[0]} & {cells[1]} & "
                    + (f"mean HAND: {hb[hk + '_median']:.3f}{hq}" if hb else "---") + " \\\\")
    dc = b5["delta_csi_M3n_minus_M1c"]
    rows.append(f"\\quad Difference in the critical success index & & & median {dc['median']:.3f} "
                f"({dc['ci95_median'][0]:.3f} to {dc['ci95_median'][1]:.3f}) \\\\")
    cap = ("The skill of the held-out predictions of the corrected run with and without the names, decomposed "
           "(exploratory; amendment 9). Murphy decomposition (Murphy, 1973) over ten bins of predicted "
           "probability (the binned form omits the within-bin variance and covariance terms, so reliability minus "
           "resolution plus uncertainty does not reproduce the Brier score exactly); within AOIs: each model's predictions shifted on the logit scale so that its mean in each AOI "
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
        rows.append(("European AOIs alone, own baseline gain", "7", "corrected", v55["europe"]["europe_lexfix"]))
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
    status = [  # 節の番号は B-11 の構成（初見の査読 第 7 回、A4）
        ("H2 verdict, $R$, readings, placebo; protocol's two-way bootstrap", P, "protocol", "3.2"),
        ("H1, H3, H4, H5; total and flexible increments; secondary splits; seven sensitivity analyses", P,
         "protocol", "3.3, C"),
        ("$p$ values of H2--H4 (computed after the critiques)", P, "6, add.~1", "3.3, S"),
        ("Rerun with the corrected reading of Norwegian names; corrected sensitivity analysis 6", CO, "4; 1",
         "2.1, 3.2, C"),
        (("Positive control on all cells; alternative measures of the baseline gain; two-way bootstrap variant; "
          "distance bands from 500\\,m; within-AOI odds ratios on all cells; fold assignments; events held out whole"),
         PH, "2", "3.3, 3.4, 4.2"),
        (("Class split of the positive control; overlap of top deciles; bound translations under other fold "
         "assignments; split-specific power ratios; interval coverage; class intervals; real minus permuted by "
         "class; operating characteristics against the real margin; interval combining fold assignments"), CB, "2",
         "S"),
        (("Positive control at settlement points and by distance band; flood contrasts and increments below "
         "500\\,m; odds ratios by class and band and within 500\\,m of water-side names"), NP, "2", "5, 6.3"),
        ("GFM labels, designated replication, rescoring, odds ratios", PH, "3", "5.2"),
        ("Check of GFM masking (planned after the GFM result had been seen)", NP, "3, add.~2", "5.2, S"),
        (("Permutation distribution; coefficients; thirteen fold assignments; test inversion; terrain oracle; "
         "within-area metrics; label variants; element groups (expected directions written after the odds ratios of "
         "three classes near settlements had been computed); settler societies excluded; etymology checks"), PH,
         "5", "3.4, 5.3, 6.2"),
        (("Ceiling; plasmode with fixed base models; flood-information controls; GFM diagnostics and V1; strict "
         "reading; sensitivity runs; AOI weighting; rounded coordinates; ladder rung T0$^-$"), PH, "6",
         "2.1, 4, 5, 6"),
        (("Permutations with recomputed composition; plasmode with refitted models (bound withdrawn in amendment 8); "
          "activation split with the corrected code; recurrence control; floodplain reading (defined after the class "
          "estimates had been seen); partly flooded cells; magnitude terciles; built-up and tree adjustment; GFM "
          "within-AOI discrimination; pairs bootstrap; descriptions (a)--(k)"), PH, "7", "3.4, 4, 5, A"),
        ("Joint tests of the near-settlement coefficients; digraph matching; recall by language", PH, "7, add.~1",
         "5.3, 6, S"),
        ("Joint test of the lexicon coefficients on the held-out predictions", PH, "7, add.~3", "3.4"),
        ("Observed increment within 500\\,m on the pooled scale (re-expression of an amendment-2 quantity)", NP,
         "--", "4.1"),
        (("Achieved equivalence bound; location groups needed; permutation test; second stage refitted in the "
         "bootstrap; operating characteristics at the margin; within-500\\,m increment on its own scale"), PH, "8",
         "3.2, 4.1, S"),
        (("No-effect reference of the bounds (fixed after bounds under earlier references had been seen) and their "
          "Monte Carlo intervals; more plasmode replicates; recurrent surface water of each cell; cleaned lexicon"),
         PH, "8", "4"),
        (("Estimates near settlements by class and variant (uncapped sample, misclassification, river size, the two "
         "label sets on the same cells); heights above main rivers; calibration and skill; audit of matched forms; "
         "reading by nearby river and lake names; second coder"), PH, "8", "5, 6, S"),
        (("Bounds with refitted intervals and over fold assignments; small and protective injected effects; effects "
          "that differ between settlements; simulated power with every model refitted; flood-mapping skill; one "
          "definition of permanent water; other subsampling seeds; bare-earth elevation model"), PH, "9",
         "3.2, 4.2, 5, S"),
        ("Description of H4 as carrying little information (a reading the protocol does not provide)", NP, "9",
         "3.3, A"),
        (("Detectable effects and further intervals near settlements; misreading over a grid of rates; the "
          "islands-and-meadows group at settlement sites; the floodplain reading by country and form; sensors; the "
          "activation split in full; distance and continent interactions; tuned penalty; country indicators; "
          "clustered injection; flood hazard maps in the baseline; areas within language regions; languages read; "
          "classes of the language-model reading"), PH, "10", "3.4, 4.2, 5, 6"),
        ("Placebo comparison of the GFM within-area evaluation (planned after that evaluation had been seen)", PH,
         "10, add.~2", "5.2"),
        (("Change of the activation split by location group; uncapped odds ratios near settlements on both label "
          "sets (made the main estimates after both versions had been seen); bounds on the estimator's scale; the GFM "
          "domain; the denominator of the labels; Nadeau--Bengio and combined intervals; census of exploratory "
          "estimates; feature codes and rounded coordinates; pilots; audit sample; permutations with rebuilt "
          "composition on the country split; digraph breakdown"), PH, "11", "2, 3.3, 3.4, 4.2, 5, 6, S"),
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
            "analysis; the simulations on the observed labels (Table~\\ref{tab:bounds}) give lower values.", "lll",
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
            "addendum 1). Recall is estimated as a domain with the design weights of the whole sample (stratum size over sample size), whereas Table~\\ref{tab:ten-valid} weights the strata by their sizes within each country, with "
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


# %% amendment B-10（初見の査読 第 6 回）: 本文の数と補足の表
P1_MEASURED_KEY, P0_MEASURED_KEY = "truth_rate_in_wet", "truth_rate_in_not"
DIL_PHI = (0.05, 0.1, 0.25, 0.5)
DIL_ORN = (1.5, 2.0, 3.0, 5.0)


def _pctx(x: float) -> str:
    """R の百分率（−1.1 % など、記号つきの数式）。"""
    return f"\\ensuremath{{{100 * x:.1f}\\%}}"


def _r3s(x: float | None) -> str:
    return missing("value") if x is None else rv(x)


def _rci(c: dict | None, what: str, key: str = "ci90_R") -> str:
    """R と区間の短い文（−0.004, 90 % interval −0.010 to 0.002）。"""
    if not isinstance(c, dict) or c.get("R") is None:
        return missing(what)
    lo, hi = c[key]
    lev = "90" if key.startswith("ci90") else "95"
    return f"{rv(c['R'])} ({lev}\\,\\% interval {rv(lo)} to {rv(hi)})"


def dilution(or_n: float, phi: float, p1: float, p0: float) -> float:
    """名の指す地面（最寄りの地面の割合 φ）だけでオッズ比 OR_n の効果を、名に最も近い全ての地面で、読み違い（p1・p0）を
    入れて観測したときのオッズ比（稀な結果の近似。B-10 追記 1 の 3）。"""
    rr = 1 + phi * (or_n - 1)
    return (p1 * rr + 1 - p1) / (p0 * rr + 1 - p0)


def _base_rates() -> dict[str, str]:
    """解析したセルの浸水の割合と、AOI ごとの浸水の割合の分布（確証の実行の held-out の表）。"""
    files = sorted((DATA / "runs" / "world_v2" / "predictions" / "primary").glob("fold[0-9].npz"))
    if not files:
        return {k: missing("predictions") for k in ("BaseRatePct", "AoiRateMed", "AoiRateQOne", "AoiRateQThree")}
    y = np.concatenate([np.load(f)["y"] for f in files]).astype(np.float64)
    aoi = np.concatenate([np.load(f)["aoi"] for f in files])
    r = pd.Series(y).groupby(aoi).mean()
    return {"BaseRatePct": pct1(y.mean()), "AoiRateMed": pct1(float(r.median())),
            "AoiRateQOne": pct1(float(r.quantile(0.25))), "AoiRateQThree": pct1(float(r.quantile(0.75)))}


def b10_macros(v66: dict | None, v67: dict | None, v55: dict | None, v56: dict | None, v57: dict | None,
               v44: dict | None, v48: dict | None, v61: dict | None, v64: dict | None, m: dict[str, str]) -> dict[str, str]:
    out: dict[str, str] = {}
    v66 = v66 or {}
    met = g._json(DATA / "runs" / "world_v2" / "metrics_primary.json") or {}
    pr = met.get("primary") or {}
    out |= ({"RLoNinetyPct": _pctx(pr["ci90_R"][0]), "RHiNinetyPct": _pctx(pr["ci90_R"][1])} if pr.get("ci90_R")
            else {k: missing("primary R") for k in ("RLoNinetyPct", "RHiNinetyPct")})
    lf = (g._json(DATA / "runs" / "world_v2_lexfix" / "metrics_primary.json") or {}).get("primary") or {}
    out["AchievedBoundLexfix"] = rv(max(abs(lf["ci90_R"][0]), abs(lf["ci90_R"][1]))) if lf.get("ci90_R") \
        else missing("lexfix R")
    n_aoi = int(m.get("NAoi", "0") or 0)
    for key, src in (("NoHydroPct", "NoHydroNExcl"), ("UntypedPct", "NUnspecifiedRiverine")):
        out[key] = pct1(int(m[src]) / n_aoi) if n_aoi and m.get(src, "").isdigit() else missing(key)
    out |= _base_rates()
    rm = g.real_margin_rate(DATA / "runs", "", "primary")
    out["PowerEquivReal"] = pct1(rm[0] / rm[1]) if rm else missing("real margin power")
    # B-7 の 11(g): 同じ activation の AOI の残差の相関、4: 他の事象のある AOI の部分集合
    rb = (v55 or {}).get("residual_between_aois") or {}
    sa, da = rb.get("same_activation") or {}, rb.get("different_activation") or {}
    out |= ({"ResActCorrNear": f"\\ensuremath{{{sa['0-50km']['corr']:.2f}}}",
             "ResActNPairsNear": str(sa["0-50km"]["n_pairs"]),
             "ResActCorrAll": f"\\ensuremath{{{sa['all']['corr']:.2f}}}", "ResActNPairs": str(sa["all"]["n_pairs"]),
             "ResActCorrOther": f"\\ensuremath{{{da['all']['corr']:.3f}}}"} if sa and da
            else {k: missing("B-7 11(g)") for k in ("ResActCorrNear", "ResActNPairsNear", "ResActCorrAll",
                                                     "ResActNPairs", "ResActCorrOther")})
    sub = (((v56 or {}).get("recur") or {}).get("subset_aois_with_other_events") or {})
    rs = sub.get("recur_settlement") or {}
    out |= ({"RecurSubNAoi": str(sub["n_aoi"]), "RecurSubR": rv(rs["R"]), "RecurSubRLoNinety": rv(rs["ci90_R"][0]),
             "RecurSubRHiNinety": rv(rs["ci90_R"][1]), "RecurSubYard": sci(sub["yardstick"], 2)} if rs
            else {k: missing("B-7 4 subset") for k in ("RecurSubNAoi", "RecurSubR", "RecurSubRLoNinety",
                                                        "RecurSubRHiNinety", "RecurSubYard")})
    b8 = (v57 or {}).get("B8") or {}
    gr = b8.get("gfm_refit") or {}
    out |= {"CalibGfmMOne": _f2(gr["M1"]["slope"]) if gr.get("M1") else missing("B8 gfm M1"),
            "CalibGfmMOnes": _f2(gr["M1s"]["slope"]) if gr.get("M1s") else missing("B8 gfm M1s")}
    tw = (v44 or {}).get(P1_MEASURED_KEY) or {}
    out |= ({"WPPrecLexK": str(tw["k"]), "WPPrecLexN": str(tw["n"])} if tw
            else {k: missing("v44") for k in ("WPPrecLexK", "WPPrecLexN")})
    # B12: HAND だけの参照と地図の技能（AOI の中央値と全てまとめた値）
    b12 = v66.get("B12") or {}
    out |= ({"MapMOneCsiPooled": _f2(b12["M1"]["csi_pooled"]), "HandMeanCsiMed": _f2(b12["HAND_mean"]["csi_median"]),
             "HandMeanCsiPooled": _f2(b12["HAND_mean"]["csi_pooled"]),
             "HandMinCsiMed": _f2(b12["HAND_min"]["csi_median"]),
             "HandMinCsiPooled": _f2(b12["HAND_min"]["csi_pooled"])} if b12
            else {k: missing("B12") for k in ("MapMOneCsiPooled", "HandMeanCsiMed", "HandMeanCsiPooled",
                                               "HandMinCsiMed", "HandMinCsiPooled")})
    # B1: 島・水辺の草地
    b1 = v66.get("B1") or {}
    for key, grp in (("IslandSiteHand", "island"), ("IslandNoAuSiteHand", "island_noau")):
        d = ((b1.get(grp) or {}).get("site_hand_mean")) or {}
        out |= ({key: f"\\ensuremath{{{d['estimate']:.1f}}}", f"{key}Lo": f"\\ensuremath{{{d['ci95'][0]:.1f}}}",
                 f"{key}Hi": f"\\ensuremath{{{d['ci95'][1]:.1f}}}"} if d
                else {f"{key}{s}": missing("B1") for s in ("", "Lo", "Hi")})
    d = ((b1.get("island_noau") or {}).get("site_hand_mean")) or {}
    out["IslandNoAuSiteNAoi"] = str(d["n_aoi"]) if d else missing("B1")
    o = (b1.get("island_noau") or {}).get("or_le500m") or {}
    out |= or_macros("IslandNoAuOR", o or None)
    out["IslandNoAuFlooded"] = str(o["n_flooded_le500m"]) if o else missing("B1")
    out["IslandNoAuFloodedUnits"] = (str(o["n_units_flooded_le500m"]) if o.get("n_units_flooded_le500m") is not None
                                     else missing("B1 units"))
    # B3・B4: 主の 500 m の推定の検出力と区間
    b3, b4 = v66.get("B3") or {}, v66.get("B4") or {}
    out |= ({"MdeNearUp": orv(b3["pairs"]["mde_or_up"]), "MdeNearDown": orv(b3["pairs"]["mde_or_down"])} if b3
            else {k: missing("B3") for k in ("MdeNearUp", "MdeNearDown")})
    jk, wb = (b4.get("jackknife_cv3") or {}), (b4.get("wild_score_bootstrap") or {})
    out |= ({"JkNearLo": orv(jk["ci95"][0]), "JkNearHi": orv(jk["ci95"][1])} if jk.get("ci95")
            else {k: missing("B4 jackknife") for k in ("JkNearLo", "JkNearHi")})
    out |= ({"WildNearLo": orv(wb["ci95"][0]), "WildNearHi": orv(wb["ci95"][1])} if wb.get("ci95")
            else {k: missing("B4 wild") for k in ("WildNearLo", "WildNearHi")})
    # 集落の近くの推定の幅（標本・ラベル・読みの変種、activation を一つずつ除く）
    vals = []
    b6 = (v57 or {}).get("B6") or {}
    vals += [b6[k]["odds_ratio"] for k in ("typed_aois", "without_no_hydro", "without_oceania_americas",
                                           "without_australia", "cleaned_reading", "uncapped") if b6.get(k)]
    vals += [(v61 or {})[k]["odds_ratio"] for k in ("without_pakistan", "without_pakistan_bangladesh")
             if (v61 or {}).get(k)]
    rnd = ((((v48 or {}).get("rounded") or {}).get("ems_all") or {}).get("not_rounded") or {}).get("le500m")
    if rnd:
        vals.append(rnd["odds_ratio"])
    for k in ("C5_near", "C7_near"):
        if isinstance((v64 or {}).get(k), dict) and (v64 or {})[k].get("odds_ratio"):
            vals.append(v64[k]["odds_ratio"])
    out |= ({"NearVarMin": orv(min(vals)), "NearVarMax": orv(max(vals))} if vals
            else {k: missing("near variants") for k in ("NearVarMin", "NearVarMax")})
    # B5: 誤分類の補正の格子（点推定の値で、定まる点の範囲）
    b5 = v66.get("B5") or {}
    pts = [r["or_corrected"] for r in b5.get("rows", []) if r["or_obs"] == "point" and r["or_corrected"] is not None]
    out |= ({"MisGridMin": orv(min(pts)), "MisGridMax": orv(max(pts))} if pts
            else {k: missing("B5") for k in ("MisGridMin", "MisGridMax")})
    # 区間の端を補正した値の最も広いもの（pairs の 95% の両端を同じ式で）
    los = [r["or_corrected"] for r in b5.get("rows", []) if r["or_obs"] == "pairs_lo" and r["or_corrected"] is not None]
    his = [r["or_corrected"] for r in b5.get("rows", []) if r["or_obs"] == "pairs_hi" and r["or_corrected"] is not None]
    out |= ({"MisGridLoMin": orv(min(los)), "MisGridHiMax": orv(max(his))} if los and his
            else {k: missing("B5") for k in ("MisGridLoMin", "MisGridHiMax")})
    # B2: 氾濫原の読み
    b2 = v66.get("B2") or {}
    loco = [v["odds_ratio"] for v in (b2.get("leave_one_country_out") or {}).values() if v.get("odds_ratio")]
    out |= ({"FpLocoMin": orv(min(loco)), "FpLocoMax": orv(max(loco))} if loco
            else {k: missing("B2") for k in ("FpLocoMin", "FpLocoMax")})
    ne = b2.get("without_EMSR097") or {}
    fpm = b2.get("main") or {}
    out["FpNoEmsrNineSeven"] = orv(ne["odds_ratio"]) if ne.get("odds_ratio") else missing("B2 EMSR097")
    out["FpEmsrNineSevenPct"] = (pct1(ne["n_flooded_near_in_EMSR097"] / fpm["n_group_flooded"])
                                 if ne and fpm.get("n_group_flooded") else missing("B2 EMSR097 share"))
    ex = b2.get("reread_without_forms") or {}
    out |= ({"FpExcl": orv(ex["odds_ratio"]), "FpExclLo": orv(ex["pairs_ci95"][0]),
             "FpExclHi": orv(ex["pairs_ci95"][1])} if ex.get("odds_ratio")
            else {k: missing("B2 excl") for k in ("FpExcl", "FpExclLo", "FpExclHi")})
    # 主・国を一つずつ除く・EMSR097 を除く・形を除くの全ての pairs 区間の下端の最大（全てが 1 を含むかを本文で）
    los = [v["pairs_ci95"][0] for v in [fpm, ne, ex, *(b2.get("leave_one_country_out") or {}).values()]
           if v.get("pairs_ci95")]
    out["FpPairsLoMax"] = orv(max(los)) if los else missing("B2 pairs")
    # B6: 国ごとの精度（一人目の符号）
    b6v = ((v66.get("B6") or {}).get("coder1") or {}).get("by_country") or {}
    pc = [v["precision"] for v in b6v.values() if v.get("precision") is not None]
    out |= ({"CtryPrecMin": pct1(min(pc)), "CtryPrecMax": pct1(max(pc)), "CtryNWord": _word(len(pc))} if pc
            else {k: missing("B6") for k in ("CtryPrecMin", "CtryPrecMax", "CtryNWord")})
    # B15: N-read の区分（名の数で）
    b15 = (v66.get("B15") or {}).get("by_names") or {}
    out |= ({"NReadNonEnvNamesPct": pct1(b15["non_environmental"]["share"]),
             "NReadWetlandNamesPct": pct1(b15["wetland"]["share"])} if b15
            else {k: missing("B15") for k in ("NReadNonEnvNamesPct", "NReadWetlandNamesPct")})
    # B13: 言語の地域の中の地域の間の関係（三つの言語をまとめた係数）
    b13 = (v66.get("B13") or {}).get("pooled") or {}
    out |= ({"BetweenCoef": f"\\ensuremath{{{b13['coef'][0]:.1f}}}",
             "BetweenCoefLo": f"\\ensuremath{{{b13['ci95'][0][0]:.1f}}}",
             "BetweenCoefHi": f"\\ensuremath{{{b13['ci95'][0][1]:.1f}}}",
             "BetweenNAoi": str(b13["n_aoi"]), "BetweenNAct": str(b13["n_activations"])} if b13
            else {k: missing("B13") for k in ("BetweenCoef", "BetweenCoefLo", "BetweenCoefHi", "BetweenNAoi",
                                               "BetweenNAct")})
    # B14: 他の言語で読んだセル
    b14 = v66.get("B14") or {}
    if b14:
        n_all = b14["n_cells_total"]
        cl = b14.get("country_level") or {}
        rg = {k: v for k, v in (b14.get("regional") or {}).items() if v["n_cells"] > 0}
        tot = sum(v["n_cells"] for v in cl.values()) + sum(v["n_cells"] for v in rg.values())
        top = sorted([(v["n_cells"], k) for k, v in cl.items()] + [(v["n_cells"], k) for k, v in rg.items()],
                     reverse=True)[:4]
        names = {"BD": "Bangladesh", "LK": "Sri Lanka", "PK": "Pakistan", "IN": "India", "MM": "Myanmar"}
        out["OtherLangPct"] = pct1(tot / n_all)
        out["OtherLangText"] = "most of them in " + ", ".join(names.get(k, k) for _, k in top[:-1]) + \
            f" and {names.get(top[-1][1], top[-1][1])}" if len(top) > 1 else (names.get(top[0][1], top[0][1]) if top
                                                                               else "none")
    else:
        out |= {k: missing("B14") for k in ("OtherLangPct", "OtherLangText")}
    # B7・B8: センサー・アフリカ
    b7 = v66.get("B7") or {}
    for key, s in (("RSensSOne", "s1"), ("RSensSar", "sar"), ("RSensOpt", "optical"), ("RSensNone", "none")):
        out[key] = _rci(((b7.get("sensor") or {}).get(s) or {}).get("R"), f"B7 {s}", "ci90")
        o = ((b7.get("sensor") or {}).get(s) or {}).get("or_le500m") or {}
        okey = key.replace("RSens", "ORSens")
        out |= ({okey: orv(o["odds_ratio"]), f"{okey}Lo": orv(o["pairs_ci95"][0]), f"{okey}Hi": orv(o["pairs_ci95"][1])}
                if o.get("odds_ratio") and o.get("pairs_log_or_sd", 0.0) <= PAIRS_UNSTABLE_SD
                else {f"{okey}{x}": missing(f"B7 OR {s}") for x in ("", "Lo", "Hi")})
    out["RNoAfrica"] = _rci((v66.get("B8") or {}).get("without_africa"), "B8", "ci90")
    # B9: 活性化をまとめた分割の Brier
    bs = (((v66.get("B9") or {}).get("brier_by_split") or {}).get("activation_merged") or {}).get("brier") or {}
    out |= ({"ActBSMOnec": f"{bs['M1c']:.4f}", "ActBSMOnes": f"{bs['M1s']:.4f}"} if {"M1c", "M1s"} <= bs.keys()
            else {k: missing("B9") for k in ("ActBSMOnec", "ActBSMOnes")})
    # B9: 物差し 5 通り × bootstrap 2 通りのうち、同等の成り立つ組の数
    act = (v66.get("B9") or {}).get("actfix_yardsticks") or {}
    eq = [((act.get(b) or {}).get("primary") or {}).get(y, {}).get("equivalent")
          for b in ("cluster", "two_way") for y in ("Y0", "Y1", "Y2", "Y3", "Y4")]
    out["ActFixEqN"] = (_word(sum(bool(e) for e in eq)) if all(e is not None for e in eq)
                        else missing("B9 yardsticks"))
    # 活性化の分割で M1c が M1s より悪い分のうち、名前（本物・置換、直したコードと凍結したコード）が取り戻した割合の最大
    ab, fr = act.get("brier") or {}, (((v66.get("B9") or {}).get("brier_by_split") or {})
                                      .get("activation_merged_frozen") or {}).get("brier") or {}
    shares = [(b["M1c"] - b[m]) / (b["M1c"] - b["M1s"]) for b, m in ((ab, "M3n"), (ab, "M3np"), (fr, "M3n"))
              if {"M1c", "M1s", m} <= b.keys() and b["M1c"] > b["M1s"]]
    out["ActRecoveredPct"] = pct1(max(shares)) if len(shares) == 3 else missing("B9 recovered")
    # B10: GFM の AOI の中の評価
    b10 = v66.get("B10") or {}
    out |= ({"GfmWithinR": rv(b10["R_within"]), "GfmWithinRLoNinety": rv(b10["ci90_R"][0]),
             "GfmWithinRHiNinety": rv(b10["ci90_R"][1])} if b10.get("R_within") is not None
            else {k: missing("B10") for k in ("GfmWithinR", "GfmWithinRLoNinety", "GfmWithinRHiNinety")})
    # 追記 2: 判定表の placebo（GFM）と、EMS の AOI の中の評価の判定
    pg, pe = b10.get("placebo") or {}, b10.get("placebo_ems_lexfix") or {}
    out |= ({"GfmWithinMinusP": rv(pg["R_minus_placebo"]), "GfmWithinMinusPLo": rv(pg["ci95_R_minus_placebo"][0]),
             "GfmWithinMinusPHi": rv(pg["ci95_R_minus_placebo"][1]), "GfmWithinRp": rv(pg["R_placebo"]),
             "GfmWithinVerdict": pg["decision"],
             "EmsWithinVerdict": pe["decision"] if pe.get("decision") else missing("B10 EMS")}
            if pg.get("R_minus_placebo") is not None
            else {k: missing("B10 placebo") for k in ("GfmWithinMinusP", "GfmWithinMinusPLo", "GfmWithinMinusPHi",
                                                       "GfmWithinRp", "GfmWithinVerdict", "EmsWithinVerdict")})
    # B11: 注入したラベルでの第 5 節の推定（形 A、名目の OR 2）
    b11 = v66.get("B11") or {}
    row = next((r for r in b11.get("rows", []) if r["form"] == "A" and abs(r["nominal_or"] - 2.0) < 1e-9), None)
    out |= ({"InjNominal": "2", "InjRealisedAll": orv(row["mean_or_ratio_to_base_all"]),
             "InjRealisedNear": orv(row["mean_or_ratio_to_base_le500m"])} if row
            else {k: missing("B11") for k in ("InjNominal", "InjRealisedAll", "InjRealisedNear")})
    # A12（追記 1 の 3）: 薄まりの換算
    p1 = (v44 or {}).get(P1_MEASURED_KEY, {}).get("rate")
    p0 = (v44 or {}).get(P0_MEASURED_KEY, {}).get("rate")
    if p1 is not None and p0 is not None:
        out |= {"DilORn": "3", "DilObsQuarter": orv(dilution(3.0, 0.25, p1, p0)),
                "DilObsTenth": orv(dilution(3.0, 0.1, p1, p0))}
    else:
        out |= {k: missing("v44 rates") for k in ("DilORn", "DilObsQuarter", "DilObsTenth")}
    # 物差しの 95 % 区間の上限（下限 DeltaTerrainLo は旧版の macro）と、低い δ の格子で読んだ上限（B-9 の C1）
    ci = pr.get("ci95_delta_M0_M1")
    out["DeltaTerrainHi"] = sci(ci[1], 2) if ci else missing("terrain gain interval")
    b12v = (v64 or {}).get("B1_B2") or {}
    bb = b12v.get("bounds") or {}
    e = ((bb.get("fixed") or {}).get("b7_recomputed") or {}).get("A")
    # 低い δ の格子（B-9 の C1）が曲線に入ったときだけ（入っていなければ旧い格子の値なので出さない）
    out["BoundLowGridA"] = _gt(e["odds_ratio_text"]) if e and b12v.get("low_ranges") else missing("v64 C1")
    # 補足の図の番号（現れる順: 図 S1 事後の確認、図 S2 地形の梯子）
    out |= {"STabFigChecks": "S1", "STabFigLadder": "S2"}
    # C（クラスタ、v67）
    out |= c10_macros(v67)
    return out


def c10_macros(v67: dict | None) -> dict[str, str]:
    from topsem.world.b10 import C_GRID  # C2 の格子（10^-4〜10^2 の 13 点）

    out: dict[str, str] = {}
    v67 = v67 or {}
    labs = []
    for key, src in (("RDist", "C1_dist"), ("RCont", "C1_cont"), ("RTuned", "C2_tuned"), ("RKn", "C3_country")):
        c = v67.get(src)
        vkey = key.replace("R", "Verdict", 1)
        if isinstance(c, dict) and c.get("R") is not None:
            out |= {key: rv(c["R"]), f"{key}LoNinety": rv(c["ci90_R"][0]), f"{key}HiNinety": rv(c["ci90_R"][1]),
                    vkey: c["decision"]}
            labs.append(c.get("decision"))
        else:
            out |= {f"{key}{s}": missing(f"v67 {src}") for s in ("", "LoNinety", "HiNinety")}
            out[vkey] = missing(f"v67 {src}")
            labs.append(None)
    # C2: 内側の交差検証が格子の端（最も強い罰則、C の最小）を選んだ fold の数
    ch = ((v67.get("stack_run") or {}).get("C2_chosen_C") or {}) if isinstance(v67.get("stack_run"), dict) else {}
    vals = [v for v in ch.values() if v is not None]
    out["TunedEdgeN"] = (_word(sum(abs(v - min(C_GRID)) < 1e-12 for v in vals)) if len(vals) == 10
                         else missing("v67 C2 chosen C"))
    cc = ((v67.get("C3_country") or {}).get("composition_after_country")
          if isinstance(v67.get("C3_country"), dict) else None)
    out |= ({"CompAfterCountry": rv(cc["R"]), "CompAfterCountryLo": rv(cc["ci95_R"][0]),
             "CompAfterCountryHi": rv(cc["ci95_R"][1])} if isinstance(cc, dict) and cc.get("R") is not None
            else {f"CompAfterCountry{s}": missing("v67 C3 composition") for s in ("", "Lo", "Hi")})
    if all(lab == "equivalent-null" for lab in labs):
        out["DistContText"] = "all kept the verdict \\emph{equivalent-null}"
    elif all(lab is not None for lab in labs):
        out["DistContText"] = "returned " + ", ".join(f"\\emph{{{lab}}}" for lab in labs) + ", respectively"
    else:
        out["DistContText"] = missing("v67 C1-C3 verdicts")
    hz = v67.get("C5_hazard")
    if isinstance(hz, dict) and hz.get("R") is not None:
        out["HazText"] = (f"the baseline gain rose from {sci(hz['yardstick_without_hazard'], 2)} to "
                          f"{sci(hz['yardstick'], 2)} (95\\,\\% interval {sci(hz['ci95_yardstick'][0], 2)} to "
                          f"{sci(hz['ci95_yardstick'][1], 2)}), and the names gave $R={rv(hz['R'])}$ (90\\,\\% interval "
                          f"{rv(hz['ci90_R'][0])} to {rv(hz['ci90_R'][1])}), which the decision table calls "
                          f"\\emph{{{hz['decision']}}}")
        same = hz.get("decision")
    else:
        out["HazText"] = missing("v67 C5")
        same = None
    ac = v67.get("C4_ac")
    out["BoundACA"] = _gt(ac["bound_text"]) if isinstance(ac, dict) and ac.get("bound_text") else missing("v67 C4")
    if same is None:
        out["RobustTenText"] = missing("v67 C5")
    elif same == "equivalent-null":
        out["RobustTenText"] = "left the verdict unchanged"
    elif same == "negative":
        out["RobustTenText"] = ("did not let the names improve the predictions (with the hazard maps they made them "
                                "slightly worse, \\emph{negative})")
    else:
        out["RobustTenText"] = f"gave, with the hazard maps, the verdict \\emph{{{same}}}"
    return out


def _tblp(label: str, caption: str, cols: str, head: str, rows: list[str]) -> str:
    """_tbl と同じ（補足の表、長い表は補足の変換で longtable にしない）。"""
    return _tbl(label, caption, cols, head, rows)


def _hand_b15() -> dict:
    """B-12 の B15: HAND だけの地図の一つの定義（セルの平均 HAND、低い順、同順位は無作為。v66 の B12、直した実行）。"""
    return ((_json("v69_b12_local.json") or {}).get("B15") or {}).get("HAND_mean") or {}


def table_ten_hand(v66: dict | None) -> str:
    b12 = (v66 or {}).get("B12") or {}
    rows = []
    for key, lab in (("HAND_mean", "HAND alone: mean HAND of the cell (ties random; the map used in every table)"),
                     ("HAND_min", "A different map: minimum HAND (ties by mean HAND, then random)"),
                     ("M1", "Terrain model (M1)"),
                     ("M1c", "Reference model (M1c)"), ("M3n", "Names model (M3n)")):
        d = b12.get(key)
        if not d:
            rows.append(f"{lab} & \\multicolumn{{6}}{{l}}{{{missing('B12')}}} \\\\")
            continue
        rows.append(f"{lab} & {_f2(d['csi_median'])} [{_f2(d['csi_q25_q75'][0])}, {_f2(d['csi_q25_q75'][1])}] & "
                    f"{_f2(d['hit_median'])} & {_f2(d['far_median'])} & {_f2(d['csi_pooled'])} & "
                    f"{_f2(d['hit_pooled'])} & {_f2(d['far_pooled'])} \\\\")
    ti = b12.get("ties_at_cut_hand_min") or {}
    cap = ("Skill of flood maps with as many mapped cells as flooded cells in each AOI (exploratory, amendment 10, "
           "item B12; corrected code; " + (f"{b12['M1']['n_aoi']}" if b12.get("M1") else "") + " AOIs with flooded "
           "and dry cells). Median over AOIs with the interquartile range, and pooled over all AOIs (true positives "
           "summed over AOIs divided by the summed union of the mapped and flooded cells). CSI: critical success "
           "index; hit: hit rate; FAR: false-alarm ratio. "
           + (f"At the cut of the minimum-HAND map a median of {pct1(ti['median_share_cells_tied'])} of the cells of "
              "an AOI are tied, so that map depends on how ties are broken; tables of earlier amendments used other "
              "tie-breaking rules and are superseded by this one." if ti else ""))
    return _tbl("tab:ten-hand", cap, "p{6.4cm}llllll",
                r"Map & CSI, median [IQR] & Hit & FAR & CSI, pooled & Hit & FAR", rows)


def _diff_ci(d: dict | None, scale: float = 1.0, fmt: str = "+.1f") -> str:
    """差と 95% 区間（{estimate, ci95}）を「+x (+a to +b)」に。"""
    if not d:
        return "---"
    return f"${scale * d['estimate']:{fmt}}$ (${scale * d['ci95'][0]:{fmt}}$ to ${scale * d['ci95'][1]:{fmt}}$)"


def table_ten_islands(v66: dict | None) -> str:
    b1 = (v66 or {}).get("B1") or {}
    rows = []
    for key, lab in (("island", "Islands and water meadows (all forms)"),
                     ("island_noau", "Without \\emph{au}, \\emph{aue} and \\emph{auen}")):
        g = b1.get(key)
        if not g:
            rows.append(f"{lab} & \\multicolumn{{6}}{{l}}{{{missing('B1')}}} \\\\")
            continue
        o = g.get("or_le500m") or {}
        e, c = _orc(o)
        fl = f"{o['n_flooded_le500m']}" if o else "---"
        if o.get("n_units_flooded_le500m") is not None:
            fl += f" ({o['n_units_flooded_le500m']})"
        rows.append(f"{lab} & {b1['n_settlements'][key]} & {_n(o['n_cells_le500m']) if o else '---'} & {fl} & "
                    f"{e} ({c}) & {_diff_ci(g.get('site_hand_mean'))} & "
                    f"{_diff_ci(g.get('site_hand_pct'), 100, '+.0f')} & {_diff_ci(g.get('all_hand_mean'))} \\\\")
    cap = ("The group of island and water-meadow names (exploratory; amendment 10, item B1). Odds ratio of flooding "
           "within 500\\,m of the settlement relative to the other names of the same AOI, given the held-out terrain "
           "prediction (the estimator of amendment 5, with a 95\\,\\% cluster-robust interval), with the flooded "
           "cells within 500\\,m and, in brackets, the location groups that hold them; difference in mean HAND (m) and "
           "in the mean HAND percentile within the AOI (percentage points) at the settlement sites (cells whose centre "
           "lies within 177\\,m of the settlement), and difference in mean HAND over all cells nearest the names, each "
           "against the other settlements or cells of the same AOI with a 95\\,\\% cluster-bootstrap interval. The "
           "group that amendment 5 recorded as river islands includes German \\emph{au}, \\emph{aue} and "
           "\\emph{auen}, which name water meadows; amendment 5 expected the sites of river-island settlements to "
           "stand higher than other settlements. The table of element groups counts every cell within 500\\,m of these "
           "settlements (one more than here, in an AOI without flooded cells, which the regression drops) and draws its "
           "bootstrap with another seed, so its counts and intervals differ slightly.")
    return _tbl("tab:ten-islands", cap, "p{3.6cm}rrrcccc",
                r"Group & Settlements & Cells $\leq500$\,m & Flooded & OR $\leq500$\,m & HAND at sites & "
                r"HAND pct.\ at sites & HAND all cells", rows)


PAIRS_UNSTABLE_SD = 3.0  # pairs の log OR の SD がこれを超えたら、再標本で分離が起きて区間が定まらないとみなす


def _or_subset_cells(o: dict | None) -> str:
    """部分集合の 500 m の OR の 4 列（OR、空き、95% 区間、浸水セル）。pairs が分離で壊れたらクラスタ頑健の区間を。"""
    if not isinstance(o, dict) or o.get("odds_ratio") is None:
        n = (o or {}).get("n_flooded_near")
        return f"--- & & too few flooded cells & {n if n is not None else '---'}"
    if o.get("pairs_log_or_sd", 0.0) > PAIRS_UNSTABLE_SD:
        iv = f"{orv(o['ci95'][0])}--{orv(o['ci95'][1])} (cluster-robust)$^{{c}}$"
    else:
        iv = f"{orv(o['pairs_ci95'][0])}--{orv(o['pairs_ci95'][1])}"
    return f"{orv(o['odds_ratio'])} & & {iv} & {o['n_flooded_near']} flooded cells"


def table_ten_robust(v66: dict | None, v67: dict | None) -> str:
    v66, v67 = v66 or {}, v67 or {}
    rows = [r"\multicolumn{5}{l}{\emph{Primary comparison in subsets (corrected code; base models of the full run)}} \\"]

    def row(lab: str, c: dict | None, what: str) -> str:
        if not isinstance(c, dict) or c.get("R") is None:
            return f"{lab} & \\multicolumn{{4}}{{l}}{{{missing(what)}}} \\\\"
        n = f"{c['n_aoi']} AOIs, {c['n_units']} groups" if "n_aoi" in c else ""
        ci90 = c.get("ci90") or c.get("ci90_R")
        ci95 = c.get("ci95") or c.get("ci95_R")
        return (f"{lab} & {rv(c['R'])} & {rv(ci90[0])} to {rv(ci90[1])} & {rv(ci95[0])} to {rv(ci95[1])} & "
                f"{c.get('decision', n)} \\\\")
    b7 = v66.get("B7") or {}
    for s, lab in (("s1", "EMS map from Sentinel-1"), ("sar", "EMS map from other radar"),
                   ("optical", "EMS map from optical images"), ("none", "Images not recorded")):
        rows.append(row(lab, ((b7.get("sensor") or {}).get(s) or {}).get("R"), f"B7 {s}"))
    for gk, gd in sorted((b7.get("generation") or {}).items()):
        rows.append(row(f"Product schema generation {gk}", gd.get("R"), f"B7 gen {gk}"))
    rows.append(row("Without the African AOIs", (v66.get("B8") or {}).get("without_africa"), "B8"))
    rows.append(r"\multicolumn{5}{l}{\emph{Odds ratio within 500\,m of water-side names in the same subsets$^{b}$}} \\")
    subsets = [(f"EMS map from {lab}", ((b7.get("sensor") or {}).get(s) or {}).get("or_le500m"))
               for s, lab in (("s1", "Sentinel-1"), ("sar", "other radar"), ("optical", "optical images"),
                              ("none", "images not recorded"))]
    subsets += [(f"Product schema generation {gk}", gd.get("or_le500m"))
                for gk, gd in sorted((b7.get("generation") or {}).items())]
    for lab, o in subsets:
        rows.append(f"{lab} & " + _or_subset_cells(o) + " \\\\")
    rows.append(r"\multicolumn{5}{l}{\emph{Alternative second stages and references (amendment 10, part C)}} \\")
    for src, lab in (("C1_dist", "Name features $\\times$ distance bands (M1cd $\\to$ M3nd)"),
                     ("C1_cont", "Name features $\\times$ continents (M1cr $\\to$ M3nr)"),
                     ("C2_tuned", "Penalty by grouped inner cross-validation (M1ct $\\to$ M3nt)"),
                     ("C3_country", "Country indicators in the reference (M1kc $\\to$ M3kn)"),
                     ("C5_hazard", "Terrain model with JRC flood hazard maps (M1hc $\\to$ M3hn)")):
        rows.append(row(lab, v67.get(src), f"v67 {src}"))
    c3 = (v67.get("C3_country") or {}).get("composition_after_country")
    rows.append("Regional composition after the country indicators (M1k $\\to$ M1kc), share of the terrain gain & "
                + (f"{rv(c3['R'])} & & {rv(c3['ci95_R'][0])} to {rv(c3['ci95_R'][1])} & \\\\" if c3
                   else f"\\multicolumn{{4}}{{l}}{{{missing('v67 C3 composition')}}} \\\\"))
    b10 = v66.get("B10") or {}
    rows.append(r"\multicolumn{5}{l}{\emph{Scored within AOIs, each AOI's flooded fraction given$^{a}$}} \\")
    for key, lab in (("placebo_ems_lexfix", "EMS labels, corrected run (amendment 9, item B4)"),
                     ("placebo", "GFM labels, every model refitted (amendment 10, item B10)")):
        c = b10.get(key)
        rows.append(row(lab, c, f"B10 {key}"))
        if isinstance(c, dict) and c.get("R_placebo") is not None:
            rows.append(f"\\quad reassigned names (placebo); real minus placebo & {rv(c['R_placebo'])} & & "
                        f"{rv(c['ci95_R_placebo'][0])} to {rv(c['ci95_R_placebo'][1])} & "
                        f"{rv(c['R_minus_placebo'])} ({rv(c['ci95_R_minus_placebo'][0])} to "
                        f"{rv(c['ci95_R_minus_placebo'][1])}) \\\\")
    cap = ("Robustness of the primary comparison and alternative estimands (exploratory, amendment 10). $R$ with 90\\,\\% "
           "and 95\\,\\% cluster-bootstrap intervals over location groups; in subsets the terrain gain is measured on "
           "the same subset. Sensors are those of the EMS post-event images (Sentinel-1 if any image is Sentinel-1). "
           "Part C reruns the second stage on the primary split with the corrected code and applies the decision table "
           "of the protocol, without replacing the confirmatory verdict; the JRC maps are estimated from a reanalysis "
           "that includes the event years. $^{a}$The logit of each model's held-out prediction is shifted within each "
           "AOI so that its mean equals the observed flooded fraction, and AOIs that are all flooded or all dry are "
           "left out; the terrain gain is measured in the same way. The placebo rows were added after the GFM result "
           "had been seen (addendum 2 of amendment 10), because the decision table requires them for an interval "
           "above zero; the last column gives $R$ minus the placebo's $R$ with its 95\\,\\% interval. $^{b}$The "
           "estimator of Sect.~5 (AOI fixed effects, distance bands and the logit of the terrain model's held-out "
           "prediction) on the cells of each subset, with a 95\\,\\% pairs cluster bootstrap over location groups "
           "(2,000 replicates); not estimated when fewer than five water-side cells within 500\\,m flooded. $^{c}$The "
           "pairs bootstrap separated in some replicates (standard deviation of the log odds ratio "
           f"above {PAIRS_UNSTABLE_SD:g}), so the cluster-robust interval is given.")
    return _tbl("tab:ten-robust", cap, "p{7.4cm}llll", r"Comparison & $R$ & 90\,\% interval & 95\,\% interval & "
                r"Verdict or sample", rows)


def table_ten_splits(v66: dict | None) -> str:
    b9 = (v66 or {}).get("B9") or {}
    rows = []
    for key, lab in (("main", "Location groups (primary)"), ("country", "Countries"), ("time", "Time (from 2024)"),
                     ("cutoff", "Language-model cut-off"), ("activation_merged", "Activations kept whole"),
                     ("activation_merged_frozen", "Activations kept whole, frozen code")):
        d = (b9.get("brier_by_split") or {}).get(key)
        if not d:
            rows.append(f"{lab} & \\multicolumn{{6}}{{l}}{{{missing('B9')}}} \\\\")
            continue
        br = d["brier"]
        rows.append(f"{lab} & {d['code']} & " + " & ".join(f"{br[k]:.4f}" if k in br else "---"
                                                          for k in ("M0", "M1", "M1s", "M1c", "M3n")) + " \\\\")
    ay = ((b9.get("actfix_yardsticks") or {}).get("cluster") or {}).get("yardsticks") or {}
    note = ("; on the activation split the terrain gain is " + ", ".join(f"{k} {sci(v, 2)}" for k, v in ay.items())
            if ay else "")
    act = b9.get("actfix_yardsticks") or {}
    if act:
        rows.append(r"\multicolumn{7}{l}{\emph{$R$ on the activation split (corrected code) under each measure of the "
                    r"terrain gain: 90\,\% interval, cluster bootstrap and two-way bootstrap}} \\")
        for yk, lab in (("Y0", "Y0: pre-specified"), ("Y1", "Y1: lower 95\\,\\% bound"), ("Y2", "Y2: recalibrated"),
                        ("Y3", "Y3: from a constant"), ("Y4", "Y4: stacked")):
            cells = []
            for bs in ("cluster", "two_way"):
                c = ((act.get(bs) or {}).get("primary") or {}).get(yk)
                cells.append(f"\\multicolumn{{2}}{{l}}{{{rv(c['ci90_R'][0])} to {rv(c['ci90_R'][1])}"
                             f"{'' if c['equivalent'] else '$^{\\dagger}$'}}}" if c else "\\multicolumn{2}{l}{---}")
            c0 = ((act.get("cluster") or {}).get("primary") or {}).get(yk) or {}
            rows.append(f"\\quad {lab} & & {rv(c0['R']) if c0 else '---'} & {cells[0]} & {cells[1]} \\\\")
    cap = ("Brier scores of the held-out predictions of each split (exploratory, amendment 10, item B9). Code: the "
           "frozen code of the confirmatory run or the corrected code" + note + "; Table~\\STabSiYardsticks{} defines the "
           "measures of the terrain gain and the bootstraps. On the activation split the reference model M1c is worse "
           "than the recalibrated terrain model M1s. Lower block: $R$ under each measure, with the 90\\,\\% interval "
           "of the cluster and the two-way bootstrap ($^{\\dagger}$: the interval does not lie within $\\pm0.05$).")
    return _tbl("tab:ten-splits", cap, "p{4.2cm}lrrrrr", r"Split & Code & M0 & M1 & M1s & M1c & M3n", rows)


def table_ten_injected(v66: dict | None) -> str:
    b11 = (v66 or {}).get("B11") or {}
    rows = []
    for r in b11.get("rows", []):
        rows.append(f"{r['form']} & {r['nominal_or']:.2f} & {orv(r['mean_or_ratio_to_base_all'])} & "
                    f"{orv(r['mean_or_ratio_to_base_le500m'])} & {r['sd_log_or_all']:.3f} & "
                    f"{r['sd_log_or_le500m']:.3f} & {r['mean_flips']:.0f} \\\\")
    if not rows:
        rows = [f"\\multicolumn{{7}}{{l}}{{{missing('B11')}}} \\\\"]
    base = b11.get("base_or") or {}
    cap = ("The estimator of Sect.~5 of the paper (AOI fixed effects, distance bands, logit of the terrain prediction) "
           "applied to labels with a known injected odds ratio (exploratory, amendment 10, item B11; 20 replicates per "
           "row; form A: all ground nearest a water-side name, form B: the part within 500\\,m). Realised: the mean "
           "estimated odds ratio divided by the estimate on the observed labels"
           + (f" ({orv(base['all'])} for all ground, {orv(base['le500m'])} within 500\\,m)" if base else "")
           + ". SD: standard deviation of the log odds ratio between replicates. Flips: mean number of labels changed.")
    return _tbl("tab:ten-injected", cap, "lrrrrrr", r"Form & Nominal OR & Realised, all ground & Realised, 500\,m & "
                r"SD, all & SD, 500\,m & Flips", rows)


def table_ten_dilution(v44: dict | None) -> str:
    p1 = (v44 or {}).get(P1_MEASURED_KEY, {}).get("rate")
    p0 = (v44 or {}).get(P0_MEASURED_KEY, {}).get("rate")
    if p1 is None or p0 is None:
        rows = [f"\\multicolumn{{5}}{{l}}{{{missing('v44')}}} \\\\"]
    else:
        rows = [f"{orn:g} & " + " & ".join(orv(dilution(orn, phi, p1, p0)) for phi in DIL_PHI) + " \\\\"
                for orn in DIL_ORN]
        rows.append(r"\multicolumn{5}{l}{\emph{Without misreading ($p_1=1$, $p_0=0$)}} \\")
        rows += [f"{orn:g} & " + " & ".join(orv(dilution(orn, phi, 1.0, 0.0)) for phi in DIL_PHI) + " \\\\"
                 for orn in DIL_ORN]
    cap = ("How an effect confined to the ground a name refers to appears in the odds ratio for all ground nearest the "
           "name (amendment 10, addendum 1, item 3; arithmetic, no data). $\\varphi$: share of the nearest ground that "
           "the named feature occupies; OR$_n$: odds ratio of flooding there. Rare-outcome approximation: risk ratio "
           "$1+\\varphi(\\mathrm{OR}_n-1)$ on the ground of truly water-related names, observed through names read with "
           + (f"precision $p_1={p1:.2f}$ and a water-related share $p_0={p0:.2f}$ among other names, both measured with "
              "the codes of the model that compiled the lexicon (Table~5 and Sect.~6.1 of the paper; $p_0$ from the "
              "names not read as water-side in the sample of Wikipedia derivations)." if p1 is not None
              else "the measured misreadings."))
    head = r"OR$_n$ & " + " & ".join(f"$\\varphi={phi:g}$" for phi in DIL_PHI)
    return _tbl("tab:ten-dilution", cap, "lrrrr", head, rows)


def table_ten_misclass(v66: dict | None) -> str:
    b5 = (v66 or {}).get("B5") or {}
    by = {}
    for r in b5.get("rows", []):
        by.setdefault((r["p1"], r["p0"]), {})[r["or_obs"]] = r["or_corrected"]
    rows = []
    for (p1, p0), d in sorted(by.items()):
        cell = [orv(d[k]) if d.get(k) is not None else "undefined" for k in ("point", "pairs_lo", "pairs_hi")]
        rows.append(f"{p1:.1f} & {p0:.2f} & " + " & ".join(cell) + " \\\\")
    if not rows:
        rows = [f"\\multicolumn{{5}}{{l}}{{{missing('B5')}}} \\\\"]
    cap = ("Odds ratio near water-side names corrected for misreading over a grid of rates (exploratory, amendment 10, "
           "item B5): $\\mathrm{OR}=(\\mathrm{OR}^*(1-p_0)-(1-p_1))/(p_1-\\mathrm{OR}^*p_0)$, with OR$^*$ the observed "
           "estimate within 500\\,m and the limits of its pairs-bootstrap 95\\,\\% interval; $p_1$: share of water-related "
           "names among those read as water-side; $p_0$: share among the others. Undefined: the observed value cannot "
           "arise from these rates.")
    return _tbl("tab:ten-misclass", cap, "rrrrr", r"$p_1$ & $p_0$ & Observed point & Observed lower & Observed upper",
                rows)


def table_ten_floodplain(v66: dict | None) -> str:
    b2 = (v66 or {}).get("B2") or {}

    def orow(lab: str, d: dict | None, extra: str = "") -> str:
        if not d or d.get("odds_ratio") is None:
            return f"{lab} & \\multicolumn{{4}}{{l}}{{{missing(lab)}}} \\\\"
        return (f"{lab} & {orv(d['odds_ratio'])} & {orv(d['pairs_ci95'][0])}--{orv(d['pairs_ci95'][1])} & "
                f"{int(d['n_group_flooded'])} & {extra or d['n_group_units']} \\\\")
    rows = [orow("All AOIs (floodplain reading)", b2.get("main")),
            orow("\\quad re-read through the same code, all forms", b2.get("reread_all_forms")),
            orow("\\quad re-read without bank, au, dike and polder forms", b2.get("reread_without_forms")),
            orow("Without activation EMSR097", b2.get("without_EMSR097"))]
    rows.append(r"\multicolumn{5}{l}{\emph{One country left out}} \\")
    for c, d in sorted((b2.get("leave_one_country_out") or {}).items(),
                       key=lambda kv: -kv[1].get("n_flooded_near_in_country", 0)):
        nf = d.get("n_flooded_near_in_country", "")
        rows.append(orow(f"\\quad without {c} ({nf} flooded cell{'' if nf == 1 else 's'})", d))
    rows.append(r"\multicolumn{5}{l}{\emph{Composition by matched form (flooded cells within 500\,m; top 12)}} \\")
    for r in (b2.get("composition_by_form") or [])[:12]:
        mark = "$^{x}$" if r["excluded_in_reread"] else ""
        rows.append(f"\\quad {r['category']} \\emph{{{r['form']}}}{mark} & & & {r['n_flooded_le500m']} & "
                    f"{r['n_cells_le500m']} cells \\\\")
    cap = ("The floodplain reading (riverine and wetland classes) within 500\\,m of settlements, given the terrain "
           "prediction (exploratory, amendment 10, item B2; defined after the class estimates had been seen). Intervals: "
           "pairs cluster bootstrap over location groups, 2,000 replicates. Countries are those of the nearest settlement. "
           "$^{x}$Forms removed in the re-reading. Columns 4--5: flooded cells and location groups with such cells (or "
           "cells, for the forms).")
    return _tbl("tab:ten-floodplain", cap, "p{7.6cm}llrr", r"Estimate & Odds ratio & 95\,\% interval & Flooded & "
                r"Groups", rows)


def table_ten_other_lang(v66: dict | None) -> str:
    b14 = (v66 or {}).get("B14") or {}
    rows = [r"\multicolumn{5}{l}{\emph{Countries whose first language is not among the lexicon blocks applied}} \\"]
    for c, d in sorted((b14.get("country_level") or {}).items(), key=lambda kv: -kv[1]["n_cells"]):
        rows.append(f"{c} ({d['first_language']}; read as {', '.join(d['lexicon_languages']) or 'none'}) & "
                    f"{_thousands(d['n_cells'])} & {_thousands(d['n_water_side'])} & {d['n_flooded_le500m']} & \\\\")
    rows.append(r"\multicolumn{5}{l}{\emph{Regions with a language other than the one read}} \\")
    for k, d in (b14.get("regional") or {}).items():
        rows.append(f"{k} & {_thousands(d['n_cells'])} & {_thousands(d['n_water_side'])} & {d['n_flooded_le500m']} & "
                    f"{d['n_aoi']} \\\\")
    if not b14:
        rows = [f"\\multicolumn{{5}}{{l}}{{{missing('B14')}}} \\\\"]
    cap = ("Analysed cells whose names are read in a language other than their own (exploratory, amendment 10, item "
           "B14). Country level: the first language of the GeoNames country record is not among the lexicon blocks "
           "applied to the country. Regions: Natural Earth first-level units (point in polygon of the cell centre).")
    return _tbl("tab:ten-other-lang", cap, "p{7.8cm}rrrr", r"Country or region & Cells & Water-side & "
                r"Flooded, 500\,m & AOIs", rows)


B13_LANGS = (("de", "German"), ("en", "English"), ("fr", "French"), ("pooled", "All three"))
NREAD_CLASSES = (("water", "water"), ("wetland", "wetland"), ("riverine", "riverine"), ("crossing", "crossing"),
                 ("lowland", "lowland"), ("highland", "highland"), ("other_environmental", "other environmental"),
                 ("non_environmental", "non-environmental"), ("unknown", "unknown"))


def table_ten_between(v66: dict | None) -> str:
    b13 = (v66 or {}).get("B13") or {}
    rows = []
    for key, lab in B13_LANGS:
        d = b13.get(key)
        if not d:
            rows.append(f"{lab} & \\multicolumn{{4}}{{l}}{{{missing('B13')}}} \\\\")
            continue
        cs = [f"${d['coef'][i]:+.2f}$ (${d['ci95'][i][0]:+.2f}$ to ${d['ci95'][i][1]:+.2f}$)" for i in (0, 1)]
        rows.append(f"{lab} & {d['n_aoi']} & {d['n_activations']} & {cs[0]} & {cs[1]} \\\\")
    cap = ("Whether areas with more water-side names flood more, within language regions (exploratory; amendment 10, "
           "item B13). Regression of the empirical logit of each AOI's flooded fraction, "
           "$\\log\\{(k+0.5)/(n-k+0.5)\\}$, on the share of its analysed cells whose nearest name is water-side and on "
           "the AOI mean of the terrain model's held-out logit, with activation fixed effects, for the AOIs of countries "
           "read with the German, English or French lexicon blocks; AOIs that are the only analysed AOI of their "
           "activation drop out with the fixed effects. Standard errors clustered by activation (CR1), 95\\,\\% "
           "intervals from the $t$ distribution with one degree of freedom fewer than the activations.")
    return _tbl("tab:ten-between", cap, "lrrcc", r"Language & AOIs & Activations & Share water-side & "
                r"Mean logit of M1", rows)


def table_ten_nread(v66: dict | None) -> str:
    b15 = (v66 or {}).get("B15") or {}
    bn, bc = b15.get("by_names") or {}, b15.get("by_cells") or {}
    rows = [f"{lab} & {_thousands(bn[k]['n'])} & {pct1(bn[k]['share'])} & {_thousands(bc[k]['n'])} & "
            f"{pct1(bc[k]['share'])} \\\\" for k, lab in NREAD_CLASSES if k in bn and k in bc]
    if not rows:
        rows = [f"\\multicolumn{{5}}{{l}}{{{missing('B15')}}} \\\\"]
    cap = ("Classes that the language-model reading (N-read) assigned to the nearest names of the analysed cells "
           "(exploratory; amendment 10, item B15), counted once per name and weighted by the number of analysed cells "
           "nearest each name. ``Unknown'': the model gave no class.")
    return _tbl("tab:ten-nread", cap, "lrrrr", r"Class & Names & Share & Cells & Share", rows)


def table_ten_valid(v66: dict | None) -> str:
    b6 = (v66 or {}).get("B6") or {}
    rows = []
    c1, c2 = (b6.get("coder1") or {}).get("by_country") or {}, (b6.get("coder2") or {}).get("by_country") or {}
    for cc in sorted(set(c1) | set(c2), key=lambda c: -(c1.get(c) or c2.get(c))["n_coded"]):
        a, b = c1.get(cc) or {}, c2.get(cc) or {}

        def f(d: dict, k: str) -> str:
            return "---" if d.get(k) is None else pct1(d[k])
        rows.append(f"{cc} & {a.get('n_coded', '---')} & {a.get('n_flagged', '---')} & {f(a, 'precision')} & "
                    f"{f(a, 'recall_weighted')} & {f(b, 'precision')} & {f(b, 'recall_weighted')} \\\\")
    for name in ("coder1", "coder2"):
        for form, d in ((b6.get(name) or {}).get("by_form") or {}).items():
            rows.append(f"Form \\emph{{{form}}} ({name.replace('coder', 'coder ')}) & & {d['n_flagged']} & "
                        f"{pct1(d['precision'])} & & & \\\\")
    if not rows:
        rows = [f"\\multicolumn{{7}}{{l}}{{{missing('B6')}}} \\\\"]
    cap = ("Precision and recall of the lexicon's water-side reading by country, against etymologies recorded on "
           "Wikipedia (exploratory, amendment 10, item B6): countries with at least ten coded names. Recall weights the two "
           "sampling strata by their sizes within the country, whereas Table~\\ref{tab:seven-languages} uses the design "
           "weights of the whole sample (stratum size over sample size); the two agree on precision and differ on "
           "recall where the sampling fraction differs between countries. Forms: the only matched forms with at least ten coded names read as water-side. "
           "Both coders are language models; coder 1 compiled the lexicon.")
    return _tbl("tab:ten-valid", cap, "lrrrrrr", r"Country & Coded & Read as water-side & Precision (1) & Recall (1) & "
                r"Precision (2) & Recall (2)", rows)


def table_ten_text() -> str:
    # B-12（A3〜A8 の組み直しと A12）: 新しい本文の節に合わせ、本文が言い、他の表に無い量だけを残す
    rows = [
        r"\multicolumn{2}{l}{\emph{Data and sample (Sects.~2--3)}} \\",
        r"\quad Flooded share of the analysed cells; AOI median [interquartile range] & \BaseRatePct; \AoiRateMed{} "
        r"[\AoiRateQOne, \AoiRateQThree] \\",
        r"\quad Evaluable cells neither flooded nor dry; cells partly permanent water or outside the image & "
        r"\NeitherCells; \DenPartialPct \\",
        r"\quad AOIs without a hydrography layer; analysed AOIs whose products describe no flood type, and those counted "
        r"riverine by the inland rule & \NoHydroNExcl{} (\NoHydroPct); \NUntypedAoi; \NUnspecifiedRiverine \\",
        r"\quad Settlements nearest the analysed cells: sections of populated places; historical, abandoned or "
        r"destroyed; coordinates on whole arc-minutes & \FcodeSectionPct; \FcodeHistPct; \RoundedSettlementsPct \\",
        r"\quad Share of the cells and of the label variance $\sum p(1-p)$: Asian AOIs; European AOIs & "
        r"\ShareCellsAsia, \ShareVarAsia; \ShareCellsEurope, \ShareVarEurope \\",
        r"\quad Effect calibrated to the margin called equivalent against the real margin (synthetic labels) & "
        r"\PowerEquivReal \\",
        r"\multicolumn{2}{l}{\emph{Uncertainty (Sect.~4)}} \\",
        r"\quad $R$ over thirteen assignments of location groups to folds: smallest; largest & \FoldRMin; \FoldRMax \\",
        r"\quad Correlation of AOI-mean residuals, AOIs of one activation within 50\,km (Fisher 95\,\%); between "
        r"activations & \ResActCorrNear{} (\ResActCorrNearLo{} to \ResActCorrNearHi; \ResActNPairsNear{} pairs); "
        r"\ResActCorrOther \\",
        r"\multicolumn{2}{l}{\emph{What assignment through settlements can carry (Sect.~5)}} \\",
        r"\quad Share of cells within 500\,m of a settlement with a water-side name & \CeilShareFiveHundred \\",
        r"\quad Largest possible $R$ from perfect prediction near water-side names within 1\,km & \CeilROneKm \\",
        r"\quad Recurrence assigned through settlements on the AOIs with other events, on the pooled scale & "
        r"\RecurPooledR \\",
        r"\multicolumn{2}{l}{\emph{Readings and labels near settlements (Sect.~6)}} \\",
        r"\quad Cells whose nearest name is read as riverside by a form such as French \emph{-sur-} & \RiversidePct \\",
        r"\quad Settlements with a namesake in Wikidata; of them named after a river; river namesakes read as "
        r"water-side by the language model; by the lexicon & \WDNamesakes; \WDRiverNamesakes; \WDNreadHits; "
        r"\WDLexHits \\",
        r"\quad Language-model reading: cells water-side; in Africa; cells answered ``unknown'' & \NReadWetPct; "
        r"\NReadWetAF; \NReadUnknownPct \\",
        r"\quad Planned estimate on the GFM labels: flooded cells within 500\,m of water-side names; location groups "
        r"holding them & \GfmCappedFlooded; \GfmCappedUnits \\",
        r"\quad Cells within 500\,m that GFM can assess: near highland names; near water-side names & "
        r"\HighAssessPct; \WsAssessPct \\",
        r"\quad Joint Wald test, class coefficients in the three distance bands: $p$ & \WaldClassP \\",
    ]
    cap = ("Quantities stated in the paper and not given in another table of the paper or of this Supplement "
           "(exploratory unless pre-specified). Intervals: 95\\,\\% unless marked. The values are typeset from the "
           "same generated file as the paper.")
    return _tbl("tab:ten-text", cap, "p{11.0cm}l", "Quantity & Value", rows)


def supp_tables_b10(v66: dict | None, v67: dict | None, v44: dict | None) -> None:
    for name, tex in (("ten_hand", table_ten_hand(v66)), ("ten_robust", table_ten_robust(v66, v67)),
                      ("ten_islands", table_ten_islands(v66)), ("ten_between", table_ten_between(v66)),
                      ("ten_nread", table_ten_nread(v66)),
                      ("ten_splits", table_ten_splits(v66)), ("ten_injected", table_ten_injected(v66)),
                      ("ten_dilution", table_ten_dilution(v44)), ("ten_misclass", table_ten_misclass(v66)),
                      ("ten_floodplain", table_ten_floodplain(v66)), ("ten_other_lang", table_ten_other_lang(v66)),
                      ("ten_valid", table_ten_valid(v66)), ("ten_text", table_ten_text())):
        (HERE / f"app_table_{name}.tex").write_text(tex, encoding="utf-8")


# %% B-11（初見の査読 第 7 回）: v68 の数と補足の表
def _or3p(prefix: str, o: dict | None, what: str) -> dict[str, str]:
    """OR と pairs の 95% 区間（{prefix}, Lo, Hi）。"""
    if not isinstance(o, dict) or o.get("odds_ratio") is None or not o.get("pairs_ci95"):
        return {f"{prefix}{s}": missing(what) for s in ("", "Lo", "Hi")}
    return {prefix: orv(o["odds_ratio"]), f"{prefix}Lo": orv(o["pairs_ci95"][0]), f"{prefix}Hi": orv(o["pairs_ci95"][1])}


def b11_macros(v68: dict | None) -> dict[str, str]:  # noqa: C901, PLR0912, PLR0915
    out: dict[str, str] = {}
    v = v68 or {}
    # B1: 活性化の分割の変化
    b1 = v.get("B1") or {}
    if b1.get("R_frozen") is not None:
        top = (b1.get("largest_unit_changes") or [{}])[0]
        out |= {"ActRFrozen": rv(b1["R_frozen"]), "ActRCorr": rv(b1["R_corrected"]),
                "ActChangeR": rv(b1["change_R"]), "ActChangeNorway": f"\\ensuremath{{{b1['change_R_norwegian_units']:.5f}}}",
                "ActChangeTop": rv(top.get("change_R", float("nan"))), "ActChangeTopCC": top.get("majority_cc", "?"),
                "ActChangeTopCountry": {"PK": "Pakistan", "BD": "Bangladesh", "IR": "Iran", "AU": "Australia",
                                        "DE": "Germany"}.get(top.get("majority_cc", ""), top.get("majority_cc", "?")),
                "ActNorwayUnits": _word(b1["n_norwegian_units"])}
    else:
        out |= {k: missing("v68 B1") for k in ("ActRFrozen", "ActRCorr", "ActChangeR", "ActChangeNorway", "ActChangeTop",
                                              "ActChangeTopCC", "ActNorwayUnits")}
    # B2: 上限なしの集落の近くの OR
    b2 = v.get("B2") or {}
    for key, src in (("NearUncap", "ems_all"), ("GfmUncap", "gfm"), ("EmsSame", "ems_same_cells"),
                     ("NotSOneGfm", "gfm_not_s1"), ("NotSOneEms", "ems_not_s1")):
        o = b2.get(src)
        out |= _or3p(f"{key}OR", o, f"v68 B2 {src}")
        out[f"{key}Flooded"] = _thousands(o["n_flooded_near"]) if isinstance(o, dict) and "n_flooded_near" in o \
            else missing(f"v68 B2 {src}")
        out[f"{key}Units"] = str(o["n_units_flooded_near"]) if isinstance(o, dict) and "n_units_flooded_near" in o \
            else missing(f"v68 B2 {src}")
    a = b2.get("ems_all") or {}
    out |= ({"NearUncapMdeUp": orv(a["mde80_up"]), "NearUncapMdeDown": orv(a["mde80_down"])} if a.get("mde80_up")
            else {k: missing("v68 B2") for k in ("NearUncapMdeUp", "NearUncapMdeDown")})
    out["NotSOneNAoi"] = str(b2["n_aoi_not_s1"]) if b2.get("n_aoi_not_s1") is not None else missing("v68 B2")
    out["GfmUncapNAoi"] = str(b2["n_aoi_gfm"]) if b2.get("n_aoi_gfm") is not None else missing("v68 B2")
    out["WetsideChangedUncapPct"] = (pct1(b2["share_wetside_changed_vs_frozen"])
                                     if b2.get("share_wetside_changed_vs_frozen") is not None else missing("v68 B2"))
    # B3: 推定量の尺度
    b3 = v.get("B3") or {}
    mg = (b3.get("margin") or {}).get("A")
    out["MarginEstA"] = mg["text"] if mg else missing("v68 B3")
    bd = b3.get("bounds") or {}
    for key, iv in (("BoundEstA", "fixed"), ("BoundFoldEstA", "fold"), ("BoundRefEstA", "refit")):
        e = ((bd.get(iv) or {}).get("b7_recomputed") or {}).get("A")
        out[key] = e["text"] if e else missing("v68 B3")
    e = ((bd.get("fixed") or {}).get("b7_recomputed") or {}).get("B")
    out["BoundEstBText"] = _gt(e["text"]) if e else missing("v68 B3")
    out["BoundACEstA"] = b3["bound_ac"]["text"] if b3.get("bound_ac") else missing("v68 B3 AC")
    cur = (b3.get("curves") or {}).get("A") or {}
    out |= ({"InjEstAtTwo": orv(cur["estimator"][cur["nominal"].index(2.0)]),
             "InjEstAtThree": orv(cur["estimator"][-1])} if 2.0 in (cur.get("nominal") or [])
            else {k: missing("v68 B3") for k in ("InjEstAtTwo", "InjEstAtThree")})
    # B4: GFM の評価域
    b4 = (v.get("B4") or {}).get("gfm_assessable_cells") or {}
    if b4.get("yardstick") is not None:
        out |= {"GfmDomYard": sci(b4["yardstick"], 2), "GfmDomYardLo": sci(b4["ci95_yardstick"][0], 2),
                "GfmDomYardHi": sci(b4["ci95_yardstick"][1], 2), "GfmDomR": rv(b4["R"]),
                "GfmDomRLoNinety": rv(b4["ci90_R"][0]), "GfmDomRHiNinety": rv(b4["ci90_R"][1]),
                "GfmDomNCells": _thousands(b4["n_cells"])}
    else:
        out |= {k: missing("v68 B4") for k in ("GfmDomYard", "GfmDomYardLo", "GfmDomYardHi", "GfmDomR",
                                               "GfmDomRLoNinety", "GfmDomRHiNinety", "GfmDomNCells")}
    sc = (v.get("B4") or {}).get("gfm_scenes_per_aoi") or {}
    out |= ({"GfmScenesMed": f"{sc['median']:g}", "GfmScenesMax": str(sc["max"]),
             "GfmScenesMultiPct": pct1(sc["share_more_than_one"])} if sc
            else {k: missing("v68 B4") for k in ("GfmScenesMed", "GfmScenesMax", "GfmScenesMultiPct")})
    # B5: 浸水の割合の分母
    b5 = v.get("B5") or {}
    cc_ = b5.get("counts_capped") or {}
    if cc_:
        out |= {"DenPartialPct": pct1(cc_["all"]["share"]), "DenPartialWetPct": pct1(cc_["le500m_wetside"]["share"]),
                "DenPartialOtherPct": pct1(cc_["le500m_other"]["share"]),
                "DenChanged": _thousands((b5.get("relabel_capped") or {}).get("n_changed_0_to_1", 0))}
    else:
        out |= {k: missing("v68 B5") for k in ("DenPartialPct", "DenPartialWetPct", "DenPartialOtherPct", "DenChanged")}
    out |= _or3p("DenOR", b5.get("or_capped"), "v68 B5")
    rr = b5.get("R_rescored") or {}
    out |= ({"DenR": rv(rr["R"]), "DenRLoNinety": rv(rr["ci90_R"][0]), "DenRHiNinety": rv(rr["ci90_R"][1])}
            if rr.get("R") is not None else {k: missing("v68 B5 R") for k in ("DenR", "DenRLoNinety", "DenRHiNinety")})
    out |= _or3p("DenORUncap", (v.get("B5u") or {}).get("or_uncapped"), "v68 B5u")
    # B6・B7: 区間
    b6 = v.get("B6") or {}
    out |= ({"NBLoNinety": rv(b6["ci90"][0]), "NBHiNinety": rv(b6["ci90"][1]), "NBMean": rv(b6["R_mean"])}
            if b6.get("ci90") else {k: missing("v68 B6") for k in ("NBLoNinety", "NBHiNinety", "NBMean")})
    b7 = (v.get("B7") or {}).get("two_way") or {}
    out |= ({"CombLoNinety": rv(b7["ci90_R_refit_ratio"][0]), "CombHiNinety": rv(b7["ci90_R_refit_ratio"][1]),
             "CombWithin": "lies within" if b7["within_margin"] else "does not lie within"}
            if b7.get("ci90_R_refit_ratio") else {k: missing("v68 B7") for k in ("CombLoNinety", "CombHiNinety",
                                                                                  "CombWithin")})
    # B8: 目録
    b8 = v.get("B8") or {}
    fl, si = b8.get("flood") or {}, b8.get("site") or {}
    ns = b8.get("near_settlement_or") or {}
    if fl:
        out |= {"CensusN": str(fl["n"]), "CensusExcl": str(fl["n_excluding_null"]),
                "CensusExclPct": pct1(fl["share_excluding_null"]), "CensusSiteN": str(si.get("n", 0)),
                "CensusSiteExcl": str(si.get("n_excluding_null", 0)), "CensusNearN": str(ns.get("n", 0)),
                "CensusNearExcl": str(ns.get("n_excluding_one", 0)),
                "CensusNearBHFive": _word(ns.get("n_bh_q_below_0_05", 0)),
                "CensusNearBHTen": _word(ns.get("n_bh_q_below_0_10", 0))}
    else:
        out |= {k: missing("v68 B8") for k in ("CensusN", "CensusExcl", "CensusExclPct", "CensusSiteN", "CensusSiteExcl",
                                               "CensusNearN", "CensusNearExcl", "CensusNearBHFive", "CensusNearBHTen")}
    # B9: 地物の符号と丸めた座標
    b9 = v.get("B9") or {}
    bs = b9.get("by_fcode_settlements") or {}
    tot = sum(x["settlements"] for x in bs.values())
    if tot:
        out |= {"FcodeSectionPct": pct1(bs.get("PPLX", {}).get("settlements", 0) / tot),
                "FcodeHistPct": pct1(sum(bs.get(c, {}).get("settlements", 0) for c in ("PPLH", "PPLQ", "PPLW")) / tot),
                "RoundedSettlementsPct": pct1(b9["rounded_share_settlements"])}
    else:
        out |= {k: missing("v68 B9") for k in ("FcodeSectionPct", "FcodeHistPct", "RoundedSettlementsPct")}
    # B10: パイロット
    b10 = {k: x for k, x in (v.get("B10") or {}).items() if isinstance(x, dict)}
    if b10:
        rn = [x["R_names"] for x in b10.values()]
        rp = [x["R_placebo"] for x in b10.values()]
        out |= {"PilotRMin": rv(min(rn)), "PilotRMax": rv(max(rn)), "PilotRpMin": rv(min(rp)), "PilotRpMax": rv(max(rp))}
    else:
        out |= {k: missing("v68 B10") for k in ("PilotRMin", "PilotRMax", "PilotRpMin", "PilotRpMax")}
    # B11: 照合の標本
    b11 = v.get("B11") or {}
    out |= ({"AuditNNew": str(b11["n_names"]), "AuditStrata": str(b11["n_strata_used"]),
             "AuditNCountries": _word(len(b11["countries"])),
             "AuditFloodCoverPct": pct1(b11["share_flooded_le500m_covered"]),
             "AuditCellCoverPct": pct1(b11["share_wetside_cells_covered"])} if b11.get("n_names")
            else {k: missing("v68 B11") for k in ("AuditNNew", "AuditStrata", "AuditNCountries", "AuditFloodCoverPct",
                                                  "AuditCellCoverPct")})
    # C1: 国外しの組み立て直した置換
    c1 = v.get("C1") if isinstance(v.get("C1"), dict) else {}

    def _ordinal(k: int) -> str:
        suf = "th" if 10 <= k % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")  # noqa: PLR2004
        return f"{k}{suf}"

    out |= ({"CtryPermPct": _ordinal(round(100 * c1["percentile"])), "CtryPermP": f"{c1['p_one_sided']:.2f}",
             "CtryPermN": str(c1["n_perms"])} if c1.get("percentile") is not None
            else {k: missing("v68 C1") for k in ("CtryPermPct", "CtryPermP", "CtryPermN")})
    return out


def table_eleven_near(v68: dict | None) -> str:
    v = v68 or {}
    b2, b5, b5u = v.get("B2") or {}, v.get("B5") or {}, v.get("B5u") or {}
    rows = []

    def row(lab: str, o: dict | None, what: str) -> str:
        if not isinstance(o, dict) or o.get("odds_ratio") is None:
            return f"{lab} & \\multicolumn{{5}}{{l}}{{{missing(what)}}} \\\\"
        mde = (f"{orv(o['mde80_down'])}, {orv(o['mde80_up'])}" if o.get("mde80_up") else "---")
        return (f"{lab} & {orv(o['odds_ratio'])} & {orv(o['pairs_ci95'][0])}--{orv(o['pairs_ci95'][1])} & "
                f"{_thousands(o.get('n_flooded_near', o.get('n_group_flooded', 0)))} & "
                f"{o.get('n_units_flooded_near', o.get('n_group_units', '---'))} & {mde} \\\\")
    rows.append(r"\multicolumn{6}{l}{\emph{Every evaluable cell (no cap), corrected-code readings}} \\")
    rows.append(row("EMS labels, 500 AOIs", b2.get("ems_all"), "v68 B2"))
    rows.append(row("GFM labels, cells that GFM can assess", b2.get("gfm"), "v68 B2"))
    rows.append(row("EMS labels on the same cells", b2.get("ems_same_cells"), "v68 B2"))
    rows.append(row(f"GFM labels, AOIs not mapped from Sentinel-1 ({b2.get('n_aoi_not_s1', '?')})",
                    b2.get("gfm_not_s1"), "v68 B2"))
    rows.append(row("EMS labels, the same AOIs and cells", b2.get("ems_not_s1"), "v68 B2"))
    rows.append(row("EMS labels, flooded fraction over observed non-permanent water", b5u.get("or_uncapped"),
                    "v68 B5u"))
    rows.append(r"\multicolumn{6}{l}{\emph{Capped table (10,000 cells per location group)}} \\")
    mn = dict(((_json("v57_b8_local.json") or {}).get("B6") or {}).get("main") or {})
    ca = (((_json("v50_gfm_diagnostics.json") or {}).get("contributing") or {}).get("ems_all_aoi")) or {}
    if mn and ca.get("n_flooded_group") == mn.get("n_group_flooded"):  # B-12 の A12: 浸水セルを持つグループの数に揃える
        mn["n_units_flooded_near"] = ca["n_units_with_flooded"]
    rows.append(row("EMS labels", mn or None, "v57 main"))
    rows.append(row("EMS labels, flooded fraction over observed non-permanent water", b5.get("or_capped"), "v68 B5"))
    cap = ("Odds ratio of flooding within 500\\,m of settlements with water-side names against other names of the same "
           "AOI, given the held-out prediction of the baseline model (exploratory; amendment 11, items B2 and B5): AOI "
           "fixed effects, distance bands and the logit of the prediction, with a 95\\,\\% pairs cluster bootstrap over "
           "location groups (2,000 replicates). Flooded: flooded cells within 500\\,m of water-side names and the "
           "location groups that hold them. MDE: the smallest odds ratios below and above one detectable with 80\\,\\% "
           "power (two-sided 5\\,\\%). The uncapped sample conditions on the predictions of the uncapped sensitivity run "
           "(the baseline model does not use names). The GFM labels are not independent of EMS where the EMS map was "
           "drawn from Sentinel-1.")
    return _tbl("tab:eleven-near", cap, "p{6.6cm}llrrl", r"Sample and labels & Odds ratio & 95\,\% interval & "
                r"Flooded & Groups & MDE", rows)


def table_eleven_scale(v68: dict | None) -> str:
    b3 = (v68 or {}).get("B3") or {}
    rows = []
    mg = (b3.get("margin") or {}).get("A")
    rows.append("Margin ($R=0.05$), all ground nearest a water-side name & "
                + (f"{orv(mg['nominal'])} & {mg['text']} \\\\" if mg else f"\\multicolumn{{2}}{{l}}{{{missing('B3')}}} \\\\"))
    for iv, lab in (("fixed", "fitted models held fixed"), ("fold", "with the spread over fold assignments"),
                    ("refit", "with the second stage refitted")):
        for f, flab in (("A", "all ground"), ("B", "within 500\\,m")):
            e = (((b3.get("bounds") or {}).get(iv) or {}).get("b7_recomputed") or {}).get(f)
            rows.append(f"Upper bound, {flab}, {lab} & "
                        + (f"{orv(e['nominal'])} & {_gt(e['text'])} \\\\" if e
                           else f"\\multicolumn{{2}}{{l}}{{{missing('B3')}}} \\\\"))
    for key, lab in (("bound_ac", "flips correlated within each settlement's ground"),
                     ("bound_het", "effects that differ between settlements (mean odds ratio)")):
        e = b3.get(key)
        rows.append(f"Upper bound, all ground, {lab} & "
                    + (f"{orv(e['nominal'])} & {_gt(e['text'])} \\\\" if e
                       else f"\\multicolumn{{2}}{{l}}{{{missing('B3 ' + key)}}} \\\\"))
    cur = b3.get("curves") or {}
    note = ""
    if cur.get("A"):
        pairs = ", ".join(f"{a:g} $\\to$ {b:.2f}" for a, b in zip(cur["A"]["nominal"][1:], cur["A"]["estimator"][1:],
                                                                    strict=True))
        note = f" Calibration curve for all ground (nominal $\\to$ estimator): {pairs}."
    cap = ("The margin and the bounds of Sect.~4 as injected (nominal) odds ratios and on the scale of the estimator of "
           "Sect.~5 (exploratory; amendment 11, item B3). The estimator was applied to labels with known injected odds "
           "ratios (Table~\\STabTenInjected); its mean odds ratio, relative to the uninjected labels, is interpolated on "
           "the log scale. Beyond the largest injected odds ratio (3) the estimator value is given as a lower limit. "
           "The no-effect reference is the mean increment of reassigned names with rebuilt composition." + note)
    return _tbl("tab:eleven-scale", cap, "p{9.2cm}ll", r"Quantity & Injected odds ratio & Estimator odds ratio", rows)


def table_eleven_checks(v68: dict | None) -> str:
    v = v68 or {}
    rows = []
    b1 = v.get("B1") or {}
    rows.append(r"\multicolumn{3}{l}{\emph{Activation split: frozen and corrected code (item B1)}} \\")
    if b1.get("R_frozen") is not None:
        rows.append(f"$R$, frozen code; corrected code & {rv(b1['R_frozen'])}; {rv(b1['R_corrected'])} & "
                    "M0 and M1 identical in the two runs \\\\")
        rows.append(f"Change of $R$; from location groups with Norwegian cells & {rv(b1['change_R'])}; "
                    f"\\ensuremath{{{b1['change_R_norwegian_units']:.5f}}} & {b1['n_norwegian_units']} groups \\\\")
        for u in (b1.get("largest_unit_changes") or [])[:3]:
            rows.append(f"\\quad unit {u['unit']} of the activation split ({u['majority_cc']}, "
                        f"{_thousands(u['n_cells'])} cells) & "
                        f"{rv(u['change_R'])} & \\\\")
    else:
        rows.append(f"\\multicolumn{{3}}{{l}}{{{missing('v68 B1')}}} \\\\")
    rows.append(r"\multicolumn{3}{l}{\emph{Intervals for $R$ on the primary split (items B6 and B7)}} \\")
    b6 = v.get("B6") or {}
    rows.append("Refitting over 13 fold assignments, Nadeau--Bengio, 90\\,\\% & "
                + (f"{rv(b6['ci90'][0])} to {rv(b6['ci90'][1])} & mean {rv(b6['R_mean'])} \\\\" if b6.get("ci90")
                   else f"\\multicolumn{{2}}{{l}}{{{missing('v68 B6')}}} \\\\"))
    for sch, lab in (("two_way", "activations $\\times$ countries"), ("groups", "location groups")):
        b7 = (v.get("B7") or {}).get(sch) or {}
        rows.append(f"Second stage refitted, ratio bootstrap, {lab}, 90\\,\\% & "
                    + (f"{rv(b7['ci90_R_refit_ratio'][0])} to {rv(b7['ci90_R_refit_ratio'][1])} & "
                       f"{'within' if b7['within_margin'] else 'not within'} $\\pm0.05$ \\\\" if b7.get("ci90_R_refit_ratio")
                       else f"\\multicolumn{{2}}{{l}}{{{missing('v68 B7')}}} \\\\"))
    c1 = v.get("C1") if isinstance(v.get("C1"), dict) else {}
    rows.append(r"\multicolumn{3}{l}{\emph{Country split: reassignments with the composition recomputed (item C1)}} \\")
    if c1.get("percentile") is not None:
        rows.append(f"$R$ of the real names (corrected code); mean of {c1['n_perms']} reassignments (2.5--97.5\\,\\%) & "
                    f"{rv(c1['R_real'])}; {rv(c1['R_perm_mean'])} ({rv(c1['R_perm_q025'])} to "
                    f"{rv(c1['R_perm_q975'])}) & \\\\")
        rows.append(f"Percentile of the real names; one-sided $p$ & {_ordinal(round(100 * c1['percentile']))}; "
                    f"{c1['p_one_sided']:.2f} & \\\\")
    else:
        rows.append(f"\\multicolumn{{3}}{{l}}{{{missing('v68 C1')}}} \\\\")
    rows.append(r"\multicolumn{3}{l}{\emph{Domain of the GFM labels: EMS-trained predictions on EMS labels (item B4)}} \\")
    for key, lab in (("gfm_aois_all_cells", "all cells of the GFM AOIs"),
                     ("gfm_assessable_cells", "cells that GFM can assess")):
        b4 = (v.get("B4") or {}).get(key) or {}
        rows.append(f"Baseline gain and $R$, {lab} & "
                    + (f"{sci(b4['yardstick'], 2)} ({sci(b4['ci95_yardstick'][0], 2)} to "
                       f"{sci(b4['ci95_yardstick'][1], 2)}); {rv(b4['R'])} & {_thousands(b4['n_cells'])} cells \\\\"
                       if b4.get("yardstick") is not None else f"\\multicolumn{{2}}{{l}}{{{missing('v68 B4')}}} \\\\"))
    b5 = v.get("B5") or {}
    rows.append(r"\multicolumn{3}{l}{\emph{Denominator of the flooded fraction (item B5)}} \\")
    cc_ = b5.get("counts_capped") or {}
    if cc_:
        rows.append(f"Cells partly unobserved or partly permanent water: all; within 500\\,m of water-side; of other "
                    f"names & {pct1(cc_['all']['share'])}; {pct1(cc_['le500m_wetside']['share'])}; "
                    f"{pct1(cc_['le500m_other']['share'])} & \\\\")
        rl = b5.get("relabel_capped") or {}
        rows.append(f"Relabelled over the observed non-permanent-water area: changed to flooded; left out & "
                    f"{_thousands(rl.get('n_changed_0_to_1', 0))}; {rl.get('n_dropped', 0)} & \\\\")
        rr = b5.get("R_rescored") or {}
        if rr.get("R") is not None:
            rows.append(f"$R$ rescored on these labels (90\\,\\%) & {rv(rr['R'])} ({rv(rr['ci90_R'][0])} to "
                        f"{rv(rr['ci90_R'][1])}) & \\\\")
    b10 = {k: x for k, x in (v.get("B10") or {}).items() if isinstance(x, dict)}
    rows.append(r"\multicolumn{3}{l}{\emph{Pilot simulations before the freeze, no injected effect (item B10)}} \\")
    for rep, x in sorted(b10.items()):
        rows.append(f"Pilot 3, {rep}: nested increment $R$ of names; of reassigned names & {rv(x['R_names'])}; "
                    f"{rv(x['R_placebo'])} & \\\\")
    if not b10:
        rows.append(f"\\multicolumn{{3}}{{l}}{{{missing('v68 B10')}}} \\\\")
    cap = ("Further checks of amendment 11 (exploratory). B1: the change of $R$ on the activation split between the "
           "frozen and the corrected code, by location group. B6: a 90\\,\\% interval from the corrected variance of "
           "Nadeau and Bengio (2003) over the 13 fold assignments (it reflects refitting, not the sampling of test "
           "cells). B7: 500 replicates that refit the second stage and recompute the baseline gain under each "
           "resampling. B4: the predictions of the corrected run, trained on EMS labels, scored on EMS labels. B5: "
           "labels with the flooded fraction measured over the observed area that is not permanent water. B10: the "
           "pilots, with synthetic labels calibrated to the real flooded fraction of each AOI, in which names and "
           "reassigned names improved the nested comparison alike. C1: on the country split, the names reassigned "
           "within each AOI 100 times with the regional composition recomputed from each reassignment (corrected "
           "code).")
    return _tbl("tab:eleven-checks", cap, "p{9.0cm}p{3.4cm}l", r"Quantity & Value & Note", rows)


def table_eleven_census(v68: dict | None, v69: dict | None = None) -> str:
    # B-12 の追記 1: 集落の近くの族の印を直した目録（v69 の B8c）があればそれを使う。BY の q も
    b8 = (v69 or {}).get("B8c") or (v68 or {}).get("B8") or {}
    rows = []
    for fam, lab in (("flood", "Name effects on flood labels or predictions"),
                     ("site", "Name effects on terrain or water descriptors (no flood labels)")):
        x = b8.get(fam) or {}
        if not x:
            rows.append(f"{lab} & \\multicolumn{{3}}{{l}}{{{missing('v68 B8')}}} \\\\")
            continue
        rows.append(f"{lab} & {x['n']} & {x['n_excluding_null']} & {pct1(x['share_excluding_null'])} \\\\")
        for k, g_ in sorted((x.get("by_kind") or {}).items()):
            rows.append(f"\\quad {k.replace('_', ' ')} & {g_['n']} & {g_['excl']} & "
                        f"{pct1(g_['excl'] / g_['n']) if g_['n'] else '---'} \\\\")
    ns = b8.get("near_settlement_or") or {}
    if ns:
        rows.append(f"Odds ratios within 500\\,m of settlements & {ns['n']} & {ns['n_excluding_one']} & "
                    f"{pct1(ns['n_excluding_one'] / ns['n'])} \\\\")
        rows.append(f"\\quad Benjamini--Hochberg $q<0.05$; $q<0.10$ & & {ns['n_bh_q_below_0_05']}; "
                    f"{ns['n_bh_q_below_0_10']} & \\\\")
        if ns.get("n_by_q_below_0_05") is not None:
            rows.append(f"\\quad Benjamini--Yekutieli $q<0.05$; $q<0.10$ & & {ns['n_by_q_below_0_05']}; "
                        f"{ns['n_by_q_below_0_10']} & \\\\")
    fixed = "B8c" in (v69 or {})
    note = (" The odds ratios within 500\\,m of settlements are those of every estimate whose key or record names a "
            "near-settlement quantity; the first version of this census flagged them by the key alone and missed "
            f"{ns.get('n_flagged_by_keys_only', '')} of them (an implementation error corrected in addendum~1 of "
            "amendment 12). The $q$ values are approximate: they are computed from $p$ values recovered from asymmetric "
            "intervals, and the estimates are not independent. The machine-readable list is v69\\_estimates.csv in the "
            "public repository.") if fixed and ns else " The machine-readable list is v68\\_estimates.csv in the public repository."
    cap = ("Census of the exploratory interval estimates of amendments 2 to 10 (amendment 11, item B8): every estimate "
           "with a 95\\,\\% interval in the verification outputs v36 to v67, collected by key rules and counted once "
           "where the same value appears in several files. Distance-band and logit coefficients, near-minus-far "
           "contrasts, intervals that are not finite or inverted, and fits that separate are excluded, and so are the "
           "controls (records of surface water, recurrence, oracles and hazard maps). Excludes: the interval excludes the "
           "no-effect value (one for odds ratios, zero otherwise). The estimates are strongly dependent (many are "
           "variants of the same quantity). The first pass counted the excluded coefficients as name effects; the rules "
           "were refined before any count was used." + note)
    return _tbl("tab:eleven-census", cap, "p{8.6cm}rrr", r"Family & Estimates & Excludes & Share", rows)


def table_eleven_gazetteer(v68: dict | None) -> str:
    v = v68 or {}
    b9 = v.get("B9") or {}
    rows = []
    bs = b9.get("by_fcode_settlements") or {}
    bc = b9.get("by_fcode_cells") or {}
    tot = sum(x["settlements"] for x in bs.values()) or 1
    for code in sorted(bs, key=lambda c: -bs[c]["settlements"]):
        x, c = bs[code], bc.get(code, {})
        rows.append(f"{code} & {_thousands(x['settlements'])} & {pct1(x['settlements'] / tot)} & {x['wetside']} & "
                    f"{_thousands(c.get('cells', 0))} \\\\")
    if b9.get("rounded_share_settlements") is not None:
        rows.append(f"\\multicolumn{{5}}{{l}}{{Settlements whose coordinates are whole arc-minutes: "
                    f"{pct1(b9['rounded_share_settlements'])}}} \\\\")
    if not rows:
        rows = [f"\\multicolumn{{5}}{{l}}{{{missing('v68 B9')}}} \\\\"]
    cap = ("The nearest settlements of the analysed cells by GeoNames feature code (exploratory; amendment 11, item B9): "
           "settlements, their share, those read as water-side and the cells nearest them. PPL: populated place; "
           "PPLA--PPLA4: seats of administrative divisions; PPLX: section of a populated place; PPLH: historical; PPLQ: "
           "abandoned; PPLW: destroyed; PPLL: locality; PPLF: farm village; PPLR: religious.")
    return _tbl("tab:eleven-gazetteer", cap, "lrrrr", r"Feature code & Settlements & Share & Water-side & Cells", rows)


def supp_tables_b11(v68: dict | None, v69: dict | None = None) -> None:
    for name, tex in (("eleven_near", table_eleven_near(v68)), ("eleven_scale", table_eleven_scale(v68)),
                      ("eleven_checks", table_eleven_checks(v68)), ("eleven_census", table_eleven_census(v68, v69)),
                      ("eleven_gazetteer", table_eleven_gazetteer(v68))):
        (HERE / f"app_table_{name}.tex").write_text(tex, encoding="utf-8")



# %% B-11（初見の査読 第 7 回）: 本文の組み直しで使う量・表・図
_CC_NAMES = {"PK": "Pakistan", "BD": "Bangladesh", "AU": "Australia", "DE": "Germany", "PT": "Portugal",
             "GB": "the United Kingdom", "IT": "Italy", "FR": "France", "ES": "Spain", "IR": "Iran", "IN": "India",
             "RO": "Romania", "PL": "Poland", "SI": "Slovenia", "SE": "Sweden", "NO": "Norway", "IE": "Ireland",
             "GR": "Greece", "HR": "Croatia", "AT": "Austria", "CZ": "Czechia", "BG": "Bulgaria", "HU": "Hungary",
             "SK": "Slovakia", "LK": "Sri Lanka", "MM": "Myanmar", "FI": "Finland", "BE": "Belgium", "NL": "the Netherlands",
             "KH": "Cambodia", "VN": "Vietnam", "MX": "Mexico", "US": "the United States", "LA": "Laos", "TH": "Thailand",
             "PH": "the Philippines", "MG": "Madagascar", "MZ": "Mozambique"}


def _country_list(codes: list[str]) -> str:
    """国の符号を英語の名に直し、最後を and でつなぐ（本文で ISO の符号を出さない）。"""
    names = [_CC_NAMES.get(c) or missing(f"country name {c}") for c in codes]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def b11_text_macros(v36: dict | None, v54: dict | None, v57: dict | None, v57b: dict | None) -> dict[str, str]:
    """B-11 の A（書き直し）で本文が引く量。既存の成果物から読むだけで、新しい計算はしない。"""
    out: dict[str, str] = {}
    met = g._json(DATA / "runs" / "world_v2" / "metrics_primary.json") or {}
    br = met.get("brier") or {}
    # 余白（ΔBS）と、基礎のモデルを較正し直すだけの改善（M1 → M1s）を並べる（査読者 2 の M2(d)）
    out["RecalGainBS"] = sci(br["M1"] - br["M1s"], 2) if {"M1", "M1s"} <= br.keys() else missing("brier M1s")
    out["PrimBSMOnec"] = f"{br['M1c']:.5f}" if "M1c" in br else missing("brier M1c")
    # H5 の範囲: プロトコルの層（意味・地形・市街地）に大陸の層を足す（編集者の 4 節の 1、査読者 2・3）
    h5 = [v["R_within"] for v in (met.get("H5_strata") or {}).values()
          if isinstance(v, dict) and v.get("R_within") is not None]
    h5 += [v["R"] for v in (met.get("by_continent") or {}).values() if isinstance(v, dict) and v.get("R") is not None]
    out |= ({"HFiveRMinAll": rv(min(h5)), "HFiveRMaxAll": rv(max(h5))} if h5
            else {k: missing("H5 strata") for k in ("HFiveRMinAll", "HFiveRMaxAll")})
    # H4 の全体の増分を活性化 × 国の二方向 bootstrap で（表 S7 の Y0 の行）
    tw = ((((v36 or {}).get("splits") or {}).get("llmcut") or {}).get("two_way") or {}).get("H4_total") or {}
    y0 = tw.get("Y0") or {}
    out |= ({"RHFourTotalTw": rv(y0["R"]), "RHFourTotalTwLo": rv(y0["ci90_R"][0]),
             "RHFourTotalTwHi": rv(y0["ci90_R"][1])} if y0.get("ci90_R")
            else {k: missing("v36 H4 total two-way") for k in ("RHFourTotalTw", "RHFourTotalTwLo", "RHFourTotalTwHi")})
    # 包含: 浸水セルが 20 未満の AOI のうち、浸水が他の種別だったもの
    wet = (((v54 or {}).get("inclusion_flow") or {}).get("excluded_by_reason") or {}).get("wet<20") or {}
    out["NExclWetOther"] = str(wet["type"].get("other", 0)) if wet.get("type") else missing("inclusion flow")
    # 言語モデルのエージェントによる批評の回数（記録を見られた回と、原稿だけの回）
    rev = REPO / "docs" / "research" / "reviews"
    blind = sorted(rev.glob("*_paper2_blind*_editor.md"))
    access = sorted(rev.glob("*_paper2_opus_referee*.md"))
    out |= {"NBlindRounds": _word(len(blind)), "NCritiqueRounds": _word(len(blind) + len(access)),
            "NAccessRounds": _word(len(access))}
    if v57b:
        out |= {"CodeTwoAgreeStatus": _pct0(v57b["status_agreement"]), "CodeTwoNBoth": str(v57b["n_both_coded"])}
    else:
        out |= {k: missing("v57b") for k in ("CodeTwoAgreeStatus", "CodeTwoNBoth")}
    m1 = (((v57 or {}).get("B8") or {}).get("primary") or {}).get("M1") or {}
    out |= ({"CalibPrimMOne": _f2(m1["slope"]), "CalibPrimMOneInt": f"\\ensuremath{{{m1['intercept']:.2f}}}"} if m1
            else {k: missing("B8 primary M1") for k in ("CalibPrimMOne", "CalibPrimMOneInt")})
    # 補足の図の番号（現れる順: S1 事後の確認、S2 地形の梯子、S3 集落の近くの他の推定、S4 名前の描くもの）
    out |= {"STabFigNearMore": "S3", "STabFigNames": "S4", "STabFigPlasmode": "S5", "STabFigSpec": "S6"}
    return out


def table_ceiling() -> None:
    """本文の表 3（第 4.1 節）: 集落を通じた割り当てが運べるもの。尺度は全て R（セル全体の Brier の変化 ÷ 基礎の利得）。"""
    rows = [
        r"\multicolumn{3}{l}{\emph{What the ground nearest water-side names contributes to $R$}} \\",
        (r"Observed contribution (90\,\% interval)$^{a}$ & \ObsRAllWet{} (\ObsRAllWetLo{} to \ObsRAllWetHi) & "
         r"\ObsRFiveHundred{} (\ObsRFiveHundredLo{} to \ObsRFiveHundredHi) \\"),
        r"Largest contribution if names predicted this ground perfectly & \CeilRAll & \CeilRFiveHundred \\",
        (r"\multicolumn{3}{l}{\emph{Other information in place of the names ($R$, 90\,\% interval; all cells)}} "
         r"\\"),
        (r"Flooding of the cell itself in other EMS events$^{b}$ & \multicolumn{2}{l}{\RecurCellR{} "
         r"(\RecurCellRLoNinety{} to \RecurCellRHiNinety)} \\"),
        (r"\quad the same, averaged around the settlement and assigned through it & \multicolumn{2}{l}{\RecurSetR{} "
         r"(\RecurSetRLoNinety{} to \RecurSetRHiNinety)} \\"),
        (r"\quad the same, on the \RecurSubNAoi{} AOIs where the record exists$^{c}$ & \multicolumn{2}{l}{\RecurSubR{} "
         r"(\RecurSubRLoNinety{} to \RecurSubRHiNinety)} \\"),
        (r"Recurrent surface water of the cell itself, 1984--2021 & \multicolumn{2}{l}{\GswCellR{} "
         r"(\GswCellRLoNinety{} to \GswCellRHiNinety)} \\"),
        (r"\quad the same, around the settlement and assigned through it & \multicolumn{2}{l}{\GswOracleR{} "
         r"(\GswOracleRLoNinety{} to \GswOracleRHiNinety)} \\"),
        r"HAND of the settlement's site & \multicolumn{2}{l}{\OracleR{} (\OracleRLoNinety{} to \OracleRHiNinety)} \\",
        (r"Flooded share around the settlement in the same event$^{d}$ & \multicolumn{2}{l}{\LabOracleR{} "
         r"(\LabOracleRLoNinety{} to \LabOracleRHiNinety)} \\"),
    ]
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[tp]\n"
           "\\caption{What assignment through settlements can register (exploratory; amendments 2 and 5 to 8). All "
           "quantities are on the scale of the pooled $R$: the change in the Brier score of all cells divided by the "
           "baseline gain. The first block concerns the ground nearest a settlement with a water-side name "
           "(\\MatchWetPct{} of the cells) and the ground within 500\\,m of such a settlement; the second replaces the "
           "name features with other information about each cell or settlement, assigned to every cell. "
           "$^{a}$~The increments on this ground (amendment 2) re-expressed on the pooled scale (not planned). "
           "$^{b}$~Known for \\RecurCovPct{} of the cells; zero with an indicator elsewhere. $^{c}$~Against that "
           "subset's baseline gain. $^{d}$~Leaks the outcome by construction (the cell itself is left out): it shows "
           "what the design can register, not a possible reading of names.}\\label{tab:ceiling}\n\\scriptsize\n"
           "\\begin{tabular}{lll}\n\\tophline\n"
           " & All ground nearest the name & Within 500\\,m of the settlement \\\\\n\\middlehline\n"
           + "\n".join(rows) + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    (HERE / "tab_ceiling.tex").write_text(tex, encoding="utf-8")


def table_bounds() -> None:
    """本文の表 4（第 4.2 節）: 余白と上限をオッズ比で。注入の名目の尺度と第 5 節の推定量の尺度を別の列に。"""
    rows = [
        (r"\multicolumn{5}{l}{\emph{Injected effects (base models fixed; the baseline model's held-out prediction "
         r"taken as the true probability)}} \\"),
        (r"Odds ratio at the margin & \PlasMarginA & \MarginEstA & \ensuremath{>}\,20 & --- \\"),
        (r"Upper bound, fitted models held fixed (Monte Carlo range)$^{a}$ & \BoundNewA{} (\BoundNewALo--\BoundNewAHi) "
         r"& \BoundEstA & \BoundNewB{} (\BoundNewBLo--\BoundNewBHi) & \BoundEstBText \\"),
        r"\quad with the spread over fold assignments & \BoundFoldA & \BoundFoldEstA & \BoundFoldB & \BoundEstBText \\",
        (r"\quad with the second stage refitted within the bootstrap & \BoundRefA & \BoundRefEstA & \BoundRefB & "
         r"\BoundEstBText \\"),
        r"\quad for truly water-related names$^{b}$ & \BoundNewCA & --- & \BoundNewCB & --- \\",
        r"\quad flips correlated within the ground of each settlement & \BoundACA & \BoundACEstA & --- & --- \\",
        r"Lower bound for effects that remove flooding & \BoundNegA & --- & \BoundNegB & --- \\",
        (r"Detectable with 80\,\% power, normal approximation (fixed; refitted) & \MdeNewA; \MdeRefA & --- & "
         r"\MdeNewB; \MdeRefB & --- \\"),
        (r"\multicolumn{5}{l}{\emph{Direct estimates: water-side against other names, given the baseline model "
         r"(95\,\% interval)}} \\"),
        (r"EMS labels, capped sample$^{c}$ & --- & \ORWetAllAdj{} (\ORWetAllAdjLo--\ORWetAllAdjHi) & --- & "
         r"\ORWetNearAdj{} (\PairsWetNearLo--\PairsWetNearHi) \\"),
        r"EMS labels, every evaluable cell$^{d}$ & --- & --- & --- & \NearUncapOR{} (\NearUncapORLo--\NearUncapORHi) \\",
        r"GFM labels, every assessable cell$^{d}$ & --- & --- & --- & \GfmUncapOR{} (\GfmUncapORLo--\GfmUncapORHi) \\",
        r"\quad EMS labels on the same cells$^{d}$ & --- & --- & --- & \EmsSameOR{} (\EmsSameORLo--\EmsSameORHi) \\",
    ]
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[tp]\n"
           "\\caption{Odds ratios of flooding that the margin and the observed increment correspond to, and the direct "
           "estimates (exploratory; amendments 6 to 12). Each column has one scale. Nominal: the odds ratio injected "
           "into the labels. Estimator: the same effect on the scale of the odds ratios of Sect.~\\ref{sec:near}, read "
           "from the estimator's values on labels with known injected effects (interpolated on the log scale; beyond the "
           "largest effect used for this conversion, an injected odds ratio of 3, written as greater than the estimator's "
           "value there). No-effect reference: the mean "
           "increment of reassigned names with the regional composition rebuilt (amendment 8); bounds with other "
           "references are in Tables~\\STabEightBounds{}, \\STabNineBounds{} and~\\STabElevenScale{} of the Supplement. $^{a}$~Monte Carlo 2.5th to 97.5th "
           "percentiles over the simulated replicates (\\PlasRepsA{} each), which reflect simulation error only. "
           "$^{b}$~Effect injected on a latent indicator of truly water-related names, drawn with the misreading rates "
           "measured against recorded etymologies. $^{c}$~Within 500\\,m, pairs cluster bootstrap over location "
           "groups (2,000 replicates); on all ground, cluster-robust interval. $^{d}$~Pairs cluster bootstrap "
           "(2,000 replicates).}\\label{tab:bounds}\n\\scriptsize\n"
           "\\begin{tabular}{lllll}\n\\tophline\n"
           " & \\multicolumn{2}{l}{All ground nearest the name} & \\multicolumn{2}{l}{Within 500\\,m of the settlement}"
           " \\\\\n & Nominal & Estimator & Nominal & Estimator \\\\\n\\middlehline\n"
           + "\n".join(rows) + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    (HERE / "tab_bounds.tex").write_text(tex, encoding="utf-8")


AMENDMENT_ROWS = (
    ("BOne", "1", ("Corrects the implementation of sensitivity analysis 6 (location groups from AOIs buffered by 10\\,km); "
                  "both runs reported.")),
    ("BTwo", "2", ("Exploratory analyses after critiques by language-model agents with access to the records: positive "
                  "control, alternative measures of the baseline gain, two-way bootstrap, distance bands, odds ratios "
                  "near settlements, fold assignments, events held out whole, injected effects.")),
    ("BThree", "3", ("Independent flood labels from GFM: label rules, odds ratios near settlements, refitting (designated "
                    "replication) and rescoring of the primary comparison.")),
    ("BFour", "4", ("Corrects the reading of the lexicon for Norway; reruns the confirmatory analysis; keeps the "
                   "confirmatory run primary if the verdict is the same; both runs reported.")),
    ("BFive", "5", ("After the first critiques of the manuscript: permutations, second-stage coefficients, thirteen fold "
                   "assignments, a terrain oracle, within-area metrics, label variants, element groups, recorded "
                   "etymologies.")),
    ("BSix", "6", ("After the second critiques: ceiling of local effects, injected effects on the observed labels, "
                  "flood-information controls, GFM diagnostics, a strict reading, the protocol's $p$ values.")),
    ("BSeven", "7", ("After the third critiques: permutations with recomputed composition, refitted base models, the "
                    "activation split with the corrected code, recurrence of flooding in other events, a floodplain "
                    "reading, partly flooded cells, residual correlation between AOIs of one event.")),
    ("BEight", "8", ("After the fourth critiques: the no-effect reference of the bounds, the achieved equivalence bound, "
                    "the second stage refitted inside the bootstrap, recurrent surface water of each cell, a cleaned "
                    "lexicon, heights above the main rivers, a second coder.")),
    ("BNine", "9", ("After the fifth critiques: bounds with refitted intervals and over fold assignments, small and "
                   "protective injected effects, effects that differ between settlements, simulated power with every "
                   "model refitted, flood-mapping skill, one definition of permanent water, other subsampling seeds, a "
                   "bare-earth elevation model.")),
    ("BTen", "10", ("After the sixth critiques: reframing and corrections; detectable effects and further intervals near "
                   "settlements, misreading over a grid of rates, the islands-and-meadows group at settlement sites, the "
                   "floodplain reading by country and form, sensors, the activation split in full, distance and "
                   "continent interactions, a tuned penalty, country indicators, spatially clustered injection, modelled "
                   "flood hazard maps as baseline.")),
    ("BEleven", "11", ("After the seventh critiques: refocusing on the scheme and the within-area estimand, and "
                      "corrections; the change of the activation split between code versions, uncapped odds ratios "
                      "near settlements on both label sets, bounds on the estimator's scale, the GFM domain, the "
                      "denominator of the labels, corrected and combined intervals, a census of the exploratory "
                      "estimates, the gazetteer's feature codes and rounded coordinates, the pilots, a new audit sample, "
                      "and permutations with rebuilt composition on the country split.")),
    ("BTwelve", "12", ("After the eighth critiques: restructuring around the pre-specified test, the uncertainty its "
                      "interval leaves out, what assignment through settlements can carry and the limits of readings "
                      "and labels; operating characteristics of two decision rules, the leak of the primary split and "
                      "the reference model, constant forecasts, a specification curve of the odds ratios near "
                      "settlements, the composition of the uncapped sample, languages recounted, an injection within "
                      "1\\,km, built-up land near settlements, flood typing and the GFM window, a second coder's "
                      "intervals and a revised audit sample.")),
)


def table_amendments(m: dict[str, str]) -> None:
    rows = [f"{k} & {m[f'AmendTime{key}'] if f'AmendTime{key}' in m else missing(f'amendment time {key}')} & {txt} \\\\"
            for key, k, txt in AMENDMENT_ROWS]
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[tp]\n"
           "\\caption{Amendments to the frozen protocol. Times are the commit times of the amendment texts in the "
           "development repository (UTC); Table~\\ref{tab:timeline} gives the commit time of every addendum.}"
           "\\label{tab:amendments}\n\\footnotesize\n\\begin{tabular}{llp{11cm}}\n\\tophline\n"
           "Amendment & Written & Content \\\\\n\\middlehline\n" + "\n".join(rows)
           + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    (HERE / "app_table_amendments.tex").write_text(tex, encoding="utf-8")


def fig_near(v68: dict | None, v57: dict | None) -> None:
    """本文の図 6（第 5 節）: 集落の近くの水系側の OR。上限なしの全セル（B-11 の B2）を主に、区間は全て pairs。"""
    v = v68 or {}
    b2, b5u = v.get("B2") or {}, v.get("B5u") or {}
    cap = ((v57 or {}).get("B6") or {}).get("main")
    cb = (((_json("v50_gfm_diagnostics.json") or {}).get("contributing") or {}).get("ems_all_aoi")) or {}
    if isinstance(cap, dict) and cb:  # 上限の表の行も、浸水セルのある場所グループの数で（他の行と同じ数え方）
        cap = cap | {"n_flooded_near": cb["n_flooded_group"], "n_units_flooded_near": cb["n_units_with_flooded"]}
    if isinstance(cap, dict):
        cap = cap | {"n_aoi": cap.get("n_aoi")}
    # B-12 の A7: 計画した GFM の推定（上限をかけた標本、amendment 3）と同じセルの EMS を、pairs の区間で（v57 の B5）
    sc = ((v57 or {}).get("B5") or {}).get("same_cells") or {}
    cg = (_json("v50_gfm_diagnostics.json") or {}).get("contributing") or {}
    pl = {}
    for key, ck in (("gfm", "gfm_gfm_cells"), ("ems", "ems_gfm_cells")):
        d, c = sc.get(key) or {}, cg.get(ck) or {}
        if d.get("pairs_ci95"):
            pl[key] = d | {"n_flooded_near": c.get("n_flooded_group"), "n_units_flooded_near": c.get("n_units_with_flooded"),
                           "n_aoi": sc.get("n_aoi")}
    blue, green, grey = "#1f6fb4", "#2b8c4a", "0.45"
    groups = [
        ("EMS labels", [("Every evaluable cell (chosen)", b2.get("ems_all"), blue, True),
                        ("Capped sample", cap, grey, False),
                        ("Observed denominator, every cell", b5u.get("or_uncapped"), blue, False)]),
        ("Cells that GFM can assess", [("GFM labels, capped sample (planned)", pl.get("gfm"), green, "D"),
                                       ("EMS labels, same cells, capped", pl.get("ems"), grey, False),
                                       ("GFM labels, every cell (chosen)", b2.get("gfm"), green, True),
                                       ("EMS labels, same cells, every cell", b2.get("ems_same_cells"), blue, False)]),
        ("AOIs not mapped from Sentinel-1", [("GFM labels", b2.get("gfm_not_s1"), green, True),
                                             ("EMS labels, same cells", b2.get("ems_not_s1"), blue, False)]),
    ]
    n = sum(len(r) + 1 for _, r in groups)
    fig, ax = plt.subplots(figsize=(6.6, 0.24 * n + 0.75), layout="constrained")
    y, ticks, labels, src = n - 1, [], [], []
    for title, rows in groups:
        ax.text(0.205, y, title, fontsize=7, fontweight="bold", va="center", ha="left")
        y -= 1
        for lab, d, col, filled in rows:
            ticks.append(y)
            labels.append(lab)
            if not isinstance(d, dict) or not d.get("pairs_ci95"):
                ax.text(1.0, y, "pending", fontsize=6, va="center", ha="center", color="0.5")
                y -= 1
                continue
            orr, ci = d["odds_ratio"], d["pairs_ci95"]
            lo, hi = max(ci[0], 0.2), min(ci[1], 5.0)
            ax.plot([lo, hi], [y, y], color=col, lw=1.4)
            ax.plot([orr], [y], "D" if filled == "D" else "o", ms=4.6 if filled == "D" else 4.2, color=col,
                    mfc=col if filled else "white", zorder=3)
            nf = d.get("n_flooded_near", d.get("n_group_flooded"))
            nu = d.get("n_units_flooded_near", d.get("n_group_units"))
            if nf is not None:
                ax.text(5.6, y, f"{int(nf):,} / {nu} / {d.get('n_aoi', '')}", fontsize=6, va="center", ha="left",
                        color="0.3", clip_on=False)
            src.append({"group": title, "estimate": lab, "odds_ratio": orr, "pairs_ci95": list(ci),
                        "flooded_near": nf, "groups_with_flooded_near": nu})
            y -= 1
    ax.text(5.6, n - 1, "flooded / groups / AOIs", fontsize=6, va="center", ha="left", color="0.3", clip_on=False)
    ax.axvline(1, color="0.4", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 5.0)
    ax.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.25", "0.5", "1", "2", "4"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel("odds ratio of flooding within 500 m of water-side names, given the baseline model's prediction",
                  fontsize=7.5)
    ax.tick_params(axis="x", labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    g.SOURCE["FigNear"] = src
    fig.savefig(HERE / "fig_near.pdf")
    plt.close(fig)




def table_eleven_headings() -> None:
    """補足の表: Wikipedia の語源の節として取った見出しの語（言語ごと、v44 の HEAD）。初見の査読 第 7 回（査読者 3 の 16）。"""
    spec = importlib.util.spec_from_file_location("v44h", REPO / "backend" / "verification" / "r_okauchi"
                                                  / "v44_wikipedia_etymology.py")
    src = (REPO / "backend" / "verification" / "r_okauchi" / "v44_wikipedia_etymology.py").read_text(encoding="utf-8")
    head_src = src[src.index("HEAD = {"):src.index("\n}\n", src.index("HEAD = {")) + 3]
    strong_src = src[src.index("STRONG = ("):src.index(")\n", src.index("STRONG = (")) + 2]
    ns: dict = {}
    exec(head_src + strong_src, ns)  # noqa: S102  定数の辞書だけを読む（モジュールは網から取るので import しない）
    del spec
    latin, other = [], []
    for lang, words in ns["HEAD"].items():
        if all(ord(ch) < 0x250 for w in words for ch in w):
            esc = ", ".join(w.replace("&", "\\&") for w in words)
            latin.append(f"{lang} & {_latex_accents(esc)} \\\\")
        else:
            other.append(lang)
    rows = latin + [f"\\multicolumn{{2}}{{p{{14cm}}}}{{Languages in other scripts, listed in the code: "
                    f"{', '.join(other)}.}} \\\\"]
    strong = ", ".join(w for w in ns["STRONG"] if all(ord(ch) < 0x250 for ch in w))
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[tp]\n"
           "\\caption{Headings under which the check against recorded etymologies took a Wikipedia section as an "
           "explanation of the name (amendment 5): the heading, lower-cased, equals one of these words or begins with it "
           "followed by a conjunction or punctuation; a heading that contains any of the following anywhere also counts: "
           + _latex_accents(strong) + ". Sections under other headings, such as history, were not taken.}"
           "\\label{tab:eleven-headings}\n\\footnotesize\n\\begin{tabular}{lp{14cm}}\n\\tophline\n"
           "Language & Heading words \\\\\n\\middlehline\n" + "\n".join(rows)
           + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    (HERE / "app_table_eleven_headings.tex").write_text(tex, encoding="utf-8")


def _latex_accents(t: str) -> str:
    """見出しの語のアクセント付きの文字を LaTeX の命令に（T1 の lmodern で組めるように）。"""
    rep = {"é": "\\'e", "è": "\\`e", "á": "\\'a", "í": "\\'i", "ó": "\\'o", "ú": "\\'u", "ñ": "\\~n",
           "ã": "\\~a", "õ": "\\~o", "ç": "\\c{c}", "ä": "\\\"a", "ö": "\\\"o", "ü": "\\\"u", "å": "\\aa{}",
           "ø": "\\o{}", "æ": "\\ae{}", "ő": "\\H{o}", "ű": "\\H{u}", "č": "\\v{c}", "š": "\\v{s}", "ž": "\\v{z}",
           "ř": "\\v{r}", "ů": "\\r{u}", "ý": "\\'y", "ě": "\\v{e}", "ł": "\\l{}", "ś": "\\'s", "ź": "\\'z",
           "ż": "\\.{z}", "ą": "\\k{a}", "ę": "\\k{e}", "ć": "\\'c", "ń": "\\'n", "ş": "\\c{s}", "ı": "\\i{}",
           "ğ": "\\u{g}", "ā": "\\=a", "ē": "\\=e", "ī": "\\=i", "ū": "\\=u",
           "ķ": "\\c{k}", "ļ": "\\c{l}", "ņ": "\\c{n}", "ģ": "\\c{g}", "ė": "\\.e", "į": "\\k{i}", "ų": "\\k{u}",
           "ț": "\\c{t}", "ș": "\\c{s}", "ă": "\\u{a}", "â": "\\^a", "î": "\\^i", "ê": "\\^e", "ô": "\\^o",
           "đ": "\\dj{}", "ạ": "a", "ừ": "u", "ố": "o", "ọ": "o", "ồ": "o", "à": "\\`a", "ò": "\\`o", "ì": "\\`i",
           "ù": "\\`u", "–": "--"}
    return "".join(rep.get(ch, ch) for ch in t)



# %% B-12（初見の査読 第 8 回）: 本文の組み直しで使う量・表・図
def b12_text_macros(v69: dict | None, v68: dict | None, v57: dict | None, v41: dict | None, v50: dict | None,
                    v54: dict | None, v44: dict | None, v43: dict | None, v36: dict | None) -> dict[str, str]:
    """B-12 の B（v69）と既存の成果物から、本文が引く量。"""
    v, w = v69 or {}, v68 or {}
    out: dict[str, str] = {}

    def need(keys: tuple[str, ...], why: str) -> dict[str, str]:
        return {k: missing(why) for k in keys}

    # 方法: どちらでもないセル、種別、GFM の時間窓
    b14 = v.get("B14") or {}
    out["NeitherCells"] = _thousands(b14["n_neither"]) if b14 else missing("v69 B14")
    b13 = v.get("B13") or {}
    if b13:
        out["NUntypedAoi"] = str(b13["aoi_with_untyped"])
        out["NIceJamAoi"] = str(sum(1 for r in b13.get("typed_other_aois", [])
                                    if any("ice" in d.lower() for d in r["descriptions"]) and r["inland_rule"]))
        out["TraceOpticalPct"] = _pct0(b13["trace_share_by_sensor"].get("optical") or 0.0)
    else:
        out |= need(("NUntypedAoi", "NIceJamAoi", "TraceOpticalPct"), "v69 B13")
    b12w = v.get("B12") or {}
    out |= ({"GfmWinAoi": str(b12w["n_aoi"]), "GfmWinSameDay": str(b12w["aoi_with_an_image_on_an_ems_day"]),
             "GfmWinAdjOnly": str(b12w["aoi_with_images_only_on_adjacent_days"])} if b12w
            else need(("GfmWinAoi", "GfmWinSameDay", "GfmWinAdjOnly"), "v69 B12"))
    # 目録（直した族と BY）
    c8 = (v.get("B8c") or {}).get("near_settlement_or") or {}
    old = (w.get("B8") or {}).get("near_settlement_or") or {}
    out |= ({"CensusNearN": str(c8["n"]), "CensusNearExcl": str(c8["n_excluding_one"]),
             "CensusNearBHMin": f"{min(r['q_bh'] for r in c8['smallest_q']):.2f}",
             "CensusNearBYMin": f"{min(r['q_by'] for r in c8['smallest_q']):.2f}",
             "CensusNearNOld": str(old.get("n", "")) or missing("v68 B8")} if c8
            else need(("CensusNearN", "CensusNearExcl", "CensusNearBHMin", "CensusNearBYMin", "CensusNearNOld"),
                      "v69 B8c"))
    # 名前の特徴の実質の数（一定の 1 列を除く）
    # 名前の特徴のうち、解析した全てのセルで一定の列を除いた数（直したコードの二段目の設計。fold 0 の訓練と評価で全セル）
    z2 = DATA / "runs" / "world_v2_b8stage2_local" / "predictions" / "primary" / "fold0_stage2.npz"
    if z2.is_file():
        zz = np.load(z2, allow_pickle=True)
        xn = np.vstack([zz["xn_tr"], zz["xn_te"]])
        n_const = int(sum(np.nanmin(xn[:, i]) == np.nanmax(xn[:, i]) for i in range(xn.shape[1])))
        out["NNameFeatEff"] = str(xn.shape[1] - n_const)
    else:
        out["NNameFeatEff"] = missing("stage-2 dump")
    bc = (((_json("v66_b10_local.json") or {}).get("B15") or {}).get("by_cells") or {}).get("unknown") or {}
    if bc.get("share") is not None:
        out["NReadUnknownPct"] = pct1(bc["share"])
    # 結果: 大陸ごとのセルとラベルの分散の割合（p(1-p) で重み）
    st = ((((v36 or {}).get("splits") or {}).get("primary") or {}).get("strata") or {})
    cont = {k.split("_")[1]: x for k, x in st.items() if k.startswith("continent_")}
    if cont:
        n_all = sum(x["n_cells"] for x in cont.values())
        var = {k: x["n_cells"] * x["wet_rate"] * (1 - x["wet_rate"]) for k, x in cont.items()}
        v_all = sum(var.values())
        out |= {"ShareCellsAsia": _pct0(cont["AS"]["n_cells"] / n_all), "ShareVarAsia": _pct0(var["AS"] / v_all),
                "ShareCellsEurope": _pct0(cont["EU"]["n_cells"] / n_all), "ShareVarEurope": _pct0(var["EU"] / v_all)}
    else:
        out |= need(("ShareCellsAsia", "ShareVarAsia", "ShareCellsEurope", "ShareVarEurope"), "v36 strata")
    llm = g._json(DATA / "runs" / "world_v2" / "metrics_llmcut.json") or {}
    fx = llm.get("flexible") or {}  # B-12 の A12: 言語モデルの分割の柔軟な比較（辞書と文字モデル、M1 → M3）の判定
    out["VerdictFlexLlm"] = fx["decision"] if fx.get("decision") else missing("llmcut flexible")
    h4 = llm.get("H4") or {}
    out |= ({"RHFourLoNinety": rv(h4["ci90_R"][0]), "RHFourHiNinety": rv(h4["ci90_R"][1])} if h4.get("ci90_R")
            else need(("RHFourLoNinety", "RHFourHiNinety"), "llmcut H4"))
    # 第 4 節: 組み合わせた区間、達成した境界、fold の割り当て、動作特性
    b7 = w.get("B7") or {}
    tw, gr = b7.get("two_way") or {}, b7.get("groups") or {}
    out |= ({"CombTwoWayLo": rv(tw["ci90_R_refit_ratio"][0]), "CombTwoWayHi": rv(tw["ci90_R_refit_ratio"][1]),
             "CombGroupsLo": rv(gr["ci90_R_refit_ratio"][0]), "CombGroupsHi": rv(gr["ci90_R_refit_ratio"][1]),
             "RatioGroupsLo": rv(gr["ci90_R_fixed_ratio"][0]), "RatioGroupsHi": rv(gr["ci90_R_fixed_ratio"][1])}
            if tw.get("ci90_R_refit_ratio") and gr.get("ci90_R_refit_ratio")
            else need(("CombTwoWayLo", "CombTwoWayHi", "CombGroupsLo", "CombGroupsHi", "RatioGroupsLo",
                       "RatioGroupsHi"), "v68 B7"))
    mx = ((v57 or {}).get("B2") or {}).get("max_over_combos") or {}
    out["AchievedBoundMax"] = rv(mx["achieved_bound_R"]) if mx else missing("v57 B2 max")
    fs = ((_json("v47_b5_cluster.json") or {}).get("fold_seeds") or {}).get("fixed_code") or {}
    out |= ({"FoldRMin": rv(fs["R_min"]), "FoldRMax": rv(fs["R_max"])} if fs.get("R")
            else need(("FoldRMin", "FoldRMax"), "v47 fold seeds"))
    b1 = v.get("B1") or {}
    se = b1.get("share_equivalent") or {}
    if se:
        c, p = se["constant"], se["proportional"]
        out |= {"OpCorr": f"{b1['corr_delta_yardstick']:.2f}",
                "OpFixedSize": _pct0(max(c["+0.05"]["fixed"], c["-0.05"]["fixed"])),
                "OpRatioSize": _pct0(max(c["+0.05"]["ratio"], c["-0.05"]["ratio"])),
                "OpFixedSizeProp": _pct0(max(p["+0.05"]["fixed"], p["-0.05"]["fixed"])),
                "OpRatioSizeProp": _pct0(max(p["+0.05"]["ratio"], p["-0.05"]["ratio"])),
                "OpPowerMin": _pct0(min(c["+0.00"]["fixed"], c["+0.00"]["ratio"], p["+0.00"]["fixed"],
                                        p["+0.00"]["ratio"]))}
    else:
        out |= need(("OpCorr", "OpFixedSize", "OpRatioSize", "OpFixedSizeProp", "OpRatioSizeProp", "OpPowerMin"),
                    "v69 B1")
    b5 = (v.get("B5") or {}).get("outlying") or []
    out |= ({"NullOutPk": _word(sum(r["top_unit_country"] == "PK" for r in b5)),
             "NullOutBd": _word(sum(r["top_unit_country"] == "BD" for r in b5))} if b5
            else need(("NullOutPk", "NullOutBd"), "v69 B5"))
    # 漏れと基準のモデル、定数の予測
    b10 = v.get("B10") or {}
    out |= ({"RefGainShare": pct1(b10["gain_M1c_over_M1s_share_of_yardstick"]),
             "RefNoNlearnShare": pct1(-b10["gain_M1c_no_nlearn_over_M1s_share_of_yardstick"]),
             "ResActCorrNearLo": f"\\ensuremath{{{b10['residual_corr']['ci95_fisher'][0]:.2f}}}",
             "ResActCorrNearHi": f"\\ensuremath{{{b10['residual_corr']['ci95_fisher'][1]:.2f}}}"} if b10
            else need(("RefGainShare", "RefNoNlearnShare", "ResActCorrNearLo", "ResActCorrNearHi"), "v69 B10"))
    b17 = v.get("B17") or {}
    act, ctry = b17.get("activation_frozen") or {}, b17.get("country") or {}
    out |= ({"ActBSMZero": f"{act['brier_M0']:.4f}", "ActBSConst": f"{act['brier_best_constant']:.4f}",
             "CountryBSMZero": f"{ctry['brier_M0']:.4f}", "CountryBSConst": f"{ctry['brier_best_constant']:.4f}"}
            if act.get("brier_M0") and ctry.get("brier_M0")
            else need(("ActBSMZero", "ActBSConst", "CountryBSMZero", "CountryBSConst"), "v69 B17"))
    b11 = v.get("B11") or {}
    if b11.get("R") is not None:
        ok = all(all(x) for x in b11["interior"].values())
        out["TunedNamesText"] = (f"$R$ was {rv(b11['R'])}" + ("" if ok else
                                 " (the optimum lay at the edge of the widened grid in some folds)"))
    else:
        out["TunedNamesText"] = missing("v69 B11")
    # 第 5 節: 再起の記録の全体の尺度、1 km の注入
    b7r = v.get("B7r") or {}
    out["RecurPooledR"] = rv(b7r["pooled_scale_R"]) if b7r else missing("v69 B7r")
    c1k = (v.get("C1") or {}) if isinstance(v.get("C1"), dict) else {}
    out["MarginOneKmText"] = (c1k["text"] if c1k.get("text") else missing("B-12 C1"))
    # 第 6 節: 読み
    out |= ({"WDNamesakes": str(v43["n_with_p138"]), "WDRiverNamesakes": str(v43["n_hydro"]),
             "WDNreadHits": str(v43["recall_nread"]["k"]), "WDLexHits": str(v43["recall_lexicon"]["k"])} if v43
            else need(("WDNamesakes", "WDRiverNamesakes", "WDNreadHits", "WDLexHits"), "v43"))
    out["RiversidePct"] = pct1((v54 or {})["riverside"]["share_cells"]) if v54 else missing("v54 riverside")
    b3 = v.get("B3") or {}
    out |= ({"OtherLangCountryPct": pct1(b3["country_level_share"]), "OtherLangRegionPct": pct1(b3["regional_share_v66"]),
             "OtherLangOldPct": pct1(b3["country_level_share_v66"] + b3["regional_share_v66"])} if b3
            else need(("OtherLangCountryPct", "OtherLangRegionPct", "OtherLangOldPct"), "v69 B3"))
    p0 = ((v44 or {}).get("truth_rate_in_not") or {})
    out["PZeroPct"] = _pct0(p0["rate"]) if p0 else missing("v44 p0")
    b21, b21b = v.get("B21") or {}, v.get("B21b") or {}
    out |= ({"AuditNTotal": str(b21b["n_names"]), "AuditNAdded": str(b21b["n_added"]),
             "AuditNNot": str(b21["n_not_read"])} if b21 and b21b
            else need(("AuditNTotal", "AuditNAdded", "AuditNNot"), "v69 B21"))
    met = g._json(DATA / "runs" / "world_v2" / "metrics_primary.json") or {}
    h5 = met.get("H5_strata") or {}
    out["LexAnyPct"] = (pct1(h5["nearest_name_in_lexicon"]["n"] / met["n_cells"])
                        if h5.get("nearest_name_in_lexicon") and met.get("n_cells") else missing("lexicon any"))
    # 集落の近くのラベル、GFM、OR
    b18 = v.get("B18") or {}
    wb = b18.get("with_built10") or {}
    out["BuiltUpBandText"] = (
        "with the share of built-up land in each cell at 10\\,m (ESA WorldCover) added, these ratios fell to "
        f"{orv(wb['band_500m_1km']['odds_ratio'])} ({orv(wb['band_500m_1km']['ci95'][0])}--"
        f"{orv(wb['band_500m_1km']['ci95'][1])}) and {orv(wb['band_gt1km']['odds_ratio'])} "
        f"({orv(wb['band_gt1km']['ci95'][0])}--{orv(wb['band_gt1km']['ci95'][1])})") if wb else missing("v69 B18")
    c1e = (v41 or {}).get("C1_ems_lexfix") or {}
    ci = c1e.get("ci90_R") or c1e.get("ci90")
    out |= ({"GfmDomRFixLo": rv(ci[0]), "GfmDomRFixHi": rv(ci[1])} if ci
            else need(("GfmDomRFixLo", "GfmDomRFixHi"), "v41 C1_ems_lexfix"))
    pl = (((v41 or {}).get("odds_gfm") or {}).get("wetside_adjusted_M1") or {}).get("terms", {}).get("wetside:le500m")
    out |= ({"GFMORWetNearPlanned": orv(pl["odds_ratio"]), "GFMORWetNearPlannedLo": orv(pl["ci95_odds_ratio"][0]),
             "GFMORWetNearPlannedHi": orv(pl["ci95_odds_ratio"][1])} if pl
            else need(("GFMORWetNearPlanned", "GFMORWetNearPlannedLo", "GFMORWetNearPlannedHi"), "v41 planned"))
    cg = ((v50 or {}).get("contributing") or {}).get("gfm_gfm_cells") or {}
    out |= ({"GfmCappedFlooded": str(cg["n_flooded_group"]), "GfmCappedUnits": str(cg["n_units_with_flooded"])} if cg
            else need(("GfmCappedFlooded", "GfmCappedUnits"), "v50 contributing"))
    b19 = v.get("B19") or {}
    rg = b19.get("ranges") or {}
    if all(k in rg for k in ("GFM_gfm_cells", "EMS_gfm_cells", "EMS_all")):  # 追記 2・3 で直した区分
        out |= {"SpecN": str(b19["n"]),
                "SpecGfmMin": orv(rg["GFM_gfm_cells"]["min"]), "SpecGfmMax": orv(rg["GFM_gfm_cells"]["max"]),
                "SpecEmsGfmMin": orv(rg["EMS_gfm_cells"]["min"]), "SpecEmsGfmMax": orv(rg["EMS_gfm_cells"]["max"]),
                "SpecEmsMin": orv(rg["EMS_all"]["min"]), "SpecEmsMax": orv(rg["EMS_all"]["max"])}
    else:
        out |= need(("SpecN", "SpecGfmMin", "SpecGfmMax", "SpecEmsGfmMin", "SpecEmsGfmMax", "SpecEmsMin", "SpecEmsMax"),
                    "v69 B19")
    b2 = v.get("B2") or {}
    if b2.get("by_country"):
        top = b2["by_country"]
        kh_mm = sum(r["share_of_flooded_near"] for r in top if r["cc"] in {"KH", "MM"})
        pkbd = b2.get("share_pk_bd") or 0.0
        nb = b2.get("without_pk_bd") or {}
        ta = (b2.get("by_activation") or [{}])[0]
        inc = (g._json(DATA / "global" / "incidence.json") or {}).get("aois") or []
        tc = sorted({c for r in inc if r.get("activation") == ta.get("activation") for c in (r.get("countries") or [])})
        tcn = {"Cambodia": "Cambodian", "Myanmar": "Myanmar", "Bangladesh": "Bangladeshi", "Pakistan": "Pakistani"}
        tlab = f"one {tcn.get(tc[0], tc[0])} activation" if len(tc) == 1 else "one activation"
        out["NearCompositionText"] = (
            f"In the sample of every evaluable cell, Cambodia and Myanmar hold {pct1(kh_mm)} of the flooded cells within "
            f"500\\,m of water-side names ({tlab} alone {pct1(ta.get('share_of_flooded_near', 0))}) and "
            f"Pakistan and Bangladesh {pct1(pkbd)}"
            + (f"; without Pakistan and Bangladesh the estimate was {orv(nb['odds_ratio'])} "
               f"({orv(nb['pairs_ci95'][0])}--{orv(nb['pairs_ci95'][1])})" if nb.get("pairs_ci95") else "")
            + " (Table~\\STabTwelveChecks{} of the Supplement).")
    else:
        out["NearCompositionText"] = missing("v69 B2")
    b6 = v.get("B6") or {}
    hg, he = b6.get("gfm_not_s1") or {}, b6.get("ems_not_s1") or {}
    sa = b6.get("share_assessable_le500m") or {}
    out |= ({"HighNAoi": str(b6["n_aoi_not_s1"]), "HighGfmOR": orv(hg["odds_ratio"]),
             "HighGfmORLo": orv(hg["ci95_odds_ratio"][0]), "HighGfmORHi": orv(hg["ci95_odds_ratio"][1]),
             "HighEmsOR": orv(he["odds_ratio"]), "HighEmsORLo": orv(he["ci95_odds_ratio"][0]),
             "HighEmsORHi": orv(he["ci95_odds_ratio"][1]), "HighAssessPct": _pct0(sa["highland"]),
             "WsAssessPct": _pct0(sa["water_side"])} if hg and he and sa
            else need(("HighNAoi", "HighGfmOR", "HighGfmORLo", "HighGfmORHi", "HighEmsOR", "HighEmsORLo",
                       "HighEmsORHi", "HighAssessPct", "WsAssessPct"), "v69 B6"))
    # 付録: 改訂の数（批評に応じたもの）
    am = sorted((REPO / "docs" / "protocol" / "amendments").glob("B-*.md"))
    n_am = len(am)
    out |= {"NAmend": _word(n_am), "NAmendCrit": _word(n_am - 2)}  # B-1 と B-4 は実装の誤りの直し
    return out


def table_reading(v54: dict | None, v44: dict | None, v57b: dict | None, v69: dict | None, met: dict | None) -> None:
    """本文の表 5（B-12 の A8）: 辞書の読みの要約（大陸ごとの照合、上位の形、語源との一致、言語）。"""
    v = v69 or {}
    rows: list[str] = []
    nr = ((v54 or {}).get("nread") or {}).get("wetside_by_continent_cells") or {}
    order = [("EU", "Europe"), ("AS", "Asia"), ("AF", "Africa"), ("OC", "Oceania")]
    am = {k: nr[k] for k in ("NA", "SA") if k in nr}
    if nr:
        def amer(key: str) -> float:
            n = sum(x["n_cells"] for x in am.values())
            return sum(x[key] * x["n_cells"] for x in am.values()) / n if n else float("nan")
        lex = [_pct0(nr[k]["lexicon"]) for k, _ in order] + [_pct0(amer("lexicon"))]
        llm = [_pct0(nr[k]["llm"]) for k, _ in order] + [_pct0(amer("llm"))]
        n_all = sum(x["n_cells"] for x in nr.values())
        lex_all = sum(x["lexicon"] * x["n_cells"] for x in nr.values()) / n_all
        llm_all = sum(x["llm"] * x["n_cells"] for x in nr.values()) / n_all
        rows.append("Nearest name water-side, lexicon & " + " & ".join(lex) + f" & {pct1(lex_all)} \\\\")
        rows.append("Nearest name water-side, language model & " + " & ".join(llm) + f" & {_pct0(llm_all)} \\\\")
    h5 = (met or {}).get("H5_strata") or {}
    if h5.get("nearest_name_in_lexicon"):
        rows.append(f"Nearest name matches any lexicon element & \\multicolumn{{5}}{{l}}{{}} & "
                    f"{pct1(h5['nearest_name_in_lexicon']['n'] / met['n_cells'])} \\\\")
    rows.append("\\middlehline")
    tf = REPO / "docs" / "verification" / "v51_top_forms.csv"
    n_ws = h5.get("nearest_name_wetside", {}).get("n")
    if tf.is_file() and n_ws:
        d = pd.read_csv(tf).head(8)
        forms = ", ".join(f"\\emph{{{('-' if r.pos == 'suffix' else '') + str(r.form) + ('-' if r.pos == 'prefix' else '')}}}"
                          f" ({r.lang}, {_pct0(r.n_cells / n_ws)})" for r in d.itertuples())
        rows.append(f"Forms supplying most water-side cells (block, share) & \\multicolumn{{6}}{{p{{11.2cm}}}}{{{forms}}} \\\\")
    rows.append("\\middlehline")
    c2 = (v.get("B9") or {})
    w44 = v44 or {}

    def pci(x: dict | None) -> str:
        return f"{_pct0(x['est'])} ({_pct0(x['ci95'][0])}--{_pct0(x['ci95'][1])})" if x else "---"

    lw2, nw2 = c2.get("lex_wetside") or {}, c2.get("nread_wetside") or {}
    lw1 = w44.get("lex_wetside") or {}
    rows.append(f"Lexicon: precision, coder 2; coder 1 & \\multicolumn{{6}}{{l}}{{{pci(lw2.get('precision'))}; "
                f"{pci(lw1.get('precision'))}}} \\\\")
    rows.append(f"Lexicon: recall, coder 2; coder 1 & \\multicolumn{{6}}{{l}}{{{pci(lw2.get('recall'))}; "
                f"{pci(lw1.get('recall'))}}} \\\\")
    p0 = w44.get("truth_rate_in_not") or {}
    if p0:
        rows.append(f"Names not read as water-side with a water-related derivation, coder 1 & \\multicolumn{{6}}{{l}}{{"
                    f"{_pct0(p0['rate'])} ({_pct0(p0['ci95'][0])}--{_pct0(p0['ci95'][1])}), {p0['k']} of {p0['n']}}} \\\\")
    rows.append(f"Language model: precision; recall, coder 2 & \\multicolumn{{6}}{{l}}{{{pci(nw2.get('precision'))}; "
                f"{pci(nw2.get('recall'))}}} \\\\")
    rows.append("\\middlehline")
    b3 = v.get("B3") or {}
    if b3:
        rows.append("Cells whose country's first language no block covers; in regions of another language & "
                    f"\\multicolumn{{6}}{{l}}{{{pct1(b3['country_level_share'])}; {pct1(b3['regional_share_v66'])}}} \\\\")
    cap = ("The lexicon reading (exploratory). Upper part: share of analysed cells whose nearest name is water-side, by "
           "the continent of the nearest settlement (Americas: North and South together). Middle: the forms that supply "
           "the most water-side cells, with the language block and the share of water-side cells (a name can match "
           "several forms). Lower part: agreement with derivations recorded on Wikipedia, with 95\\,\\% bootstrap "
           "intervals. Precision: share of names read as water-side whose recorded derivation is water-related; recall: "
           "share of water-related derivations read as water-side, the two strata of the sample weighted by their size. "
           "Coder 2 is a language model of another developer; coder 1 is the model that compiled the lexicon and is "
           "therefore not an independent validator. The sample is restricted to articles that explain the name and "
           "almost excludes Pakistan, Bangladesh and Cambodia.")
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table*}[tp]\n"
           f"\\caption{{{cap}}}\\label{{tab:reading}}\n\\footnotesize\n\\begin{{tabular}}{{p{{4.4cm}}rrrrrr}}\n\\tophline\n"
           " & Europe & Asia & Africa & Oceania & Americas & All \\\\\n\\middlehline\n"
           + "\n".join(rows) + "\n\\bottomhline\n\\end{tabular}\n\\end{table*}\n")
    (HERE / "tab_reading.tex").write_text(tex, encoding="utf-8")


def table_twelve_checks(v69: dict | None) -> str:
    """補足の表（B-12 の B）: 規則の動作特性、基準のモデル、定数の予測、帰無の外れた回、言語、包含、GFM の窓、種別、
    市街地、上限なしの標本の構成。"""
    v = v69 or {}
    rows: list[str] = []
    b1 = (v.get("B1") or {}).get("share_equivalent") or {}
    if b1:
        rows.append(r"\multicolumn{3}{l}{\emph{Share of resampled samples declared equivalent (item B1)}} \\")
        for sh, lab in (("constant", "shift uniform over cells"), ("proportional", "shift proportional to the baseline gain")):
            for r0 in ("-0.05", "+0.00", "+0.05"):
                x = b1[sh][r0]
                rows.append(f"True $R={r0.replace('+0.00', '0')}$, {lab}: protocol's rule; ratio rule & "
                            f"{_pct0(x['fixed'])}; {_pct0(x['ratio'])} & \\\\")
    b10 = v.get("B10") or {}
    if b10:
        bs = b10["brier"]
        rows.append(r"\multicolumn{3}{l}{\emph{Reference model on the primary split, corrected code (item B10)}} \\")
        rows.append(f"Brier score: M1s; M1c; M1c without the area mean of the character model's score & "
                    f"{bs['M1s']:.5f}; {bs['M1c']:.5f}; {bs['M1c_no_nlearn']:.5f} & \\\\")
        rc = b10["residual_corr"]
        rows.append(f"Correlation of sibling AOI residuals within 50\\,km (Fisher 95\\,\\%) & "
                    f"{rc['r']:.2f} ({rc['ci95_fisher'][0]:.2f} to {rc['ci95_fisher'][1]:.2f}) & {rc['n_pairs']} pairs \\\\")
    b17 = v.get("B17") or {}
    if b17:
        rows.append(r"\multicolumn{3}{l}{\emph{Position model against a constant forecast (item B17)}} \\")
        for k, lab in (("primary", "location groups held out"), ("country", "countries held out"),
                       ("time", "floods from 2024 held out"), ("llmcut", "floods after the language model's cut-off"),
                       ("activation_frozen", "activations kept whole")):
            x = b17.get(k) or {}
            if x.get("brier_M0"):
                rows.append(f"Brier score of M0; of the best constant, {lab} & {x['brier_M0']:.4f}; "
                            f"{x['brier_best_constant']:.4f} & {'worse' if x['M0_worse_than_best_constant'] else 'better'} \\\\")
    b5 = (v.get("B5") or {}).get("outlying") or []
    if b5:
        rows.append(r"\multicolumn{3}{l}{\emph{Outlying null simulations of the protocol (item B5)}} \\")
        for cc, gname in sorted({(r["top_unit_country"], r["top_unit"]) for r in b5}):
            rr = [r for r in b5 if r["top_unit"] == gname]
            rows.append(f"Replicates led by location group {gname} ({cc}); its share of $R$ & {len(rr)} & "
                        f"{min(r['top_unit_share_of_R'] for r in rr):.2f} to {max(r['top_unit_share_of_R'] for r in rr):.2f} \\\\")
    b4 = v.get("B4") or {}
    if b4:
        rows.append(r"\multicolumn{3}{l}{\emph{Inclusion criteria applied in order (item B4)}} \\")
        for s_ in b4["steps"]:
            rows.append(f"{s_['criterion']} & {s_['excluded']} excluded & {s_['remaining']} remaining \\\\")
    b13 = v.get("B13") or {}
    if b13:
        rows.append(r"\multicolumn{3}{l}{\emph{Flood typing (item B13) and the GFM window (item B12)}} \\")
        rows.append(f"AOIs with untyped flooding; counted riverine by the inland rule & {b13['aoi_with_untyped']}; "
                    f"{b13['riverine_by_inland_rule']} & with the ice-jam AOIs \\\\")
        rows.append(f"AOIs with typed flooding of another kind; descriptions & {b13['aoi_with_typed_other']} & "
                    + "; ".join(f"{d} ({n})" for d, n in b13["other_descriptions_top"][:3]) + " \\\\")
        ts = b13["trace_share_by_sensor"]
        rows.append("Flood traces as a share of the flooded area: optical; Sentinel-1; other radar & "
                    f"{pct1(ts.get('optical') or 0)}; {pct1(ts.get('s1') or 0)}; {pct1(ts.get('sar') or 0)} & \\\\")
        co = b13["coast_km"]
        rows.append(f"Distance of the AOI centroid from the coast: minimum; 10th percentile; median (km) & "
                    f"{co['min']:.1f}; {co['q10']:.1f}; {co['median']:.0f} & {co['n_le_20km']} AOIs within 20\\,km \\\\")
    b12w = v.get("B12") or {}
    if b12w:
        rows.append(f"AOIs with GFM images: with one on an EMS image day; only on adjacent days & "
                    f"{b12w['aoi_with_an_image_on_an_ems_day']}; {b12w['aoi_with_images_only_on_adjacent_days']} & "
                    f"of {b12w['n_aoi']} \\\\")
    b18 = v.get("B18") or {}
    if b18.get("with_built10"):
        rows.append(r"\multicolumn{3}{l}{\emph{Labels near settlements and built-up land at 10\,m (item B18)}} \\")
        for k, lab in (("band_500m_1km", "ground 0.5--1\\,km against within 500\\,m"),
                       ("band_gt1km", "ground beyond 1\\,km against within 500\\,m"),
                       ("wetside:le500m", "water-side names within 500\\,m (cluster-robust)")):
            a_, b_ = b18["without"][k], b18["with_built10"][k]
            rows.append(f"Odds ratio, {lab}: without; with the built-up share & {orv(a_['odds_ratio'])} "
                        f"({orv(a_['ci95'][0])}--{orv(a_['ci95'][1])}); {orv(b_['odds_ratio'])} "
                        f"({orv(b_['ci95'][0])}--{orv(b_['ci95'][1])}) & \\\\")
        pp_ = b18.get("pairs") or {}
        if pp_.get("with"):
            rows.append(f"\\quad water-side, pairs bootstrap 95\\,\\%: without; with & "
                        f"{orv(pp_['without']['pairs_ci95'][0])}--{orv(pp_['without']['pairs_ci95'][1])}; "
                        f"{orv(pp_['with']['pairs_ci95'][0])}--{orv(pp_['with']['pairs_ci95'][1])} & \\\\")
    b2 = v.get("B2") or {}
    if b2.get("by_country"):
        rows.append(r"\multicolumn{3}{l}{\emph{Every evaluable cell: where the flooded cells within 500\,m of water-side "
                    r"names are (item B2)}} \\")
        for r in b2["by_country"][:6]:
            rows.append(f"{_CC_NAMES.get(r['cc'], r['cc'])[:1].upper() + _CC_NAMES.get(r['cc'], r['cc'])[1:]} & "
                        f"{_thousands(r['flooded_near'])} ({pct1(r['share_of_flooded_near'])}) & {r['aois']} AOIs \\\\")
        for r in (b2.get("by_activation") or [])[:2]:
            rows.append(f"\\quad activation {r['activation']} & {_thousands(r['flooded_near'])} "
                        f"({pct1(r['share_of_flooded_near'])}) & {r['aois']} AOIs \\\\")
        nb = b2.get("without_pk_bd") or {}
        if nb.get("pairs_ci95"):
            rows.append(f"Odds ratio without Pakistan and Bangladesh (pairs 95\\,\\%) & {orv(nb['odds_ratio'])} "
                        f"({orv(nb['pairs_ci95'][0])}--{orv(nb['pairs_ci95'][1])}) & "
                        f"{_thousands(nb['n_flooded_near'])} flooded cells \\\\")
    b11 = v.get("B11") or {}
    if b11.get("R") is not None:
        rows.append(r"\multicolumn{3}{l}{\emph{Penalty on the name features and the composition only (item B11)}} \\")
        cs = b11["C_chosen"]
        rows.append(f"$R$; chosen $C$ for M1c and M3n (range over folds) & {rv(b11['R'])} & "
                    f"{min(cs['M1c']):.0e}--{max(cs['M1c']):.0e}; {min(cs['M3n']):.0e}--{max(cs['M3n']):.0e} \\\\")
    cap = ("Further checks of amendment 12 (exploratory). B1: shares of resampled samples in which each decision rule "
           "declares equivalence when the true $R$ is set by shifting the names increment (uniformly over cells, or in "
           "proportion to each location group's baseline gain); the protocol's rule divides the interval of the change "
           "in Brier score by the point baseline gain, the ratio rule uses the interval of $R$. B10: the reference model "
           "refitted from the second-stage design of the corrected run. B17: held-out Brier scores of the position model "
           "and of the best constant forecast (the test cells' flooded share). B5: the thirteen outlying null "
           "replicates of the protocol's simulations. B4, B12, B13: inclusion, GFM images and flood typing. B18: the "
           "model of Sect.~6.2 with and without the share of built-up land at 10\\,m (ESA WorldCover 2021). B2: countries "
           "of the flooded cells within 500\\,m of water-side names in the sample of every evaluable cell. B11: the second "
           "stage with the penalty chosen by grouped inner cross-validation for the name features and the composition "
           "only, on a grid widened until the optimum is interior.")
    return _tbl("tab:twelve-checks", cap, "p{9.4cm}p{3.6cm}l", r"Quantity & Value & Note", rows)


def fig_spec(v69: dict | None) -> None:
    """補足の図 S6（B-12 の B19、追記 2・3）: 集落の近くの水系側の OR の全ての版。色はラベルと評価の範囲。"""
    rows = ((v69 or {}).get("B19") or {}).get("rows") or []
    if not rows:
        return
    d = pd.DataFrame(rows)
    grp = {("GFM", "gfm_cells"): (0, "#1f77b4", "GFM labels"),
           ("EMS", "gfm_cells"): (1, "#d95f02", "EMS labels, cells GFM can assess"),
           ("EMS", "all"): (2, "#444444", "EMS labels, all cells")}
    d["g"] = [grp[(a, b)][0] for a, b in zip(d["labels"], d["domain"], strict=True)]
    d = d.sort_values(["g", "estimate"], ascending=[False, True]).reset_index(drop=True)
    planned = (d["file"].str.startswith("v41") & (d["path"] == "odds_gfm/wetside_adjusted_M1/terms/wetside:le500m"))
    chosen = d["path"].isin(["B2/ems_all", "B2/gfm"])
    fig, ax = plt.subplots(figsize=(7.2, 0.15 * len(d) + 1.1))
    for i, r in d.iterrows():
        _, col, _ = grp[(r["labels"], r["domain"])]
        mk = "D" if planned[i] else ("s" if chosen[i] else "o")
        ax.plot([r["ci95"][0], r["ci95"][1]], [i, i], color=col, lw=1.0)
        ax.plot([r["estimate"]], [i], marker=mk, ms=4.5 if mk != "o" else 3.5, color=col,
                mfc="white" if r["uncapped"] else col)
    ax.axvline(1.0, color="0.6", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(0.18, 2.4)
    ax.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.25", "0.5", "1", "2"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([f"am. {r['amendment']}: {r['file'].split('_')[0]} {r['path'][:52]}" for _, r in d.iterrows()],
                       fontsize=5.2)
    ax.set_ylim(-0.7, len(d) - 0.3)
    ax.set_xlabel("odds ratio of flooding within 500 m of water-side names (95 % interval)", fontsize=7)
    ax.tick_params(axis="x", labelsize=7)
    hs = [plt.Line2D([], [], color=c, marker="o", lw=1.0, ms=3.5, label=lab) for _, c, lab in grp.values()]
    hs += [plt.Line2D([], [], color="0.3", marker="D", lw=0, ms=4.5, label="planned (amendment 3)"),
           plt.Line2D([], [], color="0.3", marker="s", lw=0, ms=4.5, label="chosen after both were seen (amendment 11)"),
           plt.Line2D([], [], color="0.3", marker="o", mfc="white", lw=0, ms=3.5, label="every evaluable cell (no cap)")]
    ax.legend(handles=hs, fontsize=5.5, loc="lower left", frameon=False)
    fig.tight_layout()
    fig.savefig(HERE / "fig_spec.pdf")
    plt.close(fig)


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
    fig_plasmode(v52, v56, v59, _json("v64_b9_aggregate.json"), _json("v67_b10_aggregate.json"))
    fig_near_more(v55, v57, _json("v61_leaveout.json"))
    fig_near(_json("v68_b11_local.json"), v57)
    t = pd.read_parquet(DATA / "verification" / "v41_table.parquet",
                        columns=["aoi", "y", "ix", "iy", "p_M1c", "p_M3n", "P_d1", "lex_wetside_n1"])
    gs = g._json(DATA / "global" / "gfm_summary.json")
    ex = pick_example(t, gs) if gs else None
    ex2 = pick_example_near(t, gs) if gs else None
    m |= (fig_example([ex] + ([ex2] if ex2 and ex2 != ex else []), t) if ex else
          {"ExampleAoiName": missing("example AOI")})
    m |= timeline(m) | gfm_items()
    table_amendments(m)
    write_new_tables(v41, v42, v44, v45, v47)
    table_ceiling()
    table_bounds()
    (HERE / "app_table_gfm_compact.tex").write_text(table_gfm_compact(v41, v47, v50), encoding="utf-8")
    b6_supp_tables(v48, v50, v51, v52, v53, v54)
    (HERE / "app_table_posthoc_models.tex").write_text(table_posthoc_models(), encoding="utf-8")
    (HERE / "app_table_checks.tex").write_text(table_checks(checks), encoding="utf-8")
    supp_tables_b7(post["v37"], v43, v44, v55, v56, v52)
    v60, v61 = _json("v60_where_near.json"), _json("v61_leaveout.json")
    supp_tables_b8(v57, v57b, v59, v41, v60, v61)
    supp_tables_b9(_json("v63_b9_local.json"), _json("v64_b9_aggregate.json"), v57)
    v66, v67 = _json("v66_b10_local.json"), _json("v67_b10_aggregate.json")
    m |= b10_macros(v66, v67, v55, v56, v57, v44, v48, _json("v61_leaveout.json"), _json("v64_b9_aggregate.json"), m)
    supp_tables_b10(v66, v67, v44)
    v68 = _json("v68_b11_local.json")
    m |= b11_macros(v68)
    m |= b11_text_macros(post["v36"], v54, v57, v57b)
    table_eleven_headings()
    v69 = _json("v69_b12_local.json")
    supp_tables_b11(v68, v69)
    m |= b12_text_macros(v69, v68, v57, v41, v50, v54, v44, v43, post["v36"])
    table_reading(v54, v44, v57b, v69, g._json(DATA / "runs" / "world_v2" / "metrics_primary.json"))
    (HERE / "app_table_twelve_checks.tex").write_text(table_twelve_checks(v69), encoding="utf-8")
    fig_spec(v69)
    m |= _orp("NoPkOR", (v61 or {}).get("without_pakistan")) | _orp("NoSAsiaOR", (v61 or {}).get("without_pakistan_bangladesh"))
    if v60:  # 水系側の名の近くの浸水セルの散らばり（編集者の求め、記述）
        m |= {"WhereFlooded": str(v60["flooded_near_total"]),
              "WhereNCountries": str(v60["n_countries_with_flooded_near"]),
              "WhereNActivations": str(v60["n_activations_with_flooded_near"]),
              "WhereTopFivePct": pct1(v60["share_top5_countries"]),
              "WhereTopActPct": pct1(v60["largest_activation_share"]),
              "WhereTopCountries": _country_list([r["cc"] for r in v60["by_country"][:5]]),
              # パキスタンとバングラデシュの割合。討論の「地名学者が見ていない、特にこの 2 か国」の文に使う
              "WherePkBdPct": pct1(sum(r["share_of_flooded_near"] for r in v60["by_country"]
                                       if r["cc"] in {"BD", "PK"}))}
    else:
        m |= {k: missing("v60") for k in ("WhereFlooded", "WhereNCountries", "WhereNActivations", "WhereTopFivePct",
                                          "WhereTopActPct", "WhereTopCountries", "WherePkBdPct")}
    tables = copernicus_tables()
    m |= supplement()
    write_numbers(m)
    print(f"macros: {len(m)}; tables: {len(tables)}; missing: {len(g.MISSING)}")
    for w in sorted(set(g.MISSING)):
        print("  missing:", w)


if __name__ == "__main__":
    main()
