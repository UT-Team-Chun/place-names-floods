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
                "OracleVerdict": o["decision"]}
    else:
        out |= {k: missing("oracle") for k in ("OracleR", "OracleRLo", "OracleRHi", "OracleVerdict")}
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
def fig_main(primary: dict, extra: dict, v36: dict | None, v41: dict | None, v42: dict | None,
             v47: dict | None) -> None:
    rename = {"learned only (N-learn)": "character model only (N-learn)",
              "placebo: names permuted": "placebo: names reassigned",
              "lexicon and learned score, same floods": "lexicon and character model, same floods"}
    def cap(t: str) -> str:
        return t[:1].upper() + t[1:]

    rows = [((lab[:len(lab) - len(lab.lstrip())] + cap(rename.get(lab.strip(), lab.strip()))), r, c, st)
            for lab, r, c, st in g._forest_rows(primary, extra, v36)
            if not lab.startswith(("Total incl.", "Flexible model"))]
    v47 = v47 or {}
    lf = (v47.get("lexfix") or {}).get("primary")
    if isinstance(lf, dict):
        rows.insert(1, ("   Same, Norwegian names read correctly", lf["R"], lf["ci95_R"], "main"))
    post: list[tuple[str, float, list, str]] = []
    if v41:
        post.append(("Independent labels (GFM): rescored", v41["C1_gfm"]["R"], v41["C1_gfm"]["ci95_R"], "post"))
    c2 = v47.get("gfm_c2")
    if isinstance(c2, dict):
        if c2["decision"] == "undefined-denominator":
            post.append(("Independent labels (GFM): refitted", None, None, "post"))
        else:
            post.append(("Independent labels (GFM): refitted", c2["R"], c2["ci95_R"], "post"))
    for key, lab in (("notrace", "Flood traces counted as not flooded"),
                     ("nonwater", "Flooded share over land only")):
        lv = (v47.get("label_variants") or {}).get(key)
        if isinstance(lv, dict):
            post.append((lab, lv["R"], lv["ci95_R"], "post"))
    for run, lab in (("world_v2_actgrp", "Flood events held out whole"),):
        met = g._json(DATA / "runs" / run / "metrics_primary.json")
        if met:
            post.append((lab, met["primary"]["R"], met["primary"]["ci95_R"], "post"))
    if v42:
        s = v42["settler_excluded"]
        c95 = [v / s["yardstick_M0_M1"] for v in s["ci95_delta"]]  # 固定マージンの形（ΔBS の区間 / 物差しの点推定）
        post.append(("Without the Americas and Oceania", s["R"], c95, "post"))
    rows_all = rows + post
    fs = (v47.get("fold_seeds") or {}).get("fixed_code") or {}
    pm = v47.get("permutations")
    g.SOURCE["Fig4"] = [{"comparison": n.strip(), "R": r, "ci95_low": c[0] if c else None,
                         "ci95_high": c[1] if c else None, "group": s} for n, r, c, s in rows_all]
    n = len(rows_all) + 2
    fig, ax = plt.subplots(figsize=(6.6, 0.32 * n + 1.3))
    ax.axvspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axvline(0, color="0.4", lw=0.8)
    ys = np.arange(n)[::-1]
    style = {"main": ("#1f6fb4", "-", True), "secondary": ("0.45", "-", False),
             "undefined": ("0.45", "--", False), "post": ("#8c5a2b", "-", False)}
    rows_y = ys[:-2]  # 下の 2 段（y = 1, 0）は fold の割り当てと置換の分布
    for y, (_, r, c, s) in zip(rows_y, rows_all, strict=True):
        col, ls, filled = style[s]
        if r is None:  # R が未定義（物差しの 95 % 区間が 0 を含む）
            ax.text(0.0, y, "R undefined: terrain did not improve on position", ha="center", va="center",
                    fontsize=6, color=col, style="italic")
            continue
        if c is not None:
            ax.plot(c, [y, y], color=col, lw=1.5, ls=ls)
        ax.plot([r], [y], "o", ms=4.2, color=col, mfc=col if filled else "white")
    # 下の 2 段: fold の割り当て（点）と、置換の分布の 95 % 範囲（灰の帯）に実の名前（青の点）
    y_fs, y_pm = 1, 0
    if fs.get("R"):
        rs = np.array(list(fs["R"].values()), dtype=float)
        ax.plot(rs, np.full(len(rs), y_fs), "o", ms=3.2, color="#8c5a2b", mfc="white", alpha=0.9)
    if isinstance(pm, dict):
        ax.plot([pm["R_perm_q025"], pm["R_perm_q975"]], [y_pm, y_pm], color="0.7", lw=5, solid_capstyle="butt")
        ax.plot([pm["R_real"]], [y_pm], "o", ms=4.2, color="#8c5a2b", mfc="#8c5a2b")
    n_perm = pm["n_perms"] if isinstance(pm, dict) else 0
    labels = [r[0] for r in rows_all] + [f"{fs.get('n_available', 0)} assignments of places to folds",
                                         f"Names reassigned {n_perm} times (grey), real names (dot)"]
    ax.set_yticks(list(rows_y) + [y_fs, y_pm])
    ax.set_yticklabels(labels, fontsize=6.6)
    for y, (nm, *_) in zip(rows_y, rows_all, strict=True):
        if nm.startswith(("Countries held out", "Independent labels (GFM): rescored")):
            ax.axhline(y + 0.5, color="0.8", lw=0.6)
    ax.axhline(y_fs + 0.5, color="0.8", lw=0.6)
    ax.set_xlabel("R (improvement from names / improvement from terrain)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(-0.8, ys[0] + 0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(HERE / "fig_main.pdf")
    plt.close(fig)


# %% 図 5（本文）: 名前が記述するものと、集落の近くの浸水
def fig_names(v37: dict, v39: dict | None, v41: dict | None, v42: dict | None) -> None:
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
    # c: 500 m 以内の浸水のオッズ比（地形で調整）: EMS 全 AOI（v39）、GFM で評価できる AOI の EMS と GFM（v41）
    ax = axs[2]
    sets = []
    if v39:
        sets.append(("EMS, all AOIs", v39, blue, 0.22, True))
    if v41:
        sets.append(("EMS, radar AOIs", v41["odds_ems"], grey, 0.0, False))
        sets.append(("GFM, radar AOIs", v41["odds_gfm"], green, -0.22, True))
    for lab, src, col, dy, filled in sets:
        first = True
        for y, c in zip(ys, classes, strict=True):
            key = "wetside" if c == "wetside" else "categories"
            t = src[f"{key}_adjusted_M1"]["terms"].get(f"{c}:le500m")
            if not t or not np.isfinite(t["odds_ratio"]) or t["odds_ratio"] <= 0:
                continue
            lo, hi = (max(v, 0.05) for v in t["ci95_odds_ratio"])
            ax.plot([lo, min(hi, 20)], [y + dy, y + dy], color=col, lw=1.3)
            ax.plot([t["odds_ratio"]], [y + dy], "o", ms=3.8, color=col, mfc=col if filled else "white",
                    label=lab if first else None)
            first = False
    ax.set_xscale("log")
    ax.axvline(1, color="0.4", lw=0.8)
    ax.set_xlim(0.1, 10)
    ax.xaxis.set_major_locator(FixedLocator([0.1, 0.25, 0.5, 1, 2, 4, 10]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.1", "0.25", "0.5", "1", "2", "4", "10"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(ys)
    ax.set_yticklabels([g.CLASS_LABELS[c] for c in classes], fontsize=7)
    ax.set_xlabel("odds ratio of flooding within 500 m of the settlement,\ngiven the terrain model", fontsize=7)
    ax.legend(fontsize=5.8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=3,
              handletextpad=0.3, columnspacing=1.0)
    # d: 検定の反転（B-5 の 4）: δ ごとの模擬 ΔBS の 5 パーセンタイルと観測値
    ax = axs[3]
    if v42:
        for key, lab, col, mk in (("inversion_area_wide", "effect on all ground nearest the name", blue, "o"),
                                  ("inversion_local_500m", "effect within 500 m only", brown, "s")):
            inv = v42[key]
            grid = sorted((float(k), v) for k, v in inv["grid"].items())
            dx = np.array([k for k, _ in grid])
            q05 = 1e4 * np.array([v["q05"] for _, v in grid])
            mu = 1e4 * np.array([v["mean"] for _, v in grid])
            ax.plot(np.exp(dx), q05, mk + "-", ms=3.5, color=col, lw=1.1, label=f"{lab}, 5th percentile")
            ax.plot(np.exp(dx), mu, mk + ":", ms=2.5, color=col, lw=0.8, mfc="white")
            ax.axvline(inv["or_upper_q05"], color=col, lw=0.7, ls="--")
        obs = 1e4 * v42["inversion_area_wide"]["observed_delta_bs"]
        ax.axhline(obs, color="0.2", lw=0.9)
        ax.text(1.02, obs, "observed", fontsize=6, va="bottom", color="0.2")
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator([1, 1.5, 2, 3, 5]))
        ax.xaxis.set_major_formatter(FixedFormatter(["1", "1.5", "2", "3", "5"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlim(0.95, 5.5)
        ax.set_ylim(-0.6, 2.5)
        ax.set_xlabel("simulated odds ratio for ground\nnearest a water-side name", fontsize=7)
        ax.set_ylabel("improvement from names, ΔBS (×10$^{-4}$)\n(solid: 5th percentile, dotted: mean)",
                      fontsize=7)
        ax.legend(fontsize=5.8, frameon=False, loc="upper left")
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
                         columns=["gid", "asciiname", "cc", "lat", "lon"]).drop_duplicates("gid")
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
                         columns=["gid", "asciiname", "cc", "lat", "lon"]).drop_duplicates("gid")
    st = settlements_in(meta, gn)
    lex = g_lexicon()
    st = st.assign(wet=[any(lex.match(r.asciiname, r.cc)[c] for c in ("water", "wet", "riverine", "crossing", "low"))
                        for r in st.itertuples()])
    aspect = ny / nx
    fig, axs = plt.subplots(2, 2, figsize=(7.0, min(9.0, 0.5 + 2 * 3.0 * aspect + 0.6)), layout="constrained")
    from matplotlib.colors import ListedColormap
    lab_cmap = ListedColormap(["0.82", "#1f6fb4"])
    titles = ["(a) EMS labels", "(b) GFM labels", "(c) height above nearest drainage",
              "(d) change in flood probability from names"]
    for ax, arr in ((axs[0, 0], ems), (axs[0, 1], gf)):
        ax.imshow(arr, cmap=lab_cmap, vmin=0, vmax=1, interpolation="nearest")
    hmax = float(np.nanquantile(hand, 0.98)) or 1.0
    im = axs[1, 0].imshow(np.clip(hand, 0, hmax), cmap="YlOrBr_r", interpolation="nearest")
    cb = fig.colorbar(im, ax=axs[1, 0], shrink=0.7)
    cb.set_label(f"HAND (m; clipped at the 98th percentile, {hmax:.0f} m)", fontsize=6.5)
    lim = float(np.nanquantile(np.abs(diff), 0.99)) or 1e-3
    im = axs[1, 1].imshow(diff, cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    cb = fig.colorbar(im, ax=axs[1, 1], shrink=0.7)
    cb.set_label("M3n − M1c (probability)", fontsize=7)
    for ax in axs.ravel():
        ax.plot(st.loc[~st["wet"], "col"], st.loc[~st["wet"], "row"], "k^", ms=3, mfc="white", mew=0.6)
        ax.plot(st.loc[st["wet"], "col"], st.loc[st["wet"], "row"], "^", ms=3.6, color="#e66101", mew=0.6)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax, title in zip(axs.ravel(), titles, strict=True):
        ax.set_title(title, fontsize=8, loc="left")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(color="#1f6fb4", label="flooded"), Patch(color="0.82", label="not flooded"),
               Line2D([], [], marker="^", ls="", color="#e66101", label="settlement, water-side name"),
               Line2D([], [], marker="^", ls="", color="k", mfc="white", label="other settlement")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    fig.savefig(HERE / "fig_example.pdf")
    plt.close(fig)
    act = meta["activation"]
    country = ", ".join(act.get("countries") or [])
    tag = aoi.replace("_", "\\_")
    name = act.get("name", "")
    where = name if country and country in name else f"{name}, {country}"
    return {"ExampleAoiName": f"{tag}: {where}, {act.get('activation_time', '')[:4]}",
            "ExampleNSettlements": str(len(st)), "ExampleNWetSettlements": str(int(st['wet'].sum()))}


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
         "BFive": "B-5_blind_review_analyses.md"}
CONFIRM_IMAGE_NS = 1792961983094088000  # 確証解析のイメージ（amendment B-1 に記録したタグ）


def timeline(m: dict[str, str]) -> dict[str, str]:
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
    blind = _first_commit_time("docs/research/reviews/2026-10-26_paper2_blind1_editor.md")
    rows = [("EMS product inventory retrieved", m.get("InvRetrievedUTC", "")),
            ("GeoNames gazetteer retrieved", m.get("GeoNamesUTC", "")),
            ("Power analysis on synthetic labels aggregated", (pw.get("created_utc") or "")[:16] + "Z"),
            ("Protocol content fixed (commit " + m.get("FreezeContentCommit", "") + ")", content or ""),
            ("Protocol marked frozen (commit " + m.get("FreezeMarkCommit", "") + ", tag prereg-v2-frozen)",
             mark or ""),
            ("Container image of the confirmatory analysis built", img),
            ("Confirmatory analysis aggregated (primary split)", (conf.get("created_utc") or "")[:16] + "Z"),
            ("Amendment B-1 (sensitivity analysis 6 corrected)", out["AmendTimeBOne"]),
            ("Amendment B-2 (post hoc analyses, internal review)", out["AmendTimeBTwo"]),
            ("Amendment B-3 (independent GFM labels)", out["AmendTimeBThree"]),
            ("GFM labels built", (gs.get("built_utc") or "")[:16] + "Z"),
            ("Review by readers of the manuscript recorded", blind or ""),
            ("Amendment B-4 (Norwegian names; rerun)", out["AmendTimeBFour"]),
            ("Amendment B-5 (post hoc analyses after that review)", out["AmendTimeBFive"])]
    body = "\n".join(f"{a} & {b} \\\\" for a, b in rows)
    tex = ("% generated by figures/make_figures.py — do not edit\n\\begin{table}[t]\n"
           "\\caption{Timeline of the study (UTC). Times are those of the version-controlled repository, of the run "
           "records and of the container image (whose tag is its build time).}\\label{tab:timeline}\n\\small\n"
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
    names = {"wetland": "wetland", "island": "island and river island", "spring_well": "spring and well",
             "crossing": "crossing", "dike_reclamation": "dike and polder", "lowland": "lowland (reference)"}
    for key, lab in names.items():
        gb = v42["groups"]["by_group"][key]
        e1, c1 = _orc(gb.get("or_le500m_adj"))
        e2, c2 = _orc(gb.get("or_all_adj"))
        hm = gb.get("pc_hand_mean_all")
        hand = (f"${hm['estimate']:+.1f}$ (${hm['ci95'][0]:+.1f}$ to ${hm['ci95'][1]:+.1f}$)"
                if isinstance(hm, dict) else "---")
        ph = gb.get("p_le500m_holm")
        rows.append(f"{lab} & {gb['n_settlements']} & {_n(gb['n_cells_le500m'])} & {gb['n_flooded_le500m']} & "
                    f"{e1} ({c1}) & {e2} ({c2}) & {hand} & {'---' if ph is None else f'{ph:.2f}'} \\\\")
    cap = ("Element groups of the lexicon (post hoc, amendment B-5): odds ratios of flooding given the held-out "
           "terrain prediction within 500\\,m of the settlement and over all cells, the difference in mean HAND "
           "(m) between the cells of the group and other cells of the same AOI, and Holm-adjusted $p$ values. No "
           "cell within 500\\,m of a dike or polder name was flooded, so that odds ratio cannot be estimated.")
    return _tbl("tab:groups", cap, "lrrrcccc",
                "Group & Settlements & Cells $\\leq500$\\,m & Flooded & OR $\\leq500$\\,m & OR all cells & HAND & "
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


def main() -> None:
    m, primary, extra, post = base_macros()
    v41, v42, v43 = _json("v41_gfm_independent_labels.json"), _json("v42_blind_review_analyses.json"), \
        _json("v43_wikidata_etymology.json")
    v44, v45, v46, v47 = (_json("v44_wikipedia_etymology.json"), _json("v45_gfm_masking_check.json"),
                          _json("v46_norway_impact.json"), _json("v47_b5_cluster.json"))
    m |= freeze_record() | norway_macros(v46) | gfm_macros(v41, v45, v47) | b5_macros(v42, v47)
    m |= etym_macros(v43, v44)
    if primary:
        fig_main(primary, extra, post["v36"], v41, v42, v47)
    if post["v37"]:
        fig_names(post["v37"], post["v39"], v41, v42)
    t = pd.read_parquet(DATA / "verification" / "v41_table.parquet",
                        columns=["aoi", "y", "ix", "iy", "p_M1c", "p_M3n"])
    gs = g._json(DATA / "global" / "gfm_summary.json")
    ex = pick_example(t, gs) if gs else None
    m |= fig_example(ex, t) if ex else {"ExampleAoiName": missing("example AOI")}
    m |= timeline(m) | gfm_items()
    write_new_tables(v41, v42, v44, v45, v47)
    tables = copernicus_tables()
    write_numbers(m)
    print(f"macros: {len(m)}; tables: {len(tables)}; missing: {len(g.MISSING)}")
    for w in sorted(set(g.MISSING)):
        print("  missing:", w)


if __name__ == "__main__":
    main()
