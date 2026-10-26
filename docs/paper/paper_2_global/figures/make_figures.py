# %% Study 2（全球検定）の図表・数値マクロ生成。数値は全て成果物由来（手打ち禁止、paper_1 と同じ約束）。
# 使い方（docs/paper/paper_2_global/ から）: make figures
#
# 入力（無いものは \missing{...} の印を出してコンパイルは通す。投稿前に印が残っていないことを
# `make check` で確認する）:
#   docs/data/ems_inventory_20261025.json      EMS 洪水台帳（凍結スナップショット）
#   data/global/incidence.json                 ラベル単独の発生表（world_incidence.py、NAS から取得）
#   data/global/geonames/allCountries.zip.sha256 / retrieved_utc.txt
#   docs/protocol/v2/lexicon_hydro_v1.yaml     N-lex 辞書
#   docs/protocol/preregistration_v2.md        凍結行（FROZEN <UTC> / commit <sha>）
#   data/global/nread/nread_meta.json          N-read 推論の来歴（Miyabi run から複写）
#   data/runs/<power-null|power-sesoi>/power_summary.json
#   data/runs/<eval>/summary_<split>.json      確証解析（凍結後）
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from collections import Counter
from itertools import pairwise
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "0")  # PDF に作成日時を埋めない（再生成で差分を出さない）
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[4]
OUTDIR = Path(__file__).resolve().parent
DATA = REPO / "data"
MISSING: list[str] = []


def _json(path: Path) -> dict | None:
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def missing(what: str) -> str:
    MISSING.append(what)
    return f"\\missing{{{what}}}"


def pct(x: float) -> str:
    """百分率。10% 未満は小数 1 桁（0.6% を 1% と丸めない）。"""
    v = 100 * x
    return f"{v:.1f}\\%" if abs(v) < 10 else f"{v:.0f}\\%"


def sci(x: float, nd: int = 1) -> str:
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    if -1 <= e <= 2:  # 0.01 未満と 0.01–0.1 を同じ指数表記にそろえる（区間の両端で書式が混ざらないように）
        return f"{x:.{max(nd + 1 - e, 0)}f}"
    return f"\\ensuremath{{{x / 10**e:.{nd}f}\\times10^{{{e}}}}}"


def inventory_macros(inv: dict) -> dict[str, str]:
    acts = inv["activations"]
    with_del = [a for a in acts if a.get("del_products")]
    return {
        "NActFlood": str(inv["n_flood_activations"]),
        "NActDel": str(len(with_del)),
        "InvRetrievedUTC": inv["retrieved_utc"][:10],
    }


def _kish(inc: dict) -> float:
    """場所グループの実効数（Kish: (Σw)²/Σw²、w = グループに配分されたセル数）。"""
    w: Counter = Counter()
    for r in inc["aois"]:
        if r["included"]:
            w[r["group"]] += r["quota"]
    return sum(w.values()) ** 2 / sum(v * v for v in w.values())


# EMS の国名 → GeoNames countryInfo の国名（表記揺れ 4 件、2026-10-25 の発生表で確認）
COUNTRY_ALIASES = {"Myanmar/Burma": "Myanmar", "Netherlands": "The Netherlands",
                   "Timor-Leste": "Timor Leste", "Viet Nam": "Vietnam"}
CONTINENT_NAMES = {"EU": "Europe", "AF": "Africa", "AS": "Asia", "OC": "Oceania",
                   "NA": "NorthAmerica", "SA": "SouthAmerica"}


def continent_of_name() -> dict[str, str]:
    out = {}
    path = DATA / "global" / "geonames" / "countryInfo.txt"
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            f = line.split("\t")
            out[f[4]] = f[8]
    return out


def sample_macros(inc: dict) -> dict[str, str]:
    """標本の記述（ラベル単独の発生表）: 大陸別の採用 AOI 数と除外理由。"""
    cont = continent_of_name()
    rows = [r for r in inc["aois"] if r["included"]]
    by = Counter(cont.get(COUNTRY_ALIASES.get(r["countries"][0], r["countries"][0]), "??")
                 if r["countries"] else "??" for r in rows)
    out = {f"NAoi{name}": str(by.get(code, 0)) for code, name in CONTINENT_NAMES.items()}
    out["NAoiUnknownContinent"] = str(by.get("??", 0))
    n_cont = sum(1 for code in CONTINENT_NAMES if by.get(code, 0) > 0)
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven"}
    out["NContinentsWord"] = words.get(n_cont, str(n_cont))
    reasons = inc["summary"]["reasons"]
    out["NExclWet"] = str(reasons.get("wet<20", 0))
    out["NExclEval"] = str(reasons.get("eval<400", 0))
    return out


def incidence_macros(inc: dict) -> dict[str, str]:
    s = inc["summary"]
    included = [r for r in inc["aois"] if r["included"]]
    countries = {c for r in included for c in r["countries"]}
    years = sorted(s["years"])
    return {
        "NAoiLabelled": str(s["n_aoi_labelled"]),
        "NAoi": str(s["n_aoi_included"]),
        "NActivations": str(s["n_activations_included"]),
        "NGroups": str(s["n_location_groups"]),
        "NCountries": str(len(countries)),
        "YearMin": years[0],
        "YearMax": years[-1],
        "GroupCap": f"{s['group_cap']:,}".replace(",", "{,}"),
        "KishGroups": f"{_kish(inc):.0f}",
        "NCells": f"{s['n_cells_after_cap']:,}".replace(",", "{,}"),
        "NNoEventProducts": str(s["n_no_event_products"]),
        "NUnspecifiedRiverine": str(s["n_unspecified_as_riverine"]),
    }


def lexicon_macros(path: Path) -> dict[str, str]:
    sys.path.insert(0, str(REPO / "backend" / "packages" / "topsem"))
    from topsem.world.lexicon import Lexicon

    lx = Lexicon(path)
    n = sum(len(v) for v in lx.by_lang.values())
    return {"NLexEntries": f"{n:,}".replace(",", "{,}"), "NLexLangs": str(len(lx.by_lang)),
            "LexSha": lx.sha256[:12]}


def freeze_macros(prereg: Path) -> dict[str, str]:
    head = prereg.read_text(encoding="utf-8").splitlines()[2]
    m = re.search(r"FROZEN\*\* (\S+) / commit ([0-9a-f]{7,40})", head)
    if not m:
        return {"FreezeCommit": missing("freeze commit"), "FreezeDateUTC": missing("freeze date")}
    return {"FreezeCommit": m.group(2)[:7], "FreezeDateUTC": m.group(1)}


def power_macros(null: dict | None, sesoi: dict | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if null:
        out |= {"PowerNReps": str(null["n_reps"]),
                "PowerEquivRate": pct(null["rate_equivalent_TOST"]),
                "PowerFalsePos": pct(null["rate_positive"]),
                "PowerHalfWidth": sci(null["median_ci90_halfwidth"], 2),
                "PowerNullNegative": str(null["decision_counts"].get("negative", 0))}
    else:
        out |= {k: missing(f"power null {k}") for k in
                ("PowerNReps", "PowerEquivRate", "PowerFalsePos", "PowerHalfWidth")}
    out["PowerPosRate"] = pct(sesoi["rate_positive"]) if sesoi else missing("power sesoi")
    return out


def nread_macros(meta: dict | None) -> dict[str, str]:
    if not meta:
        return {"NReadRevision": missing("N-read revision"),
                "NReadPromptSha": missing("N-read prompt sha")}
    return {"NReadRevision": meta["revision"][:7], "NReadPromptSha": meta["prompt_sha256"][:12],
            "NReadNames": f"{meta['n_names']:,}".replace(",", "{,}"),
            "NReadParsed": f"{100 * meta['parsed_ok_rate']:.2f}\\%"}


def geonames_macros() -> dict[str, str]:
    g = DATA / "global" / "geonames"
    sha = (g / "allCountries.zip.sha256")
    utc = (g / "retrieved_utc.txt")
    if not sha.is_file():
        return {"GeoNamesSha": missing("GeoNames sha256"), "GeoNamesUTC": missing("GeoNames UTC")}
    return {"GeoNamesSha": sha.read_text(encoding="utf-8").split()[0][:12],
            "GeoNamesUTC": utc.read_text(encoding="utf-8").strip()[:10] if utc.is_file()
            else missing("GeoNames UTC")}


def design_macros() -> dict[str, str]:
    """凍結前の設計判断の根拠となった数値（v32 のセル上限比較、pilot2 の柔軟比較の揺らぎ）。"""
    out: dict[str, str] = {}
    v32 = _json(REPO / "docs" / "verification" / "v32_cap_choice.json")
    if v32:
        draft = max(v32["table"], key=lambda r: r["cap"])  # 草案規則の上限（最大の行）
        out["KishDraftRule"] = f"{draft['kish_groups']:.0f}"
        out["CellsDraftRule"] = f"{draft['n_cells'] / 1e6:.1f}"
    else:
        out["KishDraftRule"] = missing("v32 draft Kish")
        out["CellsDraftRule"] = missing("v32 draft cells")
    noise, terr = [], []
    for f in sorted((DATA / "runs").glob("pilot2_d*/sim/rep*.npz")):
        with np.load(f) as z:
            y = z["y"].astype(float)
            bs = {a: float(np.mean((z[a].astype(float) - y) ** 2)) for a in ("M0", "M1", "M3")}
        noise.append(abs(bs["M1"] - bs["M3"]))
        terr.append(bs["M0"] - bs["M1"])
    out["PilotFlexNoise"] = (f"{100 * max(noise) / float(np.median(terr)):.0f}\\%" if noise
                             else missing("pilot2 flexible noise"))
    out["PilotFlexReps"] = str(len(noise)) if noise else missing("pilot2 reps")
    return out


def channel_macros() -> dict[str, str]:
    """名前チャネルの一致度（v34、ラベル不使用）。"""
    v = _json(REPO / "docs" / "verification" / "v34_channel_agreement.json")
    if not v:
        return {k: missing(f"v34 {k}") for k in ("LexWetShare", "NReadWetShare",
                                                  "NReadAgreeMin", "NReadAgreeMax")}
    same = [x["nread_same_class"] for x in v["lexicon_to_nread"].values()
            if x["nread_same_class"] is not None]
    return {"LexWetShare": pct(v["wetside"]["lexicon_share"]),
            "NReadWetShare": pct(v["wetside"]["nread_share"]),
            "NReadAgreeMin": pct(min(same)), "NReadAgreeMax": pct(max(same)),
            "LexWetGivenNRead": pct(v["wetside"]["lexicon_given_nread"]),
            "LexWetAfrica": pct(v["by_continent"]["AF"]["lexicon_wetside"]),
            "NReadWetAfrica": pct(v["by_continent"]["AF"]["nread_wetside"])}


def write_macros(m: dict[str, str]) -> None:
    lines = ["% 自動生成（figures/make_figures.py）。手で編集しない。",
             "\\providecommand{\\missing}[1]{\\textbf{[??\\,#1]}}"]
    lines += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(m.items())]
    (OUTDIR / "summary_numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


# %% 図 1（概念図）: 検定の設計 — 入れ子の比較・placebo・梯子・3 つの読み方・分割。データ不要
def fig_design(n_langs: str) -> None:
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig = plt.figure(figsize=(7.2, 3.3))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 46)
    ax.axis("off")
    blue, grey, light = "#1f6fb4", "0.45", "#e8f1f9"

    def box(x: float, y: float, w: float, h: float, text: str, *, fc: str = "white",
            ec: str = "0.3", ls: str = "-", size: float = 6.6, weight: str = "normal") -> None:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.8",
                                    fc=fc, ec=ec, lw=0.8, ls=ls))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size,
                weight=weight, linespacing=1.25)

    def arrow(x0: float, y0: float, x1: float, y1: float, col: str = "0.3") -> None:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=8,
                                     lw=0.9, color=col))

    ax.text(1, 44, "a", fontsize=9, weight="bold")
    ax.text(4, 44, "Nested comparison within places (primary, H2)", fontsize=7.2, weight="bold")
    y, h = 30, 9
    xs = [(1, 15, "Position\n(M0)"),
          (21, 22, "+ terrain, hydrology,\nland cover and\nsettlement layout (M1)"),
          (50, 21, "+ regional composition\nof names: which names\noccur in the area (M1c)"),
          (78, 21, "+ the name of each\nplace's nearest\nsettlements (M3n)")]
    for x, w, txt in xs:
        box(x, y, w, h, txt, fc=light if x == 78 else "white", ec=blue if x == 78 else "0.3")
    for (x0, w0, _), (x1, _, _) in pairwise(xs):
        arrow(x0 + w0 + 0.8, y + h / 2, x1 - 0.8, y + h / 2)
    ax.text(18, 27.5, "terrain increment $\\Delta_T$\n(the yardstick)", fontsize=6.2,
            ha="center", va="top", color=grey)
    ax.text(74.5, 27.5, "names increment $\\Delta_N$\n$R=\\Delta_N/\\Delta_T$", fontsize=6.2,
            ha="center", va="top", color=blue)
    box(78, 16.5, 21, 7, "placebo (M3np): the same\nnames permuted among the\narea's settlements; real "
        "names\nmust do better", ec=grey, ls="--", size=6.0)
    ax.text(1, 13.5, "b", fontsize=9, weight="bold")
    ax.text(4, 13.5, "Terrain-information ladder (H3)", fontsize=7.2, weight="bold")
    rungs = ["T0\nno\nelevation", "T1\n1 km", "T2\n250 m", "T3\n30 m\nstatistics",
             "T4 = M1\n+ HAND"]
    for k, r in enumerate(rungs):
        box(1 + k * 6.8, 2, 5.8, 7.5, r, size=5.4)
    ax.text(18, 11.2, "names added at every rung, within places", fontsize=6.0, ha="center",
            color=grey)
    ax.text(37, 13.5, "c", fontsize=9, weight="bold")
    ax.text(40, 13.5, "Three readings of names", fontsize=7.2, weight="bold")
    reads = [(("lexicon (N-lex)\nwater, wetland, crossing,\nlowland and upland words\nin "
               f"{n_langs} languages; knows\nno places"), 37.5),
             ("learned score (N-learn)\ncharacter model fitted\non training regions only", 54),
             (("language model (N-read)\nreads name elements;\ntested on floods after\nits "
               "training cut-off"), 70.5)]
    for txt, x in reads:
        box(x, 2, 15, 9.5, txt, size=5.4)
    ax.text(88, 13.5, "d", fontsize=9, weight="bold")
    ax.text(91, 13.5, "Held out", fontsize=7.2, weight="bold")
    box(88, 2, 11, 9.5, "places (primary)\ncountries\nlater floods\nfloods after the\n"
        "LM cut-off", size=5.4)
    fig.savefig(OUTDIR / "fig_design.pdf")
    plt.close(fig)


# %% 図 1: 採用 AOI の地図と年分布（ラベル単独の発生表だけで描ける）
def fig_map(inc: dict) -> None:
    import shapely
    from pyogrio import raw

    coast = DATA / "global" / "coast" / "ne_10m_coastline.zip"
    lines = []
    if coast.is_file():
        _, _, g, _ = raw.read(f"/vsizip/{coast}/ne_10m_coastline.shp")
        for line in shapely.from_wkb(g):
            lines += [np.asarray(part.coords) for part in getattr(line, "geoms", [line])]
    inc_rows = [r for r in inc["aois"] if r["included"]]
    exc_rows = [r for r in inc["aois"] if not r["included"]]

    def centre(r: dict) -> tuple[float, float]:
        w, s, e, n = r["bounds"]
        return (w + e) / 2, (s + n) / 2

    def draw(ax: plt.Axes, xlim: tuple[float, float], ylim: tuple[float, float],
             scale: float) -> None:
        for xy in lines:  # 海岸線はラスタ化（PDF を軽く保つ）
            ax.plot(xy[:, 0], xy[:, 1], color="0.75", lw=0.3, rasterized=True)
        if exc_rows:
            xy = np.array([centre(r) for r in exc_rows])
            ax.scatter(xy[:, 0], xy[:, 1], s=4 * scale, c="0.6", marker="x", lw=0.5,
                       label="excluded AOI")
        xy = np.array([centre(r) for r in inc_rows])
        size = np.sqrt([r["n_eval"] for r in inc_rows]) / 8 * scale
        ax.scatter(xy[:, 0], xy[:, 1], s=size, c="#1f6fb4", alpha=0.6, lw=0,
                   label="analysed AOI")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    fig = plt.figure(figsize=(7.2, 3.7))
    ax = fig.add_axes([0.0, 0.25, 0.72, 0.74])
    draw(ax, (-180, 180), (-58, 75), 1.0)
    ax.axis("off")
    ax.legend(loc="lower left", fontsize=7, frameon=False)
    ax.add_patch(plt.Rectangle((-12, 35), 45, 36, fill=False, lw=0.6, color="0.3"))
    axe = fig.add_axes([0.72, 0.31, 0.27, 0.62])
    draw(axe, (-12, 33), (35, 71), 0.5)
    axe.set_title("Europe", fontsize=7)
    for s in axe.spines.values():
        s.set_color("0.3")
        s.set_linewidth(0.6)
    ax2 = fig.add_axes([0.08, 0.07, 0.86, 0.16])
    fig.text(0.005, 0.95, "a", fontsize=9, fontweight="bold")
    fig.text(0.005, 0.22, "b", fontsize=9, fontweight="bold")
    yrs = Counter(r["event_time"][:4] for r in inc_rows)
    ks = sorted(yrs)
    ax2.bar(ks, [yrs[k] for k in ks], color="#1f6fb4")
    ax2.set_ylabel("AOIs", fontsize=7)
    ax2.tick_params(labelsize=6)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    fig.savefig(OUTDIR / "fig_map.pdf", dpi=300)
    plt.close(fig)


# %% 図 2: 主結果（凍結後）— 名前チャネルごとの比 R と区間、SESOI 帯
SESOI = 0.05
def _row(blk: dict | None, r_key: str, ci_key: str) -> tuple[float, list] | None:
    if not blk or blk.get(r_key) is None:
        return None
    return blk[r_key], blk.get(ci_key)


def _forest_rows(primary: dict, extra: dict[str, dict]) -> list[tuple[str, float, list, str]]:
    """(表示名, R̂, 95% CI, 様式) の列。様式 = main（確証）/ secondary（副次）/ undefined（物差しが
    不確か: 地形の増分の 95% CI が 0 を含み、プロトコル上 R は定義されない）。"""
    ch = primary.get("channels", {})
    spec = [
        ("Names within places (primary, H2)", _row(primary.get("primary"), "R", "ci95_R"), "main"),
        ("   lexicon only (N-lex)", _row(ch.get("M3nx"), "R", "ci95_R"), "main"),
        ("   learned only (N-learn)", _row(ch.get("M3nl"), "R", "ci95_R"), "main"),
        ("   placebo: names permuted", _row(primary.get("placebo"), "R_placebo",
                                             "ci95_R_placebo"), "main"),
        ("Position only, within places (H1)", _row(primary.get("H1"), "ratio_to_terrain",
                                                   "ci95_ratio"), "main"),
    ]
    if (met := extra.get("country")) and met.get("primary"):
        spec.append(("Countries held out", _row(met["primary"], "R", "ci95_R"),
                     "main"))
    if (met := extra.get("temporal")) and met.get("primary"):
        undefined = met["primary"]["decision"] == "undefined-denominator"
        spec.append(("Time held out (floods from 2024)" + (" *" if undefined else ""),
                     _row(met["primary"], "R", "ci95_R"), "undefined" if undefined else "main"))
    if (met := extra.get("llmcut")) and met.get("H4"):
        spec.append(("Language model, floods after its cut-off (H4)", _row(met["H4"], "R", "ci95_R"),
                     "main"))
        spec.append(("   lexicon and learned score, same floods", _row(met.get("primary"), "R",
                                                                       "ci95_R"), "main"))
    spec += [
        ("Total incl. regional composition (secondary)", _row(primary.get("total"), "R",
                                                               "ci95_R"), "secondary"),
        ("Flexible model refitted with names (secondary)", _row(primary.get("flexible"), "R",
                                                                 "ci95_R"), "secondary"),
    ]
    return [(n, v[0], v[1], s) for n, v, s in spec if v is not None]


def fig_main(rows: list[tuple[str, float, list, str]], decision: str) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 0.36 * len(rows) + 1.0))
    ax.axvspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axvline(0, color="0.4", lw=0.8)
    ys = np.arange(len(rows))[::-1]
    style = {"main": ("#1f6fb4", "-", True), "secondary": ("0.45", "-", False),
             "undefined": ("0.45", "--", False)}
    for y, (_, r, c, s) in zip(ys, rows, strict=True):
        col, ls, filled = style[s]
        if c is not None:
            ax.plot(c, [y, y], color=col, lw=1.6, ls=ls)
        ax.plot([r], [y], "o", ms=4.5, color=col, mfc=col if filled else "white")
    for y, (n, *_) in zip(ys, rows, strict=True):  # 群の区切り: 他の分割の前、副次比較の前
        if n.startswith(("Countries held out", "Total incl.")):
            ax.axhline(y + 0.5, color="0.8", lw=0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels([n for n, *_ in rows], fontsize=7)
    ax.set_xlabel("R = name increment / terrain increment (Brier score)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.text(SESOI, ys[-1] - 0.85, "SESOI", fontsize=6, ha="center", va="top", color="0.4")
    ax.text(-SESOI, ys[-1] - 0.85, "−SESOI", fontsize=6, ha="center", va="top", color="0.4")
    ax.set_ylim(ys[-1] - 1.2, ys[0] + 0.6)
    ax.set_title(f"Pre-specified verdict (H2): {decision}", fontsize=8, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_main.pdf")
    plt.close(fig)


def fig_ladder(ladder: dict, terrain: float) -> None:
    """図 3: 地形情報の梯子（T0 → T4）での名前の場所内の増分（R 単位）と 95% CI、SESOI 帯。"""
    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    x = np.arange(len(ladder["rungs"]))
    d = np.array(ladder["delta_names"]) / terrain
    c = np.array(ladder["ci95_delta_names"]) / terrain
    ax.axhspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.errorbar(x, d, yerr=[d - c[:, 0], c[:, 1] - d], fmt="o-", color="#1f6fb4", ms=4,
                lw=1.2, capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels(["T0\nno elevation", "T1\n1 km", "T2\n250 m", "T3\n30 m stats",
                        "T4\n+HAND"], fontsize=7)
    ax.set_ylabel("names within places\n(fraction of terrain increment)", fontsize=8)
    ax.set_ylim(-0.08, 0.08)
    ax.tick_params(axis="y", labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_ladder.pdf")
    plt.close(fig)


def rv(x: float) -> str:
    """比 R の表記（小数 3 桁、数式モード）。0 でない値が 0.000 に丸まる場合だけ 4 桁にする
    （区間の端が 0 を挟むか挟まないかが読めるように）。"""
    return f"\\ensuremath{{{x:.4f}}}" if 0 < abs(x) < 0.0005 else f"\\ensuremath{{{x:.3f}}}"


def _ci_macros(prefix: str, val: float, ci95: list | None) -> dict[str, str]:
    out = {prefix: rv(val)}
    if ci95:
        out |= {f"{prefix}Lo": rv(ci95[0]), f"{prefix}Hi": rv(ci95[1])}
    return out


def detail_macros(m: dict) -> dict[str, str]:
    """主分割の詳細（記述と副次、§6）。全て metrics_primary.json から。"""
    p = m["primary"]
    d01 = p["delta_bs_M0_M1"]
    lo01 = p["ci95_delta_M0_M1"][0]  # 保守的なマージン（非劣性試験の固定マージン法の慣行、事後の確認）
    out = {"DeltaTerrainLo": sci(lo01, 2), "RConsLo": rv(p["ci90_delta"][0] / lo01),
           "RConsHi": rv(p["ci90_delta"][1] / lo01)}
    out |= {"DeltaTerrain": sci(d01, 2), "DeltaNames": sci(p["delta_bs_M1c_M3n"], 2),
           "DeltaNamesLo": sci(p["ci95_delta"][0], 2), "DeltaNamesHi": sci(p["ci95_delta"][1], 2)}
    pl = m.get("placebo", {})
    out |= _ci_macros("RPlacebo", pl["R_placebo"], pl.get("ci95_R_placebo"))
    out |= _ci_macros("RMinusRp", pl["R_minus_Rp"], pl.get("ci95_R_minus_Rp"))
    ch = m.get("channels", {})
    out |= _ci_macros("RLex", ch["M3nx"]["R"], ch["M3nx"]["ci95_R"])
    out |= _ci_macros("RLearn", ch["M3nl"]["R"], ch["M3nl"]["ci95_R"])
    out |= _ci_macros("RHOne", m["H1"]["ratio_to_terrain"], m["H1"]["ci95_ratio"])
    out |= _ci_macros("RHOneFlex", m["H1_flexible"]["ratio_to_terrain"],
                      m["H1_flexible"]["ci95_ratio"])
    out |= _ci_macros("RTotal", m["total"]["R"], m["total"]["ci95_R"])
    out["VerdictTotal"] = m["total"]["decision"]
    out |= _ci_macros("RFlex", m["flexible"]["R"], m["flexible"]["ci95_R"])
    out["VerdictFlex"] = m["flexible"]["decision"]
    b = m["brier"]
    out["RComposition"] = rv((b["M1s"] - b["M1c"]) / d01)
    out["RCompositionPos"] = rv((b["M0"] - b["M0c"]) / d01)
    lad = m["ladder"]
    out |= _ci_macros("RTZero", lad["delta_names"][0] / d01,
                      [x / d01 for x in lad["ci95_delta_names"][0]])
    out["LadderP"] = f"{lad['p_one_sided_delta0_gt_delta4']:.2f}"
    # §6: T0 での地形情報の回収率 ΔBS(T0c → T0n) / ΔBS(T0 → T4)
    out |= _ci_macros("RecoverTZero", lad["terrain_info_recovered_at_T0"], lad["ci95_recovered_at_T0"])
    # §6: 比そのものの bootstrap 区間（物差しの揺らぎを含む、記述）
    rb = p["ci95_R_ratio_bootstrap"]
    out |= {"RRatioLo": rv(rb[0]), "RRatioHi": rv(rb[1])}
    out["LadderRho"] = rv(lad["spearman_k_delta"])
    tw = m.get("two_way", {}).get("ci95_R")
    if tw:
        out |= {"RTwoWayLo": rv(tw[0]), "RTwoWayHi": rv(tw[1])}
    sm = m["secondary_metrics"]
    out |= {"CaptureMOnec": f"{100 * sm['M1c']['capture10']:.1f}\\%",
            "CaptureMThreen": f"{100 * sm['M3n']['capture10']:.1f}\\%",
            "PrAucMOnec": f"{sm['M1c']['pr_auc']:.3f}", "PrAucMThreen": f"{sm['M3n']['pr_auc']:.3f}"}
    h5 = m.get("H5_strata", {})
    for key, name in (("nearest_name_in_lexicon", "LexMatched"), ("nearest_name_wetside", "Wetside"),
                      ("low_relief", "LowRelief")):
        if key in h5:
            out[f"R{name}"] = rv(h5[key]["R_within"])
            out[f"NCells{name}"] = f"{h5[key]['n']:,}".replace(",", "{,}")
    bc = [v["R"] for v in m.get("by_continent", {}).values() if v.get("R") is not None]
    if bc:
        out |= {"RContMin": rv(min(bc)), "RContMax": rv(max(bc))}
    return out


def split_macros(extra: dict[str, dict]) -> dict[str, str]:
    """副分割（国外し・時間外し・llmcut）の判定と記述。全て metrics_<split>.json から。"""
    out: dict[str, str] = {}
    for split, name in (("country", "Country"), ("temporal", "Temporal"), ("llmcut", "HFour")):
        met = extra.get(split)
        blk = (met or {}).get("H4" if split == "llmcut" else "primary")
        if not blk:
            out |= {f"R{name}": missing(f"{split} R"), f"Verdict{name}": missing(f"{split} verdict")}
            continue
        out |= _ci_macros(f"R{name}", blk["R"], blk["ci95_R"])
        out[f"Verdict{name}"] = blk["decision"]
        out[f"NCells{name}"] = f"{met['n_cells']:,}".replace(",", "{,}")
        out[f"NUnits{name}"] = str(met["n_units"])
        # 物差し（地形の増分）と名前の増分の絶対値（R が定義されない場合の記述用）
        out |= {f"DeltaTerrain{name}": sci(blk["delta_bs_M0_M1"], 2),
                f"DeltaTerrain{name}Lo": sci(blk["ci95_delta_M0_M1"][0], 2),
                f"DeltaTerrain{name}Hi": sci(blk["ci95_delta_M0_M1"][1], 2)}
        dkey = next(k for k in blk if k.startswith("delta_bs_") and k != "delta_bs_M0_M1")
        out |= {f"DeltaNames{name}": sci(blk[dkey], 2),
                f"DeltaNames{name}Lo": sci(blk["ci95_delta"][0], 2),
                f"DeltaNames{name}Hi": sci(blk["ci95_delta"][1], 2)}
    # 名前の地域構成の効き方の分割間の比較（記述）: 位置だけのモデルへ（M0 → M0c）と地形モデルへ（M1s → M1c）
    comp_pos, comp_ter = {}, {}
    for split, met in extra.items():
        b = met["brier"]
        d = b["M0"] - b["M1"]
        comp_pos[split] = (b["M0"] - b["M0c"]) / d
        comp_ter[split] = (b["M1s"] - b["M1c"]) / d
        out[f"DeltaComp{split.capitalize()}"] = sci(b["M1s"] - b["M1c"], 2)
    if comp_pos:
        out |= {"RCompPosSplitMin": rv(min(comp_pos.values())),
                "RCompPosSplitMax": rv(max(comp_pos.values()))}
        out |= {f"RCompTerrain{s.capitalize()}": rv(v) for s, v in comp_ter.items()}
    met = extra.get("llmcut")
    if met and met.get("H4"):
        d01 = met["H4"]["delta_bs_M0_M1"]
        out |= _ci_macros("RLlmLex", met["primary"]["R"], met["primary"]["ci95_R"])
        out["VerdictLlmLex"] = met["primary"]["decision"]
        tot, tpl = met["H4_total"], met["H4_total_placebo"]
        out |= _ci_macros("RHFourTotal", tot["R"], tot["ci95_R"])
        out["VerdictHFourTotal"] = tot["decision"]
        out |= _ci_macros("RHFourTotalMinusRp", tpl["R_minus_Rp"], tpl["ci95_R_minus_Rp"])
        b = met["brier"]
        out["RHFourComp"] = rv((b["M1s"] - b["M1rc"]) / d01)
        out["RLlmLexComp"] = rv((b["M1s"] - b["M1c"]) / d01)
        fl, fpl = met["flexible"], met["flexible_placebo"]
        out |= _ci_macros("RFlexLlm", fl["R"], fl["ci95_R"])
        out["VerdictFlexLlm"] = fl["decision"]
        out |= _ci_macros("RFlexLlmPlacebo", fpl["R_placebo"], fpl["ci95_R_placebo"])
    return out


# %% 図 4（凍結後・探索的 v35 + 記述 H1）: 名前が運ぶもの — 場所ではなく地域
def contrast_macros(v: dict | None) -> dict[str, str]:
    """v35（探索的）: AOI 内の対比（百分率ポイント）と AOI 間の順位相関。"""
    if not v:
        return {k: missing(f"v35 {k}") for k in ("ContrastWetRate", "ContrastOtherRate",
                                                  "ContrastDiff", "BetweenRho")}
    w, h, b = v["wetside"], v["low_hand"], v["between_aoi"]

    def pp(x: float) -> str:
        return f"\\ensuremath{{{100 * x:+.2f}}}".replace("+", "{+}")

    out = {"ContrastWetRate": f"{100 * w['mean_rate_flagged']:.1f}\\%",
           "ContrastOtherRate": f"{100 * w['mean_rate_other']:.1f}\\%",
           "NContrastAoi": str(w["n_aoi"]), "NContrastPlaces": str(w["n_units"]),
           "NBetweenAoi": str(b["n_aoi"])}
    for key, val in (("ContrastDiff", w["raw"]), ("ContrastAdj", w["M1c"]),
                     ("ContrastHand", h["raw"])):
        out |= {key: pp(val["estimate"]), f"{key}Lo": pp(val["ci95"][0]),
                f"{key}Hi": pp(val["ci95"][1])}
    for key, val in (("BetweenRho", b["rho_obs"]), ("BetweenRhoResid", b["rho_resid_M1"])):
        out |= {key: rv(val["estimate"]), f"{key}Lo": rv(val["ci95"][0]),
                f"{key}Hi": rv(val["ci95"][1])}
    t = b["tertiles"]
    out |= {"TertileWetLow": f"{100 * t[0]['obs']:.1f}\\%",
            "TertileWetHigh": f"{100 * t[2]['obs']:.1f}\\%"}
    return out


def fig_names(v: dict, m: dict) -> None:
    """a: AOI 内の対比（地形の物差しと水辺名）、b: AOI 間の散布、c: 地域構成と場所内の増分。"""
    fig = plt.figure(figsize=(7.2, 2.9), layout="constrained")
    axs = fig.subplots(1, 3, width_ratios=[1.1, 0.85, 1.05])
    ax = axs[0]
    rows = [("terrain: low vs high HAND", v["low_hand"]["raw"], "0.35"),
            ("names: water-side vs other", v["wetside"]["raw"], "#1f6fb4"),
            ("names, net of position\nand regional composition", v["wetside"]["M0c"], "#1f6fb4"),
            ("names, net of terrain model", v["wetside"]["M1c"], "#1f6fb4")]
    ys = np.arange(len(rows))[::-1]
    for y, (_, r, col) in zip(ys, rows, strict=True):
        c = 100 * np.array(r["ci95"])
        ax.plot(c, [y, y], color=col, lw=1.6)
        ax.plot([100 * r["estimate"]], [y], "o", ms=4, color=col)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([n for n, *_ in rows], fontsize=6.5)
    ax.set_xlabel("within-AOI difference in\nflooded cells (pct. points)", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.5)
    ax = axs[1]
    b = v["between_aoi"]
    ax.scatter(100 * np.array(b["share"]), 100 * np.array(b["obs"]), s=5, c="#1f6fb4",
               alpha=0.45, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel("water-side names among\nnearest settlements (%)", fontsize=7)
    ax.set_ylabel("flooded cells in the AOI (%)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    r = b["rho_obs"]
    txt = (f"Spearman ρ = {r['estimate']:.3f}\n[{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]"
           ).replace("-", "\u2212")
    ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=6.5, ha="right", va="top")
    ax = axs[2]
    d01 = m["primary"]["delta_bs_M0_M1"]
    br = m["brier"]
    bars = [("regional composition\ngiven position", (br["M0"] - br["M0c"]) / d01, None),
            ("names, flexible model,\ngiven position", m["H1_flexible"]["ratio_to_terrain"],
             m["H1_flexible"]["ci95_ratio"]),
            ("names within places\ngiven position", m["H1"]["ratio_to_terrain"],
             m["H1"]["ci95_ratio"]),
            ("regional composition\ngiven terrain", (br["M1s"] - br["M1c"]) / d01, None),
            ("names within places\ngiven terrain", m["primary"]["R"], m["primary"]["ci95_R"])]
    ys = np.arange(len(bars))[::-1]
    for y, (_, val, c) in zip(ys, bars, strict=True):
        ax.barh(y, val, color="#9cc3e4", height=0.6)
        if c:
            ax.plot(c, [y, y], color="0.2", lw=1)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([n for n, *_ in bars], fontsize=6.5)
    ax.set_xlabel("Brier-score improvement\n(× terrain increment)", fontsize=7)
    ax.tick_params(axis="x", labelsize=6.5)
    for a, lab in zip(axs, "abc", strict=True):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.set_title(lab, fontsize=9, fontweight="bold", loc="left")
    fig.savefig(OUTDIR / "fig_names.pdf")
    plt.close(fig)


# %% SI の表（凍結後）— 全比較・分割・感度・層別・大陸別。数値は全て metrics_*.json から
SENS = [("s1_alltypes", "All flood types, not only riverine"),
        ("s2a_nocap", "No cap on cells per place"),
        ("s2b_aoicap", "At most 2,000 cells per AOI instead of per place"),
        ("s3_plts", "Names of feature classes P, L, T and S"),
        ("s4_footprint", "Products without an image footprint excluded"),
        ("s5_k1", "One seed instead of five"),
        ("s6_buffer_b1", "Places merged within 10\\,km instead of 5\\,km$^{a}$"),
        ("s7_post2015", "Events from 2016 onward only")]
TO = "\\to"
CONTINENTS = {"AF": "Africa", "AS": "Asia", "EU": "Europe", "NA": "North America",
              "OC": "Oceania", "SA": "South America"}


def _r(x: float | None) -> str:
    if x is None:
        return "---"
    s = f"{x:+.4f}" if 0 < abs(x) < 0.0005 else f"{x:+.3f}"  # 非零が 0.000 に丸まらないように
    return f"${s}$".replace("+", "{+}").replace("-", "{-}")


def _ci(c: list | None) -> str:
    return "---" if not c else f"[{_r(c[0])}, {_r(c[1])}]"


def _n(x: int) -> str:
    return f"{x:,}".replace(",", "{,}")


def _table(cols: str, head: str, rows: list[str], caption: str, label: str) -> str:
    body = "\n".join(rows)
    return (f"\\begin{{table}}[h!]\n\\centering\\footnotesize\\setlength{{\\tabcolsep}}{{4pt}}\n"
            f"\\caption{{{caption}}}\\label{{{label}}}\n"
            f"\\begin{{adjustbox}}{{max width=\\linewidth}}\n"
            f"\\begin{{tabular}}{{{cols}}}\n\\toprule\n{head} \\\\\n\\midrule\n{body}\n"
            f"\\bottomrule\n\\end{{tabular}}\n\\end{{adjustbox}}\n\\end{{table}}\n")


def _tost(blk: dict) -> str:
    """§6: TOST の成否を判定 label と独立に併記する。"""
    return "---" if "equivalent_TOST" not in blk else ("yes" if blk["equivalent_TOST"] else "no")


def _cmp_row(name: str, blk: dict | None, r: str = "R", c95: str = "ci95_R",
             c90: str | None = "ci90_R") -> str:
    if not blk:
        return f"{name} & \\missing{{{name}}} & & & & \\\\"
    v90 = _ci(blk.get(c90)) if c90 else "---"
    return (f"{name} & {_r(blk.get(r))} & {_ci(blk.get(c95))} & {v90} & "
            f"{blk.get('decision', '---')} & {_tost(blk)} \\\\")


def si_comparisons(m: dict) -> str:
    """補足表 1: 主分割の全比較（副次・チャネル・梯子を含む）。"""
    d01 = m["primary"]["delta_bs_M0_M1"]
    ch, chf = m.get("channels", {}), m.get("channels_flexible", {})
    pl, tpl, fpl = m.get("placebo", {}), m.get("total_placebo", {}), m.get("flexible_placebo", {})
    rows = [
        _cmp_row("H2: names within places, M1c$\\to$M3n", m["primary"]),
        _cmp_row("\\quad placebo, M1c$\\to$M3np", pl, "R_placebo", "ci95_R_placebo", None),
        _cmp_row("\\quad real minus placebo", pl, "R_minus_Rp", "ci95_R_minus_Rp", None),
        _cmp_row("\\quad lexicon only, M1c$\\to$M3nx", ch.get("M3nx"), c90=None),
        _cmp_row("\\quad learned score only, M1c$\\to$M3nl", ch.get("M3nl"), c90=None),
        _cmp_row("H1: position only, M0c$\\to$M0n", m.get("H1"), "ratio_to_terrain", "ci95_ratio",
                 None),
        _cmp_row("\\quad flexible, M0$\\to$M0N", m.get("H1_flexible"), "ratio_to_terrain",
                 "ci95_ratio", None),
        _cmp_row("H2t: total increment, M1s$\\to$M3s", m.get("total")),
        _cmp_row("\\quad placebo, M1s$\\to$M3ps", tpl, "R_placebo", "ci95_R_placebo", None),
        _cmp_row("Flexible refit, M1$\\to$M3", m.get("flexible")),
        _cmp_row("\\quad placebo, M1$\\to$M3p", fpl, "R_placebo", "ci95_R_placebo", None),
        _cmp_row("\\quad lexicon only, M1$\\to$M3x", chf.get("M3x"), c90=None),
        _cmp_row("\\quad learned score only, M1$\\to$M3l", chf.get("M3l"), c90=None),
    ]
    for key, label in (("ladder", "H3 ladder, within places"), ("ladder_flexible",
                                                                "H3 ladder, flexible")):
        lad = m.get(key)
        if not lad:
            continue
        for k, (rung, dlt, c) in enumerate(zip(lad["rungs"], lad["delta_names"],
                                               lad["ci95_delta_names"], strict=True)):
            name = f"{label}, T{k}" if k < 4 else f"{label}, T4 = M1"
            arrow = rung.replace("->", "$" + TO + "$")
            rows.append(f"{name} ({arrow}) & {_r(dlt / d01)} & "
                        f"{_ci([x / d01 for x in c])} & --- & --- & --- \\\\")
    lo01 = m["primary"]["ci95_delta_M0_M1"][0]
    c90 = [x / lo01 for x in m["primary"]["ci90_delta"]]
    rows.append(f"H2, conservative margin (lower 95\\% bound of the terrain increment; post hoc) & "
                f"{_r(m['primary']['delta_bs_M1c_M3n'] / lo01)} & --- & {_ci(c90)} & --- & "
                f"{'yes' if all(abs(x) < SESOI for x in c90) else 'no'} \\\\")
    tw = m.get("two_way", {})
    if tw.get("ci95_R"):
        rows.append(f"H2, two-way bootstrap ({tw['n_activations']} activations $\\times$ "
                    f"{tw['n_countries']} countries) & {_r(m['primary']['R'])} & "
                    f"{_ci(tw['ci95_R'])} & --- & --- & --- \\\\")
    return _table("p{5.9cm}rcclc", "Comparison & $R$ & 95\\% interval & 90\\% interval & Verdict & TOST",
                  rows, "All pre-specified comparisons on the primary split (location groups held "
                  "out). $R$ is the change in held-out Brier score as a fraction of the terrain "
                  "increment $\\Delta\\mathrm{BS}(\\mathrm{M0}\\to\\mathrm{M1})$; intervals are "
                  "cluster-bootstrap percentile intervals over location groups (fixed-margin "
                  "method; the two-way row resamples activations and countries instead). Verdicts are "
                  "given for the comparisons to which the protocol attaches a decision rule.",
                  "tab:si-comparisons")


def si_splits(primary: dict, extra: dict[str, dict]) -> str:
    """補足表 2: 分割ごとの主比較（llmcut は H4 と主比較の両方）。"""
    rows = []
    spec = [("Location groups held out (primary)", primary, "primary"),
            ("Countries held out", extra.get("country"), "primary"),
            ("Time held out", extra.get("temporal"), "primary"),
            ("Language-model cut-off: H4, M1rc$\\to$M3rn", extra.get("llmcut"), "H4"),
            ("Language-model cut-off: lexicon and learned score", extra.get("llmcut"), "primary")]
    for name, met, key in spec:
        blk = (met or {}).get(key)
        if not blk:
            rows.append(f"{name} & \\missing{{pending}} & & & & & & \\\\")
            continue
        rows.append(f"{name} & {_n(met['n_cells'])} & {met['n_units']} & "
                    f"{sci(blk['delta_bs_M0_M1'], 2)} & {_r(blk['R'])} & {_ci(blk['ci95_R'])} & "
                    f"{blk['decision']} & {_tost(blk)} \\\\")
    return _table("p{4.0cm}rrcrclc", "Split & Cells & Units & Terrain $\\Delta$BS & $R$ & 95\\% "
                  "interval & Verdict & TOST", rows,
                  "The primary comparison under each pre-specified split. Units are the clusters "
                  "of the bootstrap (location groups; countries for the country hold-out). The "
                  "time hold-out evaluates floods from the latest years on models trained on earlier "
                  "floods elsewhere; the language-model split evaluates floods after the model's "
                  "pre-training cut-off.", "tab:si-splits")


def si_splits_secondary(primary: dict, extra: dict[str, dict]) -> str:
    """補足表 3: 分割ごとの副次比較（全体の増分・柔軟モデル、llmcut は N-read 版も）と placebo との差。"""
    spec = [("Location groups", primary), ("Countries", extra.get("country")),
            ("Time", extra.get("temporal")), ("Language-model cut-off", extra.get("llmcut"))]
    comps = [("total", "total_placebo", "total, M1s$\\to$M3s"),
             ("flexible", "flexible_placebo", "flexible, M1$\\to$M3"),
             ("H4_total", "H4_total_placebo", "language model, total, M1s$\\to$M3rs"),
             ("H4_flexible", "H4_flexible_placebo", "language model, flexible, M1$\\to$M3r")]
    rows = []
    for split, met in spec:
        if not met:
            rows.append(f"{split} & \\missing{{pending}} & & & & \\\\")
            continue
        for key, pkey, lab in comps:
            blk, pl = met.get(key), met.get(pkey) or {}
            if not blk:
                continue
            rows.append(f"{split}: {lab} & {_r(blk['R'])} & {_ci(blk['ci95_R'])} & "
                        f"{_ci(pl.get('ci95_R_minus_Rp'))} & {blk['decision']} & {_tost(blk)} \\\\")
    return _table("p{5.6cm}rcclc", "Split: comparison & $R$ & 95\\% interval & real $-$ placebo "
                  "(95\\%) & Verdict & TOST", rows,
                  "Secondary comparisons under each split: the total increment of names over the "
                  "recalibrated terrain model (including the regional composition of names) and the "
                  "flexible model refitted with the name features, for the lexicon and learned score "
                  "and, on the language-model split, for the language-model reading. The fourth column "
                  "is the interval of $R$ minus the same comparison with names permuted among the "
                  "settlements of each AOI.", "tab:si-splits-secondary")


def si_sensitivity(runs: Path) -> str:
    """補足表 4: §7 の感度解析（主判定の腕だけ、主分割）。s6 は B-1 の修正版の直下に元の実行を注記。"""
    rows = []
    for tag, name in SENS:
        met = _json(runs / f"world_v2_{tag}" / "metrics_primary.json")
        blk = (met or {}).get("primary")
        if not blk:
            rows.append(f"{name} & \\missing{{sensitivity {tag[:2]}}} & & & & & & \\\\")
            continue
        rows.append(f"{name} & {_n(met['n_cells'])} & {met['n_units']} & "
                    f"{sci(blk['delta_bs_M0_M1'], 2)} & {_r(blk['R'])} & {_ci(blk['ci90_R'])} & "
                    f"{blk['decision']} & {_tost(blk)} \\\\")
        first = _json(runs / "world_v2_s6_buffer" / "metrics_primary.json") if tag == "s6_buffer_b1" else None
        if first and first.get("primary"):
            b = first["primary"]
            rows.append(f"\\quad as first run, folds by 5\\,km places$^{{a}}$ & {_n(first['n_cells'])} & "
                        f"{first['n_units']} & {sci(b['delta_bs_M0_M1'], 2)} & {_r(b['R'])} & "
                        f"{_ci(b['ci90_R'])} & {b['decision']} & {_tost(b)} \\\\")
    return _table("p{4.0cm}rrcrclc", "Sensitivity analysis & Cells & Places & Terrain $\\Delta$BS "
                  "& $R$ & 90\\% interval & Verdict & TOST", rows,
                  "Pre-specified sensitivity analyses (protocol section 7), each re-running the "
                  "arms of the primary verdict on the primary split and applying the same decision "
                  "table. $^{a}$The first run of this analysis drew cells by 10\\,km places but, "
                  "through an implementation defect, assigned folds and bootstrap units by the "
                  "5\\,km places of the primary analysis; it was corrected and re-run after the "
                  "first result was seen (amendment B-1), and both runs are shown.",
                  "tab:si-sensitivity")


def si_strata(m: dict) -> str:
    """補足表 4: H5 の層（記述）と大陸別（記述）。"""
    h5 = m.get("H5_strata", {})
    names = {"nearest_name_in_lexicon": "Nearest name in the lexicon",
             "not_nearest_name_in_lexicon": "Nearest name not in the lexicon",
             "nearest_name_wetside": "Nearest name water-side",
             "not_nearest_name_wetside": "Nearest name not water-side",
             "low_relief": "Slope at or below the median",
             "not_low_relief": "Slope above the median",
             "built_up_gt_10pct": "Built-up cover above 10\\%",
             "not_built_up_gt_10pct": "Built-up cover at most 10\\%"}
    rows = [f"{lab} & {_n(h5[k]['n'])} & {100 * h5[k]['wet_rate']:.1f} & {_r(h5[k]['R_within'])} \\\\"
            for k, lab in names.items() if k in h5]
    rows.append("\\midrule")
    rows += [f"{CONTINENTS.get(k, k)} & {_n(v['n'])} & --- & {_r(v['R'])} \\\\"
             for k, v in sorted(m.get("by_continent", {}).items())]
    return _table("lrrr", "Stratum & Cells & Flooded (\\%) & $R$ within stratum", rows,
                  "Descriptive strata (H5) and continents for the primary comparison. $R$ within a "
                  "stratum is the names increment in the stratum divided by the terrain increment in "
                  "the same stratum; continents are those of each cell's nearest settlement. No "
                  "intervals are given because these breakdowns are descriptive.", "tab:si-strata")


def si_metrics(m: dict) -> str:
    """補足表 5: 副次評価量。"""
    sm = m.get("secondary_metrics", {})
    names = {"M1": "M1 (terrain model)", "M3": "M3 (refitted with names)",
             "M1s": "M1s (recalibrated M1)", "M3s": "M3s (+ names)",
             "M1c": "M1c (+ regional composition)", "M3n": "M3n (+ names within places)"}
    rows = [f"{lab} & {sm[k]['brier']:.5f} & {sm[k]['log_loss']:.5f} & {sm[k]['pr_auc']:.3f} & "
            f"{100 * sm[k]['capture10']:.1f} \\\\" for k, lab in names.items() if k in sm]
    return _table("lrrrr", "Model & Brier & Log loss & PR-AUC & Top-10\\% capture (\\%)", rows,
                  "Secondary metrics on the primary split (pooled over held-out folds).",
                  "tab:si-metrics")


def sens_macros(runs: Path) -> dict[str, str]:
    """§7 の感度の判定の要約（本文用）。全て metrics_primary.json から。"""
    got = {tag: (m or {}).get("primary") for tag, _ in SENS
           if (m := _json(runs / f"world_v2_{tag}" / "metrics_primary.json"))}
    out = {"NSens": str(len(SENS)), "NSensDone": str(len(got))}
    if len(got) < len(SENS):
        MISSING.append(f"sensitivity runs {len(got)}/{len(SENS)}")
    out["NSensEquiv"] = str(sum(1 for b in got.values() if b and b["decision"] == "equivalent-null"))
    post = got.get("s7_post2015")
    if post:
        out |= _ci_macros("RSensPost", post["R"], post["ci95_R"])
        out["VerdictSensPost"] = post["decision"]
        out["TostSensPost"] = "held" if post["equivalent_TOST"] else "did not hold"
    rs = [b["R"] for b in got.values() if b and b["decision"] != "undefined-denominator"]
    if rs:  # R が定義される（物差しの 95% CI が 0 を含まない）感度だけの範囲
        out |= {"RSensMin": rv(min(rs)), "RSensMax": rv(max(rs))}
    out["NSensDefined"] = str(len(rs))
    nocap = got.get("s2a_nocap")
    if nocap:  # §7-2a: 上限なし。少数の巨大な場所が支配するので物差しが不確かになり得る（記述）
        out["VerdictNoCap"] = nocap["decision"]
        out |= {"DeltaTerrainNoCap": sci(nocap["delta_bs_M0_M1"], 2),
                "DeltaTerrainNoCapLo": sci(nocap["ci95_delta_M0_M1"][0], 2),
                "DeltaTerrainNoCapHi": sci(nocap["ci95_delta_M0_M1"][1], 2),
                "DeltaNamesNoCap": sci(nocap["delta_bs_M1c_M3n"], 2),
                "DeltaNamesNoCapLo": sci(nocap["ci95_delta"][0], 2),
                "DeltaNamesNoCapHi": sci(nocap["ci95_delta"][1], 2)}
        cells: Counter = Counter()
        for f in sorted((runs / "world_v2_s2a_nocap" / "predictions" / "primary").glob("fold*.npz")):
            u, c = np.unique(np.load(f)["unit"], return_counts=True)
            cells.update(dict(zip(u.tolist(), c.tolist(), strict=True)))
        if cells:
            n = np.array(list(cells.values()), dtype=float)
            out |= {"KishNoCap": f"{n.sum() ** 2 / (n ** 2).sum():.0f}",
                    "LargestShareNoCap": pct(n.max() / n.sum())}
    return out


def write_si_tables(primary: dict | None, extra: dict[str, dict]) -> None:
    """補足表を表ごとの別ファイルに（si.tex が対応する Note の直後で読む）。"""
    tables = {"comparisons": lambda: si_comparisons(primary),
              "splits": lambda: si_splits(primary, extra),
              "splits_secondary": lambda: si_splits_secondary(primary, extra),
              "sensitivity": lambda: si_sensitivity(DATA / "runs"),
              "strata": lambda: si_strata(primary), "metrics": lambda: si_metrics(primary)}
    incomplete = False
    for name, make in tables.items():
        text = make() if primary else f"\\missing{{table {name}}}"
        incomplete |= "\\missing{" in text
        (OUTDIR / f"si_table_{name}.tex").write_text(
            "% generated by figures/make_figures.py — do not edit\n" + text + "\n", encoding="utf-8")
    if incomplete:
        MISSING.append("SI tables (results not yet complete)")


def result_macros(primary: dict | None) -> dict[str, str]:
    if not primary or "primary" not in primary:
        return {k: missing(f"result {k}") for k in ("RHat", "RLo", "RHi", "RLoNinety",
                                                     "RHiNinety", "Verdict")}
    p = primary["primary"]

    def rp(x: float) -> str:  # 要旨用の百分率（R × 100、小数 1 桁）
        return f"\\ensuremath{{{100 * x:.1f}\\%}}"

    return {"RHat": rv(p["R"]), "RLo": rv(p["ci95_R"][0]), "RHi": rv(p["ci95_R"][1]),
            "RLoNinety": rv(p["ci90_R"][0]), "RHiNinety": rv(p["ci90_R"][1]),
            "RHatPct": rp(p["R"]), "RLoPct": rp(p["ci95_R"][0]), "RHiPct": rp(p["ci95_R"][1]),
            "Verdict": p["decision"]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--power-null", default="power_v3_null")
    ap.add_argument("--power-sesoi", default="power_v3_sesoi")
    ap.add_argument("--eval", default="world_v2")
    args = ap.parse_args()

    m: dict[str, str] = {}
    m |= inventory_macros(_json(REPO / "docs" / "data" / "ems_inventory_20261025.json"))
    inc = _json(DATA / "global" / "incidence.json")
    if inc:
        m |= incidence_macros(inc)
        m |= sample_macros(inc)
        fig_map(inc)
    else:
        m |= {k: missing(f"incidence {k}") for k in
              ("NAoiLabelled", "NAoi", "NActivations", "NGroups", "NCountries", "YearMin",
               "YearMax", "GroupCap", "KishGroups", "NCells", "NNoEventProducts",
               "NUnspecifiedRiverine")}
    m |= lexicon_macros(REPO / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml")
    fig_design(m["NLexLangs"])
    m |= design_macros()
    m |= channel_macros()
    m |= freeze_macros(REPO / "docs" / "protocol" / "preregistration_v2.md")
    m |= geonames_macros()
    m |= nread_macros(_json(DATA / "global" / "nread" / "nread_meta.json"))
    m |= power_macros(_json(DATA / "runs" / args.power_null / "power_summary.json"),
                      _json(DATA / "runs" / args.power_sesoi / "power_summary.json"))
    run = DATA / "runs" / args.eval
    primary = _json(run / "metrics_primary.json")
    m |= result_macros(primary)
    if primary:
        extra = {s: v for s in ("country", "temporal", "llmcut")
                 if (v := _json(run / f"metrics_{s}.json"))}
        m |= detail_macros(primary)
        m |= split_macros(extra)
        fig_main(_forest_rows(primary, extra), primary["primary"]["decision"])
        if primary.get("ladder"):
            fig_ladder(primary["ladder"], primary["primary"]["delta_bs_M0_M1"])
        write_si_tables(primary, extra)
        m |= sens_macros(DATA / "runs")
        diag = run / "diagnostics_primary.png"  # 補足図 1（集計 world_aggregate が描く診断図）
        if diag.is_file():
            shutil.copyfile(diag, OUTDIR / "si_fig_diagnostics.png")
        v35 = _json(REPO / "docs" / "verification" / "v35_name_class_contrast.json")
        m |= contrast_macros(v35)
        if v35:
            fig_names(v35, primary)
    write_macros(m)
    digest = hashlib.sha256((OUTDIR / "summary_numbers.tex").read_bytes()).hexdigest()[:12]
    print(f"macros: {len(m)} (summary_numbers.tex {digest}); missing: {len(MISSING)}")
    for w in MISSING:
        print(f"  missing: {w}")


if __name__ == "__main__":
    main()
