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
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

REPO = Path(__file__).resolve().parents[4]
OUTDIR = Path(__file__).resolve().parent
DATA = REPO / "data"
MISSING: list[str] = []
# 図の元データ（Source Data、B-2 追記 3）: シート名 → 行（列名 → 値）。図を描く関数が描いた値をそのまま入れる
SOURCE: dict[str, list[dict]] = {}


def write_source_data() -> None:
    """submission/source_data.xlsx: 図・パネルごとに 1 シート。値は図に描いたものと同じ（丸めない）。"""
    if not SOURCE:
        return
    import datetime as _dt

    from openpyxl import Workbook

    wb = Workbook()
    wb.remove(wb.active)
    # 中身が同じなら同じファイルになるように、作成・更新日時を固定する（再生成のたびの差分を避ける）
    wb.properties.created = wb.properties.modified = _dt.datetime(2026, 1, 1, tzinfo=_dt.UTC).replace(tzinfo=None)
    for sheet, rows in sorted(SOURCE.items()):
        ws = wb.create_sheet(sheet)
        cols = list(dict.fromkeys(k for r in rows for k in r))
        ws.append(cols)
        for r in rows:
            ws.append([r.get(c) for c in cols])
    # xlsx は zip。中のファイルの時刻も固定して書き直す
    import io
    import zipfile

    buf = io.BytesIO()
    wb.save(buf)
    src = zipfile.ZipFile(io.BytesIO(buf.getvalue()))
    out = io.BytesIO()
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for info in src.infolist():
            zi = zipfile.ZipInfo(info.filename, date_time=(2026, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            data = src.read(info.filename)
            if info.filename == "docProps/core.xml":  # 保存時に書き込まれる更新日時も固定
                data = re.sub(rb"(<dcterms:modified[^>]*>)[^<]*", rb"\g<1>2026-01-01T00:00:00Z", data)
            z.writestr(zi, data)
    (OUTDIR.parent / "submission" / "source_data.xlsx").write_bytes(out.getvalue())


def _json(path: Path) -> dict | None:
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def missing(what: str) -> str:
    MISSING.append(what)
    safe = what.replace("_", "\\_").replace("->", "$\\to$")  # 本文中でも組めるように
    return f"\\missing{{{safe}}}"


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
                "PowerNullNegative": str(null["decision_counts"].get("negative", 0)),
                "PowerNullLabelList": ", ".join(f"{k} {v}" for k, v in sorted(null["decision_counts"].items()))}
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
    ax.text(4, 44, "Nested comparison within areas (primary, H2)", fontsize=7.2, weight="bold")
    y, h = 30, 9
    xs = [(1, 15, "Position\n(M0)"),
          (21, 22, "+ terrain, hydrology,\nland cover and\nsettlement layout (M1)"),
          (50, 21, "+ regional composition\nof names: which names\noccur in the area (M1c)"),
          (78, 21, "+ the name of each\ncell's nearest\nsettlements (M3n)")]
    for x, w, txt in xs:
        box(x, y, w, h, txt, fc=light if x == 78 else "white", ec=blue if x == 78 else "0.3")
    for (x0, w0, _), (x1, _, _) in pairwise(xs):
        arrow(x0 + w0 + 0.8, y + h / 2, x1 - 0.8, y + h / 2)
    ax.text(18, 27.5, "terrain increment $\\Delta_T$\n(the yardstick)", fontsize=6.2,
            ha="center", va="top", color=grey)
    ax.text(74.5, 27.5, "names increment $\\Delta_N$\n$R=\\Delta_N/\\Delta_T$", fontsize=6.2,
            ha="center", va="top", color=blue)
    box(78, 16.5, 21, 7, "placebo (M3np): the same\nnames reassigned among the\narea's settlements; real "
        "names\nmust do better", ec=grey, ls="--", size=6.0)
    ax.text(1, 13.5, "b", fontsize=9, weight="bold")
    ax.text(4, 13.5, "Terrain-information ladder (H3)", fontsize=7.2, weight="bold")
    rungs = ["T0\nno\nelevation", "T1\n1 km", "T2\n250 m", "T3\n30 m\nstatistics",
             "T4 = M1\n+ HAND"]
    for k, r in enumerate(rungs):
        box(1 + k * 6.8, 2, 5.8, 7.5, r, size=5.4)
    ax.text(18, 11.2, "names added at every rung, within areas", fontsize=6.0, ha="center",
            color=grey)
    ax.text(37, 13.5, "c", fontsize=9, weight="bold")
    ax.text(40, 13.5, "Three readings of names", fontsize=7.2, weight="bold")
    reads = [(("lexicon (N-lex)\nwater, wetland, riverine,\ncrossing, lowland and\nhighland words\nin "
               f"{n_langs} languages; knows\nno places"), 37.5),
             ("character model (N-learn)\nlearned on training\nregions only", 54),
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
    SOURCE["Fig2a"] = [{"aoi": r["tag"], "activation": r["activation"], "lon_centre": centre(r)[0],
                        "lat_centre": centre(r)[1], "evaluable_cells": r["n_eval"], "included": bool(r["included"])}
                       for r in inc["aois"]]
    SOURCE["Fig2b"] = [{"event_year": k, "analysed_aois": yrs[k]} for k in ks]
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


def _forest_rows(primary: dict, extra: dict[str, dict],
                 v36: dict | None = None) -> list[tuple[str, float, list, str]]:
    """(表示名, R̂, 95% CI, 様式) の列。様式 = main（確証）/ secondary（副次）/ undefined（物差しが
    不確か: 地形の増分の 95% CI が 0 を含み、プロトコル上 R は定義されない）。区間は全て固定マージン
    （チャネルと H1 は v36、world_aggregate の比の bootstrap ではなく）。"""
    ch = primary.get("channels", {})

    def fm_row(comp: str, blk: dict | None, rk: str, ck: str) -> tuple[float, list] | None:
        f = fixed_margin(v36, "primary", comp)
        return (f["R"], f["ci95_R"]) if f else _row(blk, rk, ck)

    spec = [
        ("Names within areas (primary, H2)", _row(primary.get("primary"), "R", "ci95_R"), "main"),
        ("   lexicon only (N-lex)", fm_row("channel_lexicon", ch.get("M3nx"), "R", "ci95_R"), "main"),
        ("   learned only (N-learn)", fm_row("channel_learned", ch.get("M3nl"), "R", "ci95_R"), "main"),
        ("   placebo: names permuted", _row(primary.get("placebo"), "R_placebo",
                                             "ci95_R_placebo"), "main"),
        ("Position only, within areas (H1)", fm_row("H1", primary.get("H1"), "ratio_to_terrain",
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
        if met.get("H4_total"):  # 副次の判定で唯一の正（positive-small）も同じ図に（B-2 追記 3）
            spec.append(("   language model incl. its composition (secondary)",
                         _row(met["H4_total"], "R", "ci95_R"), "secondary"))
    spec += [
        ("Total incl. regional composition (secondary)", _row(primary.get("total"), "R",
                                                               "ci95_R"), "secondary"),
        ("Flexible model refitted with names (secondary)", _row(primary.get("flexible"), "R",
                                                                 "ci95_R"), "secondary"),
    ]
    return [(n, v[0], v[1], s) for n, v, s in spec if v is not None]


def fig_main(rows: list[tuple[str, float, list, str]], decision: str) -> None:
    SOURCE["Fig3"] = [{"comparison": n.strip(), "R": r, "ci95_low": c[0] if c else None,
                       "ci95_high": c[1] if c else None, "style": s} for n, r, c, s in rows]
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
    SOURCE["Fig5"] = [{"rung": rg, "R": float(d[k]), "ci95_low": float(c[k, 0]), "ci95_high": float(c[k, 1])}
                      for k, rg in enumerate(ladder["rungs"])]
    ax.axhspan(-SESOI, SESOI, color="0.92", lw=0)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.errorbar(x, d, yerr=[d - c[:, 0], c[:, 1] - d], fmt="o-", color="#1f6fb4", ms=4,
                lw=1.2, capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels(["T0\nno elevation", "T1\n1 km", "T2\n250 m", "T3\n30 m stats",
                        "T4\n+HAND"], fontsize=7)
    ax.set_ylabel("names within areas\n(fraction of terrain increment)", fontsize=8)
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
    if 0 < abs(x) < 5e-5:  # 4 桁でも 0.0000 に丸まる値は指数表記で（0 を挟むかが読めるように）
        e = int(np.floor(np.log10(abs(x))))
        return f"\\ensuremath{{{x / 10**e:.0f}\\times10^{{{e}}}}}"
    return f"\\ensuremath{{{x:.4f}}}" if 0 < abs(x) < 0.0005 else f"\\ensuremath{{{x:.3f}}}"


def _ci_macros(prefix: str, val: float, ci95: list | None) -> dict[str, str]:
    out = {prefix: rv(val)}
    if ci95:
        out |= {f"{prefix}Lo": rv(ci95[0]), f"{prefix}Hi": rv(ci95[1])}
    return out


def fixed_margin(v36: dict | None, split: str, comp: str) -> dict | None:
    """v36（物差しの頑健性）から、登録どおりの物差し Y0 での固定マージンの R と区間・比の bootstrap 区間。
    world_aggregate はチャネル・H1 の区間を比の bootstrap で出しているので、判定の比較とそろえるために使う。"""
    sp = ((v36 or {}).get("splits") or {}).get(split, {}).get("cluster", {})
    row = sp.get(comp)
    if not row:
        return None
    y0 = sp["yardsticks"]["Y0"]
    return {"R": row["delta_bs"] / y0, "ci95_R": [x / y0 for x in row["ci95_delta"]],
            "ci90_R": row["Y0"]["ci90_R"], "ratio95": row["ci95_R_ratio_bootstrap"]}


def detail_macros(m: dict, v36: dict | None = None) -> dict[str, str]:
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
    # チャネルと H1 は固定マージンの区間（v36）。無ければ world_aggregate の比の bootstrap
    for name, comp, blk, rk, ck in (("RLex", "channel_lexicon", ch["M3nx"], "R", "ci95_R"),
                                    ("RLearn", "channel_learned", ch["M3nl"], "R", "ci95_R"),
                                    ("RHOne", "H1", m["H1"], "ratio_to_terrain", "ci95_ratio"),
                                    ("RHOneFlex", "H1_flexible", m["H1_flexible"], "ratio_to_terrain",
                                     "ci95_ratio")):
        f = fixed_margin(v36, "primary", comp)
        out |= _ci_macros(name, f["R"], f["ci95_R"]) if f else _ci_macros(name, blk[rk], blk[ck])
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
    # B-2 追記 1: 辞書の被覆の層（辞書に項目の無い言語の国 = ノルウェーだけ）。言語別の層の「no」と同じセル
    sno = ((v36 or {}).get("splits", {}).get("primary", {}).get("strata") or {}).get("language_no")
    if sno and sno.get("ci90_R"):
        out |= {"RStratNorway": rv(sno["R"]), "RStratNorwayLo": rv(sno["ci90_R"][0]),
                "RStratNorwayHi": rv(sno["ci90_R"][1]), "NCellsStratNorway": _n(int(sno["n_cells"]))}
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
    # 副次比較で判定 label のついたもの（全体の増分・柔軟モデル、llmcut は N-read 版も）の数と positive の数
    sec = [blk for met_ in [extra.get(s) for s in ("country", "temporal", "llmcut")] if met_
           for k in ("total", "flexible", "H4_total", "H4_flexible") if (blk := met_.get(k))]
    out["NSecondaryVerdicts"] = str(len(sec))
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


# %% 事後の解析（amendment B-2: v36–v39 と追加の実行）— 頑健性の段落・図 5・補足表。確証判定は変えない
CLASS_LABELS = {"wetside": "water-side (any)", "water": "water", "wet": "wetland",
                "riverine": "riverine", "crossing": "crossing", "low": "lowland", "high": "highland"}
BAND_LABELS = {"le250m": "\\leq250\\,m", "250_500m": "250--500\\,m", "500m_1km": "0.5--1\\,km",
               "gt1km": ">1\\,km"}


WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
         "eleven", "twelve")


def _word(n: int | str) -> str:
    n = int(n)
    return WORDS[n] if 0 <= n < len(WORDS) else str(n)


def _sgn(x: float, nd: int = 1) -> str:
    """符号つきの数（数式モード、マイナス記号で）。"""
    return f"\\ensuremath{{{x:+.{nd}f}}}".replace("+", "{+}")


def _pp(x: float, nd: int = 1) -> str:
    """百分率ポイント（符号つき）。0 でない値が ±0.0 に丸まるときは桁を足す（区間が 0 を挟むかが読めるように）。"""
    v = 100 * x
    while nd < 3 and v != 0 and round(abs(v), nd) == 0:
        nd += 1
    return f"\\ensuremath{{{v:+.{nd}f}}}".replace("+", "{+}").replace("-", "{-}")


def _pp_ci(prefix: str, r: dict | None, nd: int = 1) -> dict[str, str]:
    if not r:
        return {prefix: missing(prefix)}
    return {prefix: _pp(r["estimate"], nd), f"{prefix}Lo": _pp(r["ci95"][0], nd),
            f"{prefix}Hi": _pp(r["ci95"][1], nd)}


def _or(x: float) -> str:
    return f"{x:.2f}"


def _closer(x: float) -> str:
    """log(1 + 距離) の差 → 何 % 近いか（正 = 近い）。"""
    return f"{100 * (1 - np.exp(x)):.0f}\\%"


def positive_control_macros(v: dict | None) -> dict[str, str]:
    """v37: 正の対照（ラベル非依存）・集落からの距離・距離帯別の主比較・上位 10% の重なり・効果量の翻訳。
    セル単位は B-2 で登録した尺度、集落の点と距離帯は登録外（v37 を書く段階で追加）。"""
    if not v:
        return {"PCWetHand": missing("v37")}
    pc = v["positive_control"]
    st, ce, bands = pc["settlements"], pc["cells"], pc["bands"]
    out = {"NSettlementPoints": _n(pc["n_settlement_points"]),
           "NPCAoi": str(st["wetside"]["low_hand"]["n_aoi"])}

    def mean_pct(r: dict, key: str) -> str:
        return f"{100 * r[key]:.1f}\\%"

    # 集落の点（登録外）
    out |= _pp_ci("PCWetHand", st["wetside"]["low_hand"])
    out |= {"PCWetHandGroupPct": mean_pct(st["wetside"]["low_hand"], "mean_group"),
            "PCWetHandOtherPct": mean_pct(st["wetside"]["low_hand"], "mean_other")}
    out |= _pp_ci("PCWetHalf", st["wetside"]["low_half"])
    out |= _pp_ci("PCHighHand", st["high"]["low_hand"])
    out |= _pp_ci("PCHighHalf", st["high"]["low_half"])
    r = st["wetside"]["hand_mean"]
    out |= {"PCWetHandMean": _sgn(r['estimate'], 1),
            "PCWetHandMeanLo": _sgn(r['ci95'][0], 1),
            "PCWetHandMeanHi": _sgn(r['ci95'][1], 1),
            "PCWetHandMeanGroup": f"{r['mean_group']:.1f}", "PCWetHandMeanOther": f"{r['mean_other']:.1f}"}
    r = st["high"]["hand_mean"]
    out |= {"PCHighHandMean": _sgn(r['estimate'], 1),
            "PCHighHandMeanLo": _sgn(r['ci95'][0], 1),
            "PCHighHandMeanHi": _sgn(r['ci95'][1], 1)}
    r = st["wetside"]["log_river"]
    out |= {"PCWetRiverCloser": _closer(r["estimate"]), "PCWetRiverCloserLo": _closer(r["ci95"][1]),
            "PCWetRiverCloserHi": _closer(r["ci95"][0])}
    out["PCCrossRiverCloser"] = _closer(st["crossing"]["log_river"]["estimate"])
    nr = pc.get("settlements_not_rounded", {}).get("wetside", {}).get("hand_mean")
    if nr:
        out |= {"PCWetHandMeanNR": _sgn(nr['estimate'], 1),
                "PCWetHandMeanNRLo": _sgn(nr['ci95'][0], 1),
                "PCWetHandMeanNRHi": _sgn(nr['ci95'][1], 1)}
    out |= {"ShareCellsRounded": pct(pc["share_cells_minute_rounded"]),
            "ShareSettlementsRounded": pct(pc["share_settlement_points_minute_rounded"])}
    # 全セル（B-2 で登録した尺度）
    out |= _pp_ci("PCWetHandCells", ce["wetside"]["low_hand"])
    out |= {"PCWetHandCellsGroupPct": mean_pct(ce["wetside"]["low_hand"], "mean_group"),
            "PCWetHandCellsOtherPct": mean_pct(ce["wetside"]["low_hand"], "mean_other")}
    out |= _pp_ci("PCWetHalfCells", ce["wetside"]["low_half"])
    r = ce["wetside"]["log_river"]
    out |= {"PCWetRiverCells": _sgn(r['estimate'], 2),
            "PCWetRiverCellsLo": _sgn(r['ci95'][0], 2),
            "PCWetRiverCellsHi": _sgn(r['ci95'][1], 2)}
    r = ce["wetside"]["hand_mean"]
    out |= {"PCWetHandMeanCells": _sgn(r['estimate'], 1),
            "PCWetHandMeanCellsLo": _sgn(r['ci95'][0], 1),
            "PCWetHandMeanCellsHi": _sgn(r['ci95'][1], 1)}
    out |= _pp_ci("PCHighHalfCells", ce["high"]["low_half"])
    # B-2 追記 3: 全セルでの区分別（河川・渡しの名は低い土地を示す。まとめた null は水・低地の名による）
    for c, nm in (("riverine", "River"), ("crossing", "Cross"), ("water", "Water"), ("low", "Low")):
        if ce.get(c) and ce[c].get("low_half"):
            out |= _pp_ci(f"PC{nm}HalfCells", ce[c]["low_half"])
    r = ce["high"]["hand_mean"]
    out |= {"PCHighHandMeanCells": _sgn(r['estimate'], 1),
            "PCHighHandMeanCellsLo": _sgn(r['ci95'][0], 1),
            "PCHighHandMeanCellsHi": _sgn(r['ci95'][1], 1)}
    out |= _pp_ci("PCWetHandNear", bands["le250m"]["low_hand"])
    out |= _pp_ci("PCWetHandFar", bands["gt1km"]["low_hand"])
    d = v["nearest_settlement_distance_m"]
    out |= {"DistQOne": f"{d['0.25'] / 1000:.1f}", "DistMedian": f"{d['0.5'] / 1000:.1f}",
            "DistQThree": f"{d['0.75'] / 1000:.1f}", "ShareWithinKm": pct(v["share_within_1km"]),
            "ShareWetside": pct(v["share_by_group"]["wetside"])}
    pd_ = v.get("primary_by_distance", {})
    for key, name in (("le250m", "Near"), ("le500m", "NearFive"), ("gt1km", "Far")):
        b = pd_.get(key)
        if b:
            out |= {f"R{name}": rv(b["R"]), f"R{name}LoNinety": rv(b["ci90_R"][0]),
                    f"R{name}HiNinety": rv(b["ci90_R"][1]), f"NCells{name}": _n(b["n_cells"])}
    b = v.get("primary_by_distance_wetside", {}).get("le500m")
    if b:
        out |= {"RWetNear": rv(b["R"]), "RWetNearLoNinety": rv(b["ci90_R"][0]),
                "RWetNearHiNinety": rv(b["ci90_R"][1]), "NCellsWetNear": _n(b["n_cells"]),
                "NUnitsWetNear": str(b["n_units"])}
    fc = v["flood_contrast"]
    near = fc["by_band"]["capped"]["le250m"]["y"]
    out |= _pp_ci("ContrastNear", near)
    out |= {"NContrastNearAoi": str(near["n_aoi"]), "NContrastNearPlaces": str(near["n_units"])}
    out |= _pp_ci("ContrastNearResid", fc["by_band"]["capped"]["le250m"]["resid_M1c"])
    out |= _pp_ci("ContrastNearUncapped", fc["by_band"]["uncapped"]["le250m"]["y"])
    out |= _pp_ci("ContrastFarUncapped", fc["by_band"]["uncapped"]["gt1km"]["y"])
    out |= _pp_ci("ContrastFarUncappedResid", fc["by_band"]["uncapped"]["gt1km"]["resid_M1c"])
    out |= _pp_ci("ContrastUncapped", fc["by_group"]["uncapped"]["wetside"]["y"])
    out |= _pp_ci("ContrastUncappedResid", fc["by_group"]["uncapped"]["wetside"]["resid_M1c"])
    out |= _pp_ci("ContrastHandUncapped", fc["hand_yardstick"]["uncapped"])
    out |= _pp_ci("ContrastRiverine", fc["by_group"]["capped"]["riverine"]["y"])
    td = v["top_decile_overlap_M1c_M3n"]
    out |= {"TopDecileOverlap": pct(td["pooled"]), "TopDecileWithin": pct(td["within_aoi_median"])}
    cls = v.get("primary_by_class", {})
    for key, name in (("water", "Water"), ("high", "High"), ("wetside", "Wetside"), ("wet", "Wetland"),
                      ("riverine", "Riverine"), ("crossing", "Crossing"), ("low", "Lowland")):
        b = cls.get(key)
        if b and b.get("ci90_R"):
            out |= {f"RClass{name}": rv(b["R"]), f"RClass{name}LoNinety": rv(b["ci90_R"][0]),
                    f"RClass{name}HiNinety": rv(b["ci90_R"][1]),
                    f"DeltaBSClass{name}": sci(b["delta_bs_names"], 2),
                    f"DeltaBSClass{name}Lo": sci(b["ci90_delta"][0], 2),
                    f"DeltaBSClass{name}Hi": sci(b["ci90_delta"][1], 2),
                    f"NCellsClass{name}": _n(int(b["n_cells"])),
                    f"DeltaBSClass{name}LoNinetyFive": sci(b.get("ci95_delta", [np.nan])[0], 2) if b.get("ci95_delta") else "---",
                    f"DeltaBSClass{name}HiNinetyFive": sci(b.get("ci95_delta", [np.nan, np.nan])[1], 2) if b.get("ci95_delta") else "---",
                    f"RealMinusPermClass{name}": sci(b["real_minus_permuted_delta_bs"], 2) if "real_minus_permuted_delta_bs" in b else "---",
                    f"RealMinusPermClass{name}Lo": sci(b["ci95_real_minus_permuted"][0], 2) if "ci95_real_minus_permuted" in b else "---",
                    f"RealMinusPermClass{name}Hi": sci(b["ci95_real_minus_permuted"][1], 2) if "ci95_real_minus_permuted" in b else "---",
                    f"EqClass{name}": "equivalent" if b["equivalent"] else "not equivalent"}
    eqs = [b["equivalent"] for k, b in cls.items() if k != "wetside" and b.get("ci90_R")]
    if eqs:
        out["NClassEq"] = _word(sum(eqs))
        out["NClassTotal"] = _word(len(eqs))
    nc = v.get("near_settlement_clean")
    if nc:
        out["ShareNearRemoved"] = pct(nc["share_cells_removed_within_500m"])
        for key, name in (("primary", "NearClean"), ("primary_wetside", "WetNearClean")):
            b = nc[key]
            out |= {f"R{name}": rv(b["R"]), f"R{name}LoNinety": rv(b["ci90_R"][0]),
                    f"R{name}HiNinety": rv(b["ci90_R"][1])}
        out |= _pp_ci("ContrastNearClean", nc["flood_contrast"]["capped"]["y"])
    cov = v.get("null_coverage")
    if cov:
        out |= {"NullCoverNinety": pct(cov["coverage90"]), "NullCoverNinetyFive": pct(cov["coverage95"]),
                "NullSdRatio": f"{cov['sd_R_between_reps'] / cov['median_bootstrap_se_R']:.1f}"}
        if "coverage90_zero" in cov:  # B-2 追記 3
            out |= {"NullCoverZeroNinety": pct(cov["coverage90_zero"]),
                    "NullCoverZeroNinetyFive": pct(cov["coverage95_zero"]),
                    "NullRobustSd": f"{cov['robust_sd_R']:.4f}", "NullMedianSe": f"{cov['median_bootstrap_se_R']:.4f}",
                    "NullSdBetween": f"{cov['sd_R_between_reps']:.4f}",
                    "NullStdSd": f"{cov['sd_standardised_deviation']:.2f}", "NullMeanR": rv(cov["mean_R"])}
    e = v.get("effect_size_translation", {})
    if "share_wetside_cells_within_500m" in e:
        out["ShareWetsideWithinFiveHundred"] = pct(e["share_wetside_cells_within_500m"])
    br = e.get("by_run") or {}
    names_run = {"world_v2": "Reg", "world_v2_fs1": "FsOne", "world_v2_fs2": "FsTwo",
                 "world_v2_s6_buffer_b1": "BufferB", "world_v2_actgrp": "ActGrp"}
    fold_or = []
    for run, nm in names_run.items():
        h = (br.get(run) or {}).get("homogeneous") or {}
        if h.get("odds_ratio"):
            out[f"ORUpperRun{nm}"] = f"{h['odds_ratio']:.1f}"
            if run in ("world_v2", "world_v2_fs1", "world_v2_fs2"):
                fold_or.append(h["odds_ratio"])
    if len(fold_or) == 3:  # 確証の fold と付け替え 2 回
        out["ORUpperFoldMin"] = f"{min(fold_or):.1f}"
        out["ORUpperFoldMax"] = f"{max(fold_or):.1f}"
    # 集落の近くだけの効果（B-2 追記 3 の 1）: 格子の外は「e^{δmax} を超える」
    loc_vals = []
    for run, nm in names_run.items():
        lo = (br.get(run) or {}).get("local_500m") or {}
        if lo.get("odds_ratio"):
            out[f"ORLocalUpperRun{nm}"] = f"{lo['odds_ratio']:.1f}"
            loc_vals.append(lo["odds_ratio"])
        elif lo.get("beyond_grid"):
            out[f"ORLocalUpperRun{nm}"] = f"$>${np.exp(lo['delta_max']):.0f}"
            loc_vals.append(float("inf"))
    if loc_vals:
        out["ORLocalUpperRunMin"] = f"{min(loc_vals):.1f}"
        mx = max(loc_vals)
        out["ORLocalUpperRunMax"] = "more than 20" if mx == float("inf") else f"{mx:.1f}"
    c0 = (e.get("curve") or {}).get("0.0")
    if c0:  # 合成ラベルでの地形の増分（δ = 0 の回の平均）。実データの物差しより小さいので翻訳は ΔBS の尺度で
        out["SimYZero"] = sci(c0["mean_yardstick"], 2)
    # δ の格子がそろってから（一様な効果は 7 点、集落の近くだけの効果は 4 点。粗い格子の補間は使わない）
    for kind, prefix, curve, n_grid in (("homogeneous", "", "curve", 7), ("local_500m", "Local", "curve_local", 4)):
        tr = e.get(kind) or {}
        grid_ok = len(e.get(curve, {})) >= n_grid
        for key, name in (("margin", "Margin"), ("achieved_upper90", "Upper"), ("two_way_upper90", "TwUpper")):
            r = tr.get(key)
            if not (r and grid_ok):
                out[f"ORAt{prefix}{name}"] = missing(f"odds ratio {kind} {key}")
                continue
            if r["beyond_grid"]:  # 最大の δ でも届かない: e^{δmax} を超える、とだけ書く
                out[f"ORAt{prefix}{name}"] = f"$>${np.exp(r['delta_max']):.0f}"
                continue
            lo, hi = r["odds_ratio_mc_range"]
            out[f"ORAt{prefix}{name}"] = f"{r['odds_ratio']:.1f}"
            out[f"ORAt{prefix}{name}Lo"] = f"{lo:.1f}" if lo is not None else missing("mc range")
            out[f"ORAt{prefix}{name}Hi"] = f"{hi:.1f}" if hi is not None else missing("mc range")
            out[f"DeltaAt{prefix}{name}"] = f"{r['delta']:.2f}"
    return out


def odds_macros(v: dict | None) -> dict[str, str]:
    """v39: AOI の固定効果つきロジスティック回帰のオッズ比。全セル（距離帯なし）は B-2 で登録、区分 × 距離帯は登録外。"""
    if not v:
        return {"ORWetNearAdj": missing("v39")}
    out = {"NOddsAoi": str(v["wetside_adjusted_M1"]["n_aoi"]), "NOddsGroups": str(v["wetside_adjusted_M1"]["n_groups"]),
           "NOddsDropped": str(len(v["aoi_dropped_no_variation"]))}

    def put(name: str, t: dict) -> None:
        out[name] = _or(t["odds_ratio"])
        out[f"{name}Lo"] = _or(t["ci95_odds_ratio"][0])
        out[f"{name}Hi"] = _or(t["ci95_odds_ratio"][1])

    for key, name in (("wetside_adjusted_M1", "Adj"), ("wetside_raw", "Raw")):
        for band, bname in (("le500m", "Near"), ("500m_1km", "Mid"), ("gt1km", "Far")):
            put(f"ORWet{bname}{name}", v[key]["terms"][f"wetside:{band}"])
        for band, bname in (("500m_1km", "Mid"), ("gt1km", "Far")):
            put(f"ORBand{bname}{name}", v[key]["terms"][f"band_{band}"])
    t = v["categories_adjusted_M1"]["terms"]
    out["NOddsComparisons"] = str(sum(1 for k in t if ":" in k))
    put("ORRiverNearAdj", t["riverine:le500m"])
    # B-2 追記 3: 区分別の値（集落から 500 m 以内は登録外、全セルは登録）
    ta = (v.get("categories_allcells_adjusted_M1") or {}).get("terms", {})
    near = []
    for c, nm in (("water", "Water"), ("wet", "Wetland"), ("riverine", "Riverine"), ("crossing", "Crossing"),
                  ("low", "Lowland"), ("high", "Highland")):
        if f"{c}:le500m" in t:
            put(f"ORClass{nm}NearAdj", t[f"{c}:le500m"])
            near.append(t[f"{c}:le500m"])
        if f"{c}:all" in ta:
            put(f"ORClass{nm}AllAdj", ta[f"{c}:all"])
    if near:  # 6 区分の幅（片側の区分だけを挙げないため）と、区間が 1 を含む区分の数
        out["ORClassNearAdjMin"] = _or(min(x["odds_ratio"] for x in near))
        out["ORClassNearAdjMax"] = _or(max(x["odds_ratio"] for x in near))
        out["NClassNearIncludeOne"] = _word(sum(x["ci95_odds_ratio"][0] <= 1 <= x["ci95_odds_ratio"][1]
                                                for x in near))
        out["NClassNear"] = _word(len(near))
    cb = v.get("cells_by_class_and_band") or {}
    for c, nm in (("water", "Water"), ("wet", "Wetland"), ("riverine", "Riverine"), ("crossing", "Crossing"),
                  ("low", "Lowland")):
        n = (cb.get(c) or {}).get("le500m")
        if n is not None:
            out[f"NCellsClass{nm}Near"] = _n(int(n))
    out["ORTerrainLogit"] = _or(v["wetside_adjusted_M1"]["terms"]["logit_M1"]["odds_ratio"])
    for key, name in (("wetside_allcells_adjusted_M1", "ORWetAllAdj"), ("wetside_allcells_raw", "ORWetAllRaw")):
        if key in v:
            put(name, v[key]["terms"]["wetside:all"])
    if "wetside_clean_adjusted_M1" in v:
        put("ORWetNearCleanAdj", v["wetside_clean_adjusted_M1"]["terms"]["wetside:le500m"])
        out["ShareCleanKept"] = pct(v["clean_share_cells_kept"])
    bc = v.get("bootstrap_check_wetside")
    if bc:
        b = bc["adjusted_M1"]["wetside:le500m"]["odds_ratio_ci95_bootstrap"]
        out |= {"ORWetNearAdjBootLo": _or(b[0]), "ORWetNearAdjBootHi": _or(b[1])}
    return out


def sensor_macros(v: dict | None) -> dict[str, str]:
    """v38: 洪水後の画像のセンサー（レーダー・光学）と、群ごとの主比較。"""
    if not v:
        return {"NAoiSar": missing("v38")}
    c, inc = v["aoi_class_counts"], v["increments"]
    out = {"NAoiSar": str(c["sar"]), "NAoiOptical": str(c["optical"]),
           "NAoiSensorUnknown": str(c["unknown"])}
    gens = v.get("by_schema_generation") or {}
    for g, nm in (("0", "Zero"), ("1", "One"), ("2", "Two"), ("3", "Three")):
        if g in gens:  # 解析した AOI のスキーマ世代ごとの数（B-2 追記 1、v38）
            out[f"NAoiGen{nm}"] = str(sum(gens[g].values()))
    for key, name in (("sar", "Sar"), ("optical", "Optical")):
        b = inc[key]
        out |= {f"R{name}": rv(b["R"]), f"R{name}LoNinety": rv(b["ci90_R"][0]),
                f"R{name}HiNinety": rv(b["ci90_R"][1]), f"NUnits{name}": str(b["n_units"])}
    return out


def robustness_macros(v: dict | None) -> dict[str, str]:
    """v36: 物差し 5 通り × bootstrap 2 通りの同等性、地域構成の増分、Holm、決着しなかった比較、
    activation で束ねた分割・中心化した名前の特徴・N-read の全体の増分と地域構成（B-2 追記 2）。"""
    if not v:
        return {"NYardEqPrimary": missing("v36")}
    sp = v["splits"]
    out: dict[str, str] = {}

    def eq_count(key: str, comp: str) -> tuple[int, int]:
        n = k = 0
        r = sp.get(key)
        if not r:
            return 0, 0
        for kind in ("cluster", "two_way"):
            row = r[kind].get(comp)
            if not row:
                continue
            for y in r[kind]["yardsticks"]:
                n += 1
                k += bool(row[y]["equivalent"])
        return k, n

    def def_count(key: str, comp: str) -> int:
        """物差しが定義された（物差しが正で R の区間がある）組の数。"""
        r = sp.get(key) or {}
        return sum(bool((r.get(kind, {}).get(comp) or {}).get(y, {}).get("ci90_R"))
                   for kind in ("cluster", "two_way") for y in (r.get(kind) or {}).get("yardsticks", {}))

    k, n = eq_count("primary", "primary")
    out |= {"NYardEqPrimary": str(k), "NYardTotal": str(n), "NYardTotalWord": _word(n),
            "NYardEqPrimaryWord": _word(k)}
    for key, comp, name in (("country", "primary", "Country"), ("llmcut", "H4", "HFour"),
                            ("primary_actgrp", "primary", "ActGrp"), ("primary_dm", "demeaned", "Dm"),
                            ("primary_dm_actgrp", "demeaned", "DmActGrp"), ("country_dm", "demeaned", "DmCountry"),
                            ("temporal_dm", "demeaned", "DmTemporal"), ("llmcut_dm", "demeaned", "DmLlm")):
        kk, nn = eq_count(key, comp)
        if nn:
            out[f"NYardEq{name}"] = str(kk)
            out[f"NYardEq{name}Word"] = _word(kk)
            nd = def_count(key, comp)
            out[f"NYardDef{name}"] = str(nd)
            out[f"NYardDef{name}Word"] = _word(nd)
        else:
            out[f"NYardEq{name}"] = missing(f"yardsticks {key}")
            out[f"NYardEq{name}Word"] = missing(f"yardsticks {key}")

    def r90(key: str, kind: str, comp: str, y: str, name: str) -> dict[str, str]:
        r = sp.get(key)
        e = ((r or {}).get(kind, {}).get(comp) or {}).get(y)
        if not e or not e["ci90_R"]:
            return {f"{name}Lo": "---", f"{name}Hi": "---"}
        return {f"{name}": rv(e["R"]), f"{name}Lo": rv(e["ci90_R"][0]), f"{name}Hi": rv(e["ci90_R"][1])}

    out |= r90("primary", "two_way", "primary", "Y1", "RPrimTwoWayYOne")
    out |= r90("primary", "two_way", "primary", "Y0", "RPrimTwoWay")
    out |= r90("country", "cluster", "primary", "Y1", "RCountryYOne")
    out |= r90("llmcut", "cluster", "H4", "Y1", "RHFourYOne")
    out |= r90("llmcut", "cluster", "H4", "Y3", "RHFourYThree")
    out |= r90("primary_actgrp", "cluster", "primary", "Y3", "RActGrpYThree")
    out |= r90("primary_actgrp", "cluster", "primary", "Y2", "RActGrpYTwo")
    for y in ("Y2", "Y3", "Y4"):
        nm = {"Y2": "YTwo", "Y3": "YThree", "Y4": "YFour"}[y]
        out |= r90("llmcut", "cluster", "H4_total", y, f"RHFourTotal{nm}")
    tw = (sp.get("llmcut", {}).get("two_way", {}).get("H4_total") or {})
    if tw:
        out |= {"DeltaHFourTotalTwLo": sci(tw["ci90_delta"][0], 2), "DeltaHFourTotalTwHi": sci(tw["ci90_delta"][1], 2)}
    out |= {"YOnePrimary": sci(sp["primary"]["cluster"]["yardsticks"]["Y1"], 2),
            "YThreePrimary": sci(sp["primary"]["cluster"]["yardsticks"]["Y3"], 2),
            "YFourPrimary": sci(sp["primary"]["cluster"]["yardsticks"]["Y4"], 2)}
    out["YTwoPrimary"] = (sci(sp["primary"]["cluster"]["yardsticks"]["Y2"], 2)
                          if "Y2" in sp["primary"]["cluster"]["yardsticks"] else missing("Y2 (recal)"))
    tw = sp["primary"]["two_way_axes"]
    out |= {"NTwoWayActivations": str(tw["n_activations"]), "NTwoWayCountries": str(tw["n_countries"]),
            "NTwoWayActivationsLlm": str(sp["llmcut"]["two_way_axes"]["n_activations"])}
    # 地域構成の増分（再較正だけの基準から）。区間は ΔBS の 95% を登録の物差しで割ったもの
    for key in (*SPLITS_ALL, "primary_actgrp"):
        r = sp.get(key)
        if not r:
            continue
        comp = r["cluster"]["composition"]
        y0 = r["cluster"]["yardsticks"]["Y0"]
        name = "ActGrp" if key == "primary_actgrp" else key.capitalize()
        for ck, mname in (("M0s->M0c", "CompPosRecal"), ("M1s->M1c", "CompTerrain"),
                          ("M1s->M1rc", "CompLlm")):
            c = comp.get(ck)
            if c:
                out[f"R{mname}{name}"] = rv(c["relative_to_Y0"])
                out[f"R{mname}{name}Lo"] = rv(c["ci95_delta"][0] / y0)
                out[f"R{mname}{name}Hi"] = rv(c["ci95_delta"][1] / y0)
                if ck == "M1s->M1rc":  # 定数からの物差しでも（指摘 B4）
                    y3 = r["cluster"]["yardsticks"]["Y3"]
                    out["RCompLlmYThree"] = rv(c["delta_bs"] / y3)
                    out["DeltaCompLlm"] = sci(c["delta_bs"], 2)
                    out["DeltaCompLlmLo"] = sci(c["ci95_delta"][0], 2)
                    out["DeltaCompLlmHi"] = sci(c["ci95_delta"][1], 2)
                    tw = (r.get("two_way", {}).get("composition") or {}).get("M1s->M1rc")
                    if tw:
                        out["DeltaCompLlmTwLo"] = sci(tw["ci95_delta"][0], 2)
                        out["DeltaCompLlmTwHi"] = sci(tw["ci95_delta"][1], 2)
            elif ck == "M0s->M0c":  # 再較正の基準は全分割にあるはず
                out.setdefault(f"R{mname}{name}", missing(f"{ck} {key}"))
    ag = sp.get("primary_actgrp")
    if ag:  # 束ねた分割で地域構成が基準を悪くした大きさ（ΔBS とその区間）
        c = ag["cluster"]["composition"].get("M1s->M1c")
        if c:
            out |= {"DeltaCompActGrp": sci(c["delta_bs"], 2), "DeltaCompActGrpLo": sci(c["ci95_delta"][0], 2),
                    "DeltaCompActGrpHi": sci(c["ci95_delta"][1], 2)}
    h = v["holm_H3_H4"]["adjusted"]
    out["HolmHFour"] = "<0.001" if h.get("H4_TOST", 1) < 0.001 else f"{h['H4_TOST']:.3f}"
    out["HolmHThree"] = f"{h['H3_delta0_gt_delta4']:.2f}"
    for k, name in (("primary/total", "Total"), ("primary/flexible", "Flex")):
        inc = v["inconclusive"].get(k)
        if inc and inc["units_needed_to_exclude_SESOI"]:
            out[f"NPlacesNeeded{name}"] = _n(int(round(inc["units_needed_to_exclude_SESOI"], -2)))
    sp_runs = v.get("spread_primary", {})
    # B-2 追記 3: 国外し分割の real − permuted（下限が 0 を僅かに上回る）、fold の付け替えによる R の幅と SE の比
    cp = ((sp.get("country") or {}).get("cluster") or {}).get("primary") or {}
    if "R_minus_Rp_Y0" in cp:
        out |= {"RCountryRealMinusPlacebo": rv(cp["R_minus_Rp_Y0"]),
                "RCountryRealMinusPlaceboLo": rv(cp["ci95_R_minus_Rp_Y0"][0]),
                "RCountryRealMinusPlaceboHi": rv(cp["ci95_R_minus_Rp_Y0"][1])}
    draws = [sp_runs.get(r) for r in ("world_v2", "world_v2_fs1", "world_v2_fs2")]
    if all(draws):
        se = (draws[0]["ci90_R"][1] - draws[0]["ci90_R"][0]) / (2 * 1.645)
        rs = [d["R"] for d in draws]
        out["FoldRangeOverSe"] = f"{(max(rs) - min(rs)) / se:.1f}"
        sd = float(np.std(rs, ddof=1))
        half = 1.645 * float(np.sqrt(se ** 2 + sd ** 2))
        out |= {"RCombinedLoNinety": rv(draws[0]["R"] - half), "RCombinedHiNinety": rv(draws[0]["R"] + half),
                "RFoldSd": f"{sd:.4f}", "RBootSe": f"{se:.4f}"}
    for run, name in (("world_v2_fs1", "FoldSeedOne"), ("world_v2_fs2", "FoldSeedTwo"),
                      ("world_v2_actgrp", "ActGrp")):
        r = sp_runs.get(run)
        out |= ({f"R{name}": rv(r["R"]), f"R{name}LoNinety": rv(r["ci90_R"][0]),
                 f"R{name}HiNinety": rv(r["ci90_R"][1]), f"Verdict{name}": r["decision"],
                 f"YZero{name}": sci(r["delta_bs_M0_M1"], 2)}
                if r else {f"{k}{name}{s}": missing(run) for k, s in
                           (("R", ""), ("R", "LoNinety"), ("R", "HiNinety"), ("Verdict", ""), ("YZero", ""))})
    ag = sp_runs.get("world_v2_actgrp")
    if ag:
        out |= {"NUnitsActGrp": str(ag["n_units"]),
                "BrierActGrpMOnes": f"{ag['brier']['M1s']:.4f}", "BrierActGrpMOnec": f"{ag['brier']['M1c']:.4f}",
                "BrierActGrpMThreen": f"{ag['brier']['M3n']:.4f}", "BrierActGrpMZero": f"{ag['brier']['M0']:.4f}",
                "RCompTerrainActGrp": rv((ag["brier"]["M1s"] - ag["brier"]["M1c"]) / ag["delta_bs_M0_M1"])}
        c = ag.get("ci95_R_minus_Rp")
        if c:
            out |= {"RActGrpRealMinusPlaceboLo": rv(c[0]), "RActGrpRealMinusPlaceboHi": rv(c[1]),
                    "RActGrpRealMinusPlacebo": rv(ag["R_minus_Rp"])}
    agc = (sp.get("primary_actgrp") or {}).get("brier", {})
    if "const" in agc:
        out["BrierActGrpConst"] = f"{agc['const']:.4f}"
    for key, name in (("primary", "Primary"), ("primary_actgrp", "ActGrp"), ("country", "Country"),
                      ("temporal", "Temporal"), ("llmcut", "Llm")):
        pt = (sp.get(key) or {}).get("prespecified_two_way") or {}
        if pt.get("ci95_R"):
            out[f"RPreTwoWay{name}Lo"] = rv(pt["ci95_R"][0])
            out[f"RPreTwoWay{name}Hi"] = rv(pt["ci95_R"][1])
    for key, name in (("country", "Country"), ("llmcut", "Llm"), ("primary", "Primary")):
        ys = ((sp.get(key) or {}).get("two_way") or {}).get("yardsticks") or {}
        if "Y1" in ys:  # 二方向 bootstrap での地形の増分の 95% 下限（負なら、その bootstrap では R が未定義）
            out[f"YOneTwoWay{name}"] = sci(ys["Y1"], 2)
    pt = (sp.get("primary") or {}).get("prespecified_two_way") or {}
    if pt.get("n_countries"):  # 登録どおりの二方向 bootstrap はセルの最寄り集落の国で数える
        out["NCellCountriesPreTwoWay"] = str(pt["n_countries"])
    # AOI 内で中心化した名前の特徴（M1s → M3d）
    for key, name in (("primary_dm", "Dm"), ("primary_dm_actgrp", "DmActGrp"), ("country_dm", "DmCountry"),
                      ("llmcut_dm", "DmLlm"), ("temporal_dm", "DmTemporal")):
        r = sp.get(key)
        row = (r or {}).get("cluster", {}).get("demeaned")
        if row:
            e = row["Y0"]
            out |= {f"R{name}": rv(e["R"]), f"R{name}LoNinety": rv(e["ci90_R"][0]),
                    f"R{name}HiNinety": rv(e["ci90_R"][1]),
                    # Brier の尺度（時間外しは登録の規則で R が未定義なので、こちらで書く）
                    f"DeltaBS{name}": sci(row["delta_bs"], 2), f"DeltaBS{name}Lo": sci(row["ci95_delta"][0], 2),
                    f"DeltaBS{name}Hi": sci(row["ci95_delta"][1], 2)}
            pl = r["cluster"].get("demeaned_placebo")
            if pl:
                out[f"R{name}Placebo"] = rv(pl["Y0"]["R"])
            if row.get("decision_Y0"):  # 凍結した判定表の label（B-2 追記 3）
                out[f"Verdict{name}"] = row["decision_Y0"]
                out[f"R{name}MinusRp"] = rv(row["R_minus_Rp_Y0"])
                out[f"R{name}MinusRpLo"] = rv(row["ci95_R_minus_Rp_Y0"][0])
                out[f"R{name}MinusRpHi"] = rv(row["ci95_R_minus_Rp_Y0"][1])
        else:
            out[f"R{name}"] = missing(f"demeaned {key}")
    cl = v.get("country_leakage")
    if cl:
        n = cl["test_cells_nearest_settlement_in_training_country"]
        out |= {"NCountryCrossCells": _n(n), "ShareCountryCross": pct(n / cl["n_cells"]),
                "NCountryActSpan": str(cl["activations_spanning_folds"])}
    return out


SPLITS_ALL = ("primary", "country", "temporal", "llmcut")


def hgb_macros(v: dict | None) -> dict[str, str]:
    """v40: 早期終了の反復回数（主分割 fold 0、5 シード）。"""
    if not v:
        return {"HGBItersMZero": missing("v40")}
    m0, m1 = v["arms"]["M0"]["n_iter"], v["arms"]["M1"]["n_iter"]
    out = {"HGBItersMOneMin": str(min(m1)), "HGBItersMOneMax": str(max(m1)),
           "HGBTrainCells": _n(v["n_train"])}
    out["HGBItersMZero"] = (f"all {max(m0)}" if min(m0) == max(m0) else f"{min(m0)}--{max(m0)}")
    return out


def power_split_macros(runs: Path) -> dict[str, str]:
    """副分割の検出力（事後、B-2）: δ = 0 で同等と判定された割合、δ = 1.2 で正と判定された割合（Wilson 区間）。"""
    out: dict[str, str] = {}
    for k, name in (("c", "Country"), ("t", "Temporal"), ("l", "Llm")):
        nul = _json(runs / f"power_v3{k}_null" / "power_summary.json")
        ses = _json(runs / f"power_v3{k}_sesoi" / "power_summary.json")
        # 20 回なので整数の百分率で
        out[f"PowerEquiv{name}"] = (f"{100 * nul['rate_equivalent_TOST']:.0f}\\%" if nul
                                    else missing(f"power {k} null"))
        out[f"PowerPos{name}"] = f"{100 * ses['rate_positive']:.0f}\\%" if ses else missing(f"power {k} sesoi")
        if nul:
            out[f"PowerNReps{name}"] = str(nul["n_reps"])
            out[f"PowerUndefNull{name}"] = str(nul["decision_counts"].get("undefined-denominator", 0))
            out[f"PowerNullLabelNull{name}"] = str(nul["decision_counts"].get("equivalent-null", 0))
        if ses:  # δ = 1.2 は主分割で R ≈ 0.05 に合わせた値。他の分割での平均 R（余白との対応がずれる量）
            out[f"PowerMeanR{name}"] = rv(ses["mean_R"])
            out[f"PowerNRepsSesoi{name}"] = str(ses["n_reps"])
            # δ = 1.2 でも同等と判定された割合と equivalent-null の数（B-2 追記 2 の 7、判定の分布）
            out[f"PowerEquivSesoi{name}"] = f"{100 * ses['rate_equivalent_TOST']:.0f}\\%"
            out[f"PowerNullLabelSesoi{name}"] = str(ses["decision_counts"].get("equivalent-null", 0))
    # B-2 追記 4: 実データの余白で δ = 1.2 の回が同等と判定される割合
    for k, split, name in (("", "primary", "Primary"), ("c", "country", "Country"), ("l", "llmcut", "Llm")):
        rm = real_margin_rate(runs, k, split)
        if rm:
            lo, hi = _wilson(rm[0] / rm[1], rm[1])
            out |= {f"PowerRealEq{name}": f"{rm[0]}", f"PowerRealN{name}": f"{rm[1]}",
                    f"PowerRealEqPct{name}": f"{100 * rm[0] / rm[1]:.0f}\\%",
                    f"PowerRealEqHi{name}": f"{100 * hi:.0f}\\%", f"PowerRealEqLo{name}": f"{100 * max(lo, 0):.0f}\\%",
                    f"PowerRealNonpos{name}": str(rm[2])}
    # B-2 追記 3: 各分割の余白（0.05 × 物差し）と名前の ΔBS の 95% 上限を、δ = 1.2 の平均 ΔBS と比べる
    for k, split, name in (("", "primary", "Primary"), ("c", "country", "Country"), ("t", "temporal", "Temporal"),
                           ("l", "llmcut", "Llm")):
        r = split_ratios(runs, k, split)
        if not r:
            continue
        out[f"PowerSimDBS{name}"] = sci(r["sim"], 2)
        out[f"PowerUpperRatio{name}"] = f"{r['upper']:.2f}"
        if r["margin"] is not None:
            out[f"PowerMarginRatio{name}"] = f"{r['margin']:.1f}"
            out[f"MarginDBS{name}"] = sci(r["margin_dbs"], 2)
    return out


def split_ratios(runs: Path, k: str, split: str) -> dict | None:
    """B-2 追記 3: 実データの余白と観測の ΔBS の 95% 上限を、δ = 1.2 の合成ラベルの平均 ΔBS で割った比。"""
    ses = _json(runs / f"power_v3{k}_sesoi" / "power_summary.json")
    met = _json(runs / "world_v2" / f"metrics_{split}.json")
    reps = (ses or {}).get("reps") or []
    if not (reps and met and "delta_bs" in reps[0]):
        return None
    sim = float(np.mean([r["delta_bs"] for r in reps]))
    pm = met["primary"]
    defined = pm["decision"] != "undefined-denominator"
    return {"sim": sim, "upper": pm["ci95_delta"][1] / sim,
            "margin": pm["margin_delta_bs"] / sim if defined else None,
            "margin_dbs": pm["margin_delta_bs"] if defined else None}


def _fig4_source(v37: dict, v39: dict | None, classes: list[str]) -> None:
    """図 4 の元データ（Source Data）。a: 正の対照、b: 距離帯、c: オッズ比、d: 合成ラベルの曲線と目標。"""
    pc = v37["positive_control"]
    rows = []
    for c in classes:
        for scale in ("settlements", "cells"):
            r = (pc.get(scale, {}).get(c) or {}).get("low_hand")
            if r:
                rows.append({"class": c, "scale": scale, "difference_pp": 100 * r["estimate"],
                             "ci95_low_pp": 100 * r["ci95"][0], "ci95_high_pp": 100 * r["ci95"][1]})
    SOURCE["Fig4a"] = rows
    rows = []
    for b in BAND_LABELS:
        for series, r in (("HAND<=2m (no labels)", pc["bands"][b]["low_hand"]),
                          ("flooded, observed", v37["flood_contrast"]["by_band"]["capped"][b]["y"]),
                          ("flooded, net of terrain model", v37["flood_contrast"]["by_band"]["capped"][b]["resid_M1c"])):
            rows.append({"distance_band": b, "series": series, "difference_pp": 100 * r["estimate"],
                         "ci95_low_pp": 100 * r["ci95"][0], "ci95_high_pp": 100 * r["ci95"][1]})
    SOURCE["Fig4b"] = rows
    rows = []
    if v39:
        for c in classes:
            key = "wetside" if c == "wetside" else "categories"
            for adj in ("adjusted_M1", "raw"):
                tt = v39[f"{key}_{adj}"]["terms"].get(f"{c}:le500m")
                if tt:
                    rows.append({"class": c, "model": "given terrain model" if adj == "adjusted_M1" else "unadjusted",
                                 "odds_ratio": tt["odds_ratio"], "ci95_low": tt["ci95_odds_ratio"][0],
                                 "ci95_high": tt["ci95_odds_ratio"][1]})
    SOURCE["Fig4c"] = rows
    e = v37.get("effect_size_translation", {})
    rows = []
    for key, lab in (("curve", "effect at every distance"), ("curve_local", "effect within 500 m only")):
        for d, cv in sorted((float(k), v) for k, v in e.get(key, {}).items()):
            rows.append({"simulation": lab, "delta": d, "odds_ratio": float(np.exp(d)), "replicates": cv["n_reps"],
                         "mean_delta_bs": cv["mean_delta_bs"], "mc95_low": cv["mc_ci95_mean"][0],
                         "mc95_high": cv["mc_ci95_mean"][1]})
    for k, v in e.get("targets", {}).items():
        rows.append({"simulation": f"target: {k}", "mean_delta_bs": v})
    SOURCE["Fig4d"] = rows


def fig_names(v37: dict, v39: dict | None) -> None:
    """図 5: 名前が記録しているもの。a 集落の点での地形（ラベル非依存）、b 集落からの距離による減衰、
    c 集落から 500 m 以内の浸水のオッズ比（AOI 固定効果）、d 効果量の翻訳（合成ラベル）。"""
    fig = plt.figure(figsize=(7.2, 5.4), layout="constrained")
    axs = fig.subplots(2, 2).ravel()
    classes = list(CLASS_LABELS)
    ys = np.arange(len(classes))[::-1]
    blue, brown, grey = "#1f6fb4", "#8c5a2b", "0.45"
    # a: 集落の点（集落を含むセル）と全セルでの HAND ≤ 2 m の割合の差
    ax = axs[0]
    pc = v37["positive_control"]
    _fig4_source(v37, v39, classes)
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
    ax.set_yticklabels([CLASS_LABELS[c] for c in classes], fontsize=7)
    ax.set_xlabel("difference in share of cells containing ground\nwith HAND ≤ 2 m (pct. points; filled: at the "
                  "settlement, open: all cells)", fontsize=7)
    # b: 水系側の名と他の名の差を、最寄りの集落までの距離帯ごとに
    ax = axs[1]
    bands = list(BAND_LABELS)
    x = np.arange(len(bands))
    series = [("HAND ≤ 2 m (no labels)", [pc["bands"][b]["low_hand"] for b in bands], grey, "s", True),
              ("flooded, observed", [v37["flood_contrast"]["by_band"]["capped"][b]["y"] for b in bands],
               blue, "o", True),
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
    # c: 集落から 500 m 以内の浸水のオッズ比（AOI 固定効果、地形モデルで調整 / 調整なし）
    ax = axs[2]
    if v39:
        for y, c in zip(ys, classes, strict=True):
            col = brown if c == "high" else blue
            key = "wetside" if c == "wetside" else "categories"
            for adj, dy, filled in (("adjusted_M1", 0.15, True), ("raw", -0.15, False)):
                t = v39[f"{key}_{adj}"]["terms"][f"{c}:le500m"]
                ax.plot(t["ci95_odds_ratio"], [y + dy, y + dy], color=col if filled else grey, lw=1.4)
                ax.plot([t["odds_ratio"]], [y + dy], "o", ms=4, color=col if filled else grey,
                        mfc=col if filled else "white")
        ax.set_xscale("log")
        ax.axvline(1, color="0.4", lw=0.8)
        ax.set_xlim(0.3, 4.5)
        ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4]))
        ax.xaxis.set_major_formatter(FixedFormatter(["0.5", "1", "2", "4"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_yticks(ys)
        ax.set_yticklabels([CLASS_LABELS[c] for c in classes], fontsize=7)
        ax.set_xlabel("odds ratio of flooding within 500 m of the settlement\n(filled: given the "
                      "terrain model, open: unadjusted)", fontsize=7)
    # d: 合成ラベルでの名前効果 δ と平均 ΔBS(M1c → M3n)。一様な効果と、集落から 500 m 以内だけの効果
    ax = axs[3]
    e = v37.get("effect_size_translation", {})
    xmax = 0.0
    for key, lab, col, mk in (("curve", "effect at every distance", blue, "o"),
                              ("curve_local", "effect within 500 m only", brown, "s")):
        curve = sorted((float(k), c) for k, c in e.get(key, {}).items())
        if not curve:
            continue
        dx = np.array([k for k, _ in curve])
        mu = 1e4 * np.array([c["mean_delta_bs"] for _, c in curve])
        lo = 1e4 * np.array([c["mc_ci95_mean"][0] for _, c in curve])
        hi = 1e4 * np.array([c["mc_ci95_mean"][1] for _, c in curve])
        ax.errorbar(dx, mu, yerr=[mu - lo, hi - mu], fmt=mk, ms=3.5, color=col, lw=1.2, capsize=0,
                    label=lab)
        # 点の間は、δ の逆算（v37 の invert）と同じくデータの尺度で線形補間した線を描く
        xs = np.linspace(dx.min(), dx.max(), 400)
        ax.plot(xs, np.interp(xs, dx, mu), color=col, lw=1.0)
        xmax = max(xmax, float(dx.max()))
    tg = e.get("targets", {})
    if tg:
        # 線の説明は曲線の通らない帯（一様な効果の曲線の右、局所的な効果の曲線の上）に置く
        for key, lab, ls, va, x0 in (("margin", "equivalence margin", "--", "bottom", 1.75),
                                     ("two_way_upper90", "90% bound (two-way)", "-.", "bottom", 0.66),
                                     ("achieved_upper90", "90% bound (places)", ":", "top", 0.66)):
            if key in tg:
                ax.axhline(1e4 * tg[key], color="0.3", ls=ls, lw=0.9)
                ax.text(x0, 1e4 * tg[key], lab, fontsize=5.5, va=va, ha="left", color="0.25")
        # 観測された上限（0.2–0.4×10^-4）と余裕幅（3.8×10^-4）の両方で曲線との交点が読めるように、0 の近くは線形・上は対数
        ax.set_yscale("symlog", linthresh=0.1, linscale=0.6)
        ax.yaxis.set_major_locator(FixedLocator([0, 0.1, 0.3, 1, 3, 10]))
        ax.yaxis.set_major_formatter(FixedFormatter(["0", "0.1", "0.3", "1", "3", "10"]))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_ylim(-0.05, 1e4 * 4 * tg.get("margin", 1e-4))
        ax.set_xlabel("simulated name effect δ (log-odds; odds ratio e$^δ$ on top)", fontsize=7)
        ax.set_ylabel("ΔBS from names within areas\n(×10$^{-4}$, mean, 95% MC interval)", fontsize=7)
        ax.legend(fontsize=6, frameon=False, loc="upper left")
        top = ax.secondary_xaxis("top", functions=(np.exp, lambda v: np.log(np.clip(v, 1e-9, None))))
        top.set_xticks([t for t in (1, 2, 4, 8, 16) if np.log(t) <= xmax + 1e-9])
        top.tick_params(labelsize=6.5)
    for a, lab in zip(axs, "abcd", strict=True):
        a.tick_params(labelsize=6.5)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.set_title(lab, fontsize=9, fontweight="bold", loc="left")
    fig.savefig(OUTDIR / "fig_names.pdf")
    plt.close(fig)


def si_yardsticks(v: dict) -> str:
    """補足表: 物差し 5 通り × bootstrap 2 通りの固定マージン TOST（事後、v36）。崩れた組も全て載せる。

    行が多く横置きの 1 頁に収まらないので、確証の比較（H2・H4）と、地域内で中心化した比較（M3d）の 2 枚に分ける。
    """
    sp = v["splits"]
    names = {"Y0": "pre-specified", "Y1": "lower 95\\% bound", "Y2": "recalibrated",
             "Y3": "from a constant", "Y4": "stacked"}
    ys = list(names)
    head = "Split, comparison & Bootstrap & " + " & ".join(f"{k}: {n}" for k, n in names.items())
    spec_a = [("primary", "primary", "Location groups, H2"),
              ("primary_actgrp", "primary", "Flood events grouped, H2$^{e}$"),
              ("country", "primary", "Countries, H2"),
              ("temporal", "primary", "Time, H2$^{u}$"),
              ("llmcut", "primary", "LM cut-off, H2"),
              ("llmcut", "H4", "LM cut-off, H4"),
              ("llmcut", "H4_total", "LM cut-off, H4 total$^{s}$")]
    spec_b = [("primary_dm", "demeaned", "Location groups, demeaned"),
              ("primary_dm_actgrp", "demeaned", "Flood events grouped, demeaned$^{e}$"),
              ("country_dm", "demeaned", "Countries, demeaned"),
              ("temporal_dm", "demeaned", "Time, demeaned$^{u}$"),
              ("llmcut_dm", "demeaned", "LM cut-off, demeaned")]

    def rows_for(spec: list[tuple[str, str, str]]) -> list[str]:
        rows = []
        for key, comp, lab in spec:
            r = sp.get(key)
            if not r or comp not in r["cluster"]:
                continue
            for kind, kname in (("cluster", "cluster"), ("two_way", "two-way")):
                row = r[kind].get(comp)
                ysv = r[kind]["yardsticks"]
                cells = []
                for y in ys:
                    e = (row or {}).get(y)
                    if e is None or y not in ysv:
                        cells.append("---")
                    elif not e["ci90_R"]:
                        cells.append("undefined")
                    else:
                        mark = "" if e["equivalent"] else "$^{\\dagger}$"
                        cells.append(f"{_ci(e['ci90_R'])}{mark}")
                rows.append(f"{lab if kind == 'cluster' else ''} & {kname} & " + " & ".join(cells) + " \\\\")
            for kind, kname in (("cluster", "cluster"), ("two_way", "two-way")):  # Y1 は bootstrap ごとに違う
                yv = r[kind]["yardsticks"]
                rows.append(f"\\quad yardstick $\\Delta$BS & {kname} & " + " & ".join(
                    sci(yv[y], 2) if y in yv else "---" for y in ys) + " \\\\")
        return rows

    defs = ("Each cell is the 90\\% interval of $\\Delta\\mathrm{BS}$ divided by the yardstick; "
            "equivalence holds when it lies within $\\pm0.05$ ($^{\\dagger}$: does not hold). Y0: the "
            "pre-specified $\\widehat{\\Delta\\mathrm{BS}}(\\mathrm{M0}\\to\\mathrm{M1})$; Y1: the lower end "
            "of its 95\\% interval under the same bootstrap; Y2: $\\widehat{\\Delta\\mathrm{BS}}(\\mathrm{M0s}"
            "\\to\\mathrm{M1s})$, both base models recalibrated; Y3: $\\widehat{\\Delta\\mathrm{BS}}$ from a "
            "constant prediction (the training-side flood rate; for the single-fold time splits the "
            "evaluation-side rate, which favours the constant) to M1s; Y4: $\\widehat{\\Delta\\mathrm{BS}}"
            "(\\mathrm{M0c}\\to\\mathrm{M1c})$. Cluster: location groups (countries for the country hold-out), "
            "4,000 replicates; two-way: activations $\\times$ countries (each AOI counted in the country of "
            "most of its nearest settlements), 4,000 replicates. Undefined: the yardstick is not positive. "
            "$^{u}$The pre-specified rule leaves $R$ undefined on this split because the 95\\% interval of the "
            "terrain increment includes zero. $^{e}$All areas of each flood event on the same side of every split.")
    a = _table("p{3.3cm}l" + "c" * len(ys), head, rows_for(spec_a),
               "Equivalence of the within-area increment under alternative yardsticks (post hoc, amendment "
               "B-2): the pre-specified comparisons. " + defs + " $^{s}$Secondary comparison (total increment of "
               "the language-model reading over M1s). The pre-specified two-way bootstrap of the ratio is in "
               "Supplementary Table~\\ref{tab:si-registered}; the comparison with name features centred within "
               "each AOI is in Supplementary Table~\\ref{tab:si-yardsticks-b}.",
               "tab:si-yardsticks", landscape=True)
    b = _table("p{3.3cm}l" + "c" * len(ys), head, rows_for(spec_b),
               "As Supplementary Table~\\ref{tab:si-yardsticks}, for the name features centred within each AOI "
               "and added to the recalibrated terrain model (M1s$\\to$M3d), without the regional composition "
               "(post hoc, amendment B-2). " + defs,
               "tab:si-yardsticks-b", landscape=True)
    return a + "\n" + b


def si_registered(v: dict) -> str:
    """補足表: 登録した量のうち他の表に無いもの（B-2 追記 3）。凍結した判定表の label（物差し Y0）、比の bootstrap、
    登録どおりの二方向 bootstrap（比・セルの最寄り集落の国・1,000 回）、同等性が成り立つ最小の物差し、real − permuted。"""
    sp = v["splits"]
    spec = [("primary", "primary", "Location groups, H2"),
            ("country", "primary", "Countries, H2"),
            ("temporal", "primary", "Time, H2$^{u}$"),
            ("llmcut", "primary", "LM cut-off, H2"),
            ("llmcut", "H4", "LM cut-off, H4"),
            ("primary_actgrp", "primary", "Flood events grouped, H2$^{e}$"),
            ("primary_dm", "demeaned", "Location groups, centred$^{d}$"),
            ("primary_dm_actgrp", "demeaned", "Flood events grouped, centred$^{d,e}$"),
            ("country_dm", "demeaned", "Countries, centred$^{d}$"),
            ("temporal_dm", "demeaned", "Time, centred$^{d,u}$"),
            ("llmcut_dm", "demeaned", "LM cut-off, centred$^{d}$")]
    rows = []
    for key, comp, lab in spec:
        r = sp.get(key)
        row = ((r or {}).get("cluster") or {}).get(comp)
        if not row:
            continue
        pt = (r.get("prespecified_two_way") or {}).get("ci95_R") if comp == "primary" else None
        rmp = (_r(row["R_minus_Rp_Y0"]) + " " + _ci(row["ci95_R_minus_Rp_Y0"])) if "R_minus_Rp_Y0" in row else "---"
        rows.append(f"{lab} & {row.get('decision_Y0', '---')} & {_ci(row['ci95_R_ratio_bootstrap'])} & "
                    f"{_ci(pt) if pt else '---'} & {sci(row['smallest_yardstick_for_equivalence'], 2)} / "
                    f"{sci(r['cluster']['yardsticks']['Y0'], 2)} & {rmp} \\\\")
    return _table("p{4.2cm}lcccc", "Split, comparison & Label (Y0) & Ratio bootstrap, 95\\% & Pre-specified two-way, "
                  "95\\% & Smallest yardstick / Y0 & Real $-$ permuted, 95\\%", rows,
                  "Registered quantities not shown elsewhere (post hoc, amendment B-2). Label: the pre-specified "
                  "decision table applied with the pre-specified yardstick (Y0), with the placebo condition on the "
                  "difference between real and permuted names; it reproduces the confirmatory labels. Ratio bootstrap: "
                  "95\\% interval of $R$ from resampling the ratio of the two increments (cluster bootstrap over "
                  "location groups, countries for the country hold-out). Pre-specified two-way: the protocol's "
                  "robustness check, resampling activations and the country of each cell's nearest settlement "
                  "(1,000 replicates, ratio). Smallest yardstick: the smallest terrain increment against which the "
                  "90\\% interval of $\\Delta\\mathrm{BS}$ would lie within $\\pm0.05$ of it, beside Y0. "
                  "$^{u}$The pre-specified rule leaves $R$ undefined on this split. $^{e}$All areas of each flood "
                  "event on the same side of every split. $^{d}$Name features centred within each AOI "
                  "(M1s$\\to$M3d). ---: the protocol's two-way bootstrap is defined for the primary comparison "
                  "of each split only.",
                  "tab:si-registered")


def si_ratio(v: dict) -> str:
    """補足表: 比の bootstrap による R の 95% 区間を、v36 の全ての比較 × 分割について（B-2 追記 1・3、登録）。"""
    sp = v["splits"]
    comps = [("primary", "Names within areas, M1c$\\to$M3n (H2)"), ("placebo", "Names permuted, M1c$\\to$M3np"),
             ("channel_lexicon", "Lexicon only, M1c$\\to$M3nx"), ("channel_learned", "Learned score only, M1c$\\to$M3nl"),
             ("H1", "Position only, M0c$\\to$M0n (H1)"), ("H1_flexible", "Position only, flexible, M0$\\to$M0N"),
             ("total", "Total, M1s$\\to$M3s"), ("flexible", "Flexible, M1$\\to$M3"),
             ("H4", "Language model, M1rc$\\to$M3rn (H4)"), ("H4_placebo", "Language model, permuted"),
             ("H4_total", "Language model, total, M1s$\\to$M3rs"), ("H4_flexible", "Language model, flexible"),
             ("demeaned", "Names centred within AOIs, M1s$\\to$M3d")]
    splits = [("primary", "Location groups"), ("country", "Countries"), ("temporal", "Time$^{u}$"),
              ("llmcut", "LM cut-off"), ("primary_actgrp", "Events grouped")]
    rows = []
    for key, lab in comps:
        cells = []
        for s, _ in splits:
            src = sp.get(f"{s}_dm" if key == "demeaned" and s != "primary_actgrp" else
                         ("primary_dm_actgrp" if key == "demeaned" else s)) or {}
            row = (src.get("cluster") or {}).get(key) or {}
            cells.append(_ci(row.get("ci95_R_ratio_bootstrap")) if row else "---")
        if any(c != "---" for c in cells):
            rows.append(f"{lab} & " + " & ".join(cells) + " \\\\")
    return _table("p{5.2cm}ccccc", "Comparison & " + " & ".join(n for _, n in splits), rows,
                  "Ratio-bootstrap 95\\% intervals of $R$ for every comparison and split (post hoc, registered in "
                  "amendment B-2): the two increments are resampled together (cluster bootstrap over location "
                  "groups, countries for the country hold-out, 4,000 replicates) and their ratio taken in each "
                  "replicate. The decision table uses fixed-margin intervals instead (Methods), because the ratio "
                  "mixes the uncertainty of the denominator into $R$. $^{u}$The pre-specified rule leaves $R$ "
                  "undefined on this split. ---: comparison not run on this split.",
                  "tab:si-ratio")


def si_posthoc(v36: dict | None, v37: dict | None, v38: dict | None) -> str:
    """補足表: 主比較の事後の層・部分集合・分け方（v36–v38）。登録外の行に印。"""
    rows = []

    def row(lab: str, b: dict | None, y0: str = "---") -> None:
        if not b or b.get("R") is None:
            rows.append(f"{lab} & \\missing{{pending}} & & & & & \\\\")
            return
        eq = b.get("equivalent")
        if eq is None and b.get("ci90_R"):
            eq = -SESOI < b["ci90_R"][0] and b["ci90_R"][1] < SESOI
        rows.append(f"{lab} & {_n(b['n_cells'])} & {b['n_units']} & {y0} & {_r(b['R'])} & "
                    f"{_ci(b.get('ci90_R'))} & {'equivalent' if eq else 'not equivalent'} \\\\")

    if v37:
        rows.append("\\multicolumn{7}{l}{\\emph{Distance from the cell to its nearest settlement}} \\\\")
        for k, lab in (("le500m", "\\quad up to 500\\,m"), ("le250m", "\\quad up to 250\\,m$^{a}$"),
                       ("250_500m", "\\quad 250--500\\,m$^{a}$"),
                       ("500m_1km", "\\quad 0.5--1\\,km"), ("gt1km", "\\quad more than 1\\,km")):
            b = v37["primary_by_distance"].get(k)
            row(lab, b, sci(b["yardstick_M0_M1"], 2) if b else "---")
        b = v37["primary_by_distance_wetside"].get("le500m")
        row("\\quad up to 500\\,m, nearest name water-side$^{a}$", b, sci(b["yardstick_M0_M1"], 2) if b else "---")
        nc = v37.get("near_settlement_clean")
        if nc:
            b = nc["primary"]
            row("\\quad up to 500\\,m, built-up and tree cover removed$^{a}$", b, sci(b["yardstick_M0_M1"], 2))
            b = nc["primary_wetside"]
            row("\\quad up to 500\\,m, water-side, built-up and tree cover removed$^{a}$", b,
                sci(b["yardstick_M0_M1"], 2))
        rows.append("\\multicolumn{7}{l}{\\emph{Class of the nearest name}} \\\\")
        for c, lab in CLASS_LABELS.items():
            b = v37.get("primary_by_class", {}).get(c)
            row(f"\\quad {lab}", b, sci(b["yardstick_M0_M1"], 2) if b else "---")
    if v38:
        rows.append("\\multicolumn{7}{l}{\\emph{Post-event imagery of the delineation}} \\\\")
        for k, lab in (("sar", "\\quad radar only"), ("optical", "\\quad optical only"),
                       ("unknown", "\\quad not recorded (mostly schema generations 0--1)")):
            b = v38["increments"].get(k)
            row(lab, b, sci(b["delta_bs_M0_M1"], 2) if b else "---")
    if v36:
        rows.append("\\multicolumn{7}{l}{\\emph{Fold assignment, grouping and seeds (verdict of the decision table)}} \\\\")
        for run, lab in (("world_v2", "\\quad pre-specified folds"),
                         ("world_v2_fs1", "\\quad folds reassigned (seed 1)"),
                         ("world_v2_fs2", "\\quad folds reassigned (seed 2)"),
                         ("world_v2_actgrp", "\\quad AOIs of the same activation grouped"),
                         ("world_v2_s5_k1", "\\quad one seed instead of five")):
            r = v36.get("spread_primary", {}).get(run)
            met = _json(DATA / "runs" / run / "metrics_primary.json")
            if r and met:
                tost = "TOST holds" if r.get("equivalent_TOST") else "TOST fails"
                rows.append(f"{lab} & {_n(met['n_cells'])} & {met['n_units']} & {sci(r['delta_bs_M0_M1'], 2)} & "
                            f"{_r(r['R'])} & {_ci(r['ci90_R'])} & {r['decision']} ({tost}) \\\\")
            else:
                rows.append(f"{lab} & \\missing{{pending}} & & & & & \\\\")
    return _table("p{5.4cm}rrcrcl", "Subset or variant & Cells & Units & Terrain $\\Delta$BS & $R$ & "
                  "90\\% interval & Verdict", rows,
                  "The primary comparison in subsets of cells and under variants of the fold assignment (post "
                  "hoc). For subsets, $R$ is the names increment in the subset divided by the terrain increment "
                  "in the same subset (as for the H5 strata); the 90\\% interval comes from the cluster bootstrap "
                  "over location groups restricted to the subset, and the last column gives the equivalence test "
                  "alone. For variants, the last column gives the label of the pre-specified decision table and "
                  "the equivalence test; the run with one seed instead of five is the pre-specified sensitivity "
                  "analysis 5 (protocol \\S7-5). Built-up and tree cover removed: cells with ESA CCI built-up cover above "
                  "10\\% or tree cover above 50\\%. Distance bands up to 500\\,m, 0.5--1\\,km and beyond 1\\,km, "
                  "the classes, the sensors, the fold reassignments and the activation grouping were registered in "
                  "amendment B-2 before they were run. The terrain increment within a class differs between classes "
                  "by up to about seventeen-fold, so the verdicts by class depend strongly on their own denominators; "
                  "Supplementary Table~\\ref{tab:si-class} compares classes on the Brier scale. "
                  "$^{a}$Not registered (added while the analysis was written, or a subset of such an analysis).",
                  "tab:si-posthoc")


def si_class(primary: dict | None, v: dict) -> str:
    """補足表: 最寄り名の区分別の主比較の増分を Brier の尺度で（B-2 で登録した記述の量。95% 区間と
    real − permuted は追記 4 で登録）。区分ごとの物差しが大きく違うので、R ではなく ΔBS で並べる。"""

    def e4(x: float) -> str:  # 10^-4 単位、有効 3 桁
        s = f"{x * 1e4:+#.3g}".rstrip(".")
        return sci(x, 2) if "e" in s else "$" + s.replace("-", "{-}").replace("+", "{+}") + "$"

    def ci(c: list | None) -> str:
        return "---" if not c else f"[{e4(c[0])}, {e4(c[1])}]"

    rows = []
    p = (primary or {}).get("primary")
    if p:
        rows.append(f"all cells (primary) & {_n(primary['n_cells'])} & {primary['n_units']} & "
                    f"{e4(p['delta_bs_M1c_M3n'])} & {ci(p['ci90_delta'])} & {ci(p['ci95_delta'])} & "
                    f"---$^{{a}}$ & {p['delta_bs_M0_M1'] * 1e4:.3g} \\\\")
    for c, lab in CLASS_LABELS.items():
        b = v.get("primary_by_class", {}).get(c)
        if not b or "ci95_delta" not in b:
            rows.append(f"{lab} & \\missing{{v37 class {c}}} & & & & & & \\\\")
            continue
        rows.append(f"{lab} & {_n(int(b['n_cells']))} & {b['n_units']} & {e4(b['delta_bs_names'])} & "
                    f"{ci(b['ci90_delta'])} & {ci(b['ci95_delta'])} & "
                    f"{e4(b['real_minus_permuted_delta_bs'])} {ci(b['ci95_real_minus_permuted'])} & "
                    f"{b['yardstick_M0_M1'] * 1e4:.3g} \\\\")
    margin = f"{p['margin_delta_bs'] * 1e4:.2f}" if p else "\\missing{margin}"
    return _table("lrrrccll", "Nearest name & Cells & Units & $\\Delta$BS & 90\\% interval & 95\\% interval & "
                  "Real${}-{}$permuted (95\\%) & Terrain $\\Delta$BS", rows,
                  "The names increment within areas (M1c$\\to$M3n) by class of the nearest name, on the Brier scale "
                  "(post hoc; registered in amendment B-2 as descriptive quantities, the 95\\% intervals and the "
                  "difference between real and permuted names in its fourth addendum). All Brier-score differences "
                  "are in units of $10^{-4}$; positive values mean that names improved held-out predictions. "
                  "Intervals come from the cluster bootstrap over location groups restricted to the cells of each "
                  "class (classes overlap: a name can belong to several). Terrain $\\Delta$BS is the terrain "
                  "increment (M0$\\to$M1) in the same cells, the denominator of the class verdicts in "
                  f"Supplementary Table~\\ref{{tab:si-posthoc}}; the pooled equivalence margin is {margin}. "
                  "$^{a}$For all cells, the difference between real and permuted names is given as a ratio in "
                  "Supplementary Table~\\ref{tab:si-registered}.",
                  "tab:si-class")


def si_positive_control(v: dict) -> str:
    """補足表: 正の対照（ラベル非依存、v37）。全セルの水系側の群と高地は B-2 で登録、区分の分割と集落の点は登録外。"""
    pc = v["positive_control"]

    def num(v: float, nd: int) -> str:
        while nd < 3 and v != 0 and round(abs(v), nd) == 0:  # ±0.0 と書かない
            nd += 1
        return f"{v:+.{nd}f}"

    def f(r: dict | None, scale: float = 100.0, nd: int = 1) -> str:
        if not r:
            return "---"
        s = (f"{num(scale * r['estimate'], nd)} [{num(scale * r['ci95'][0], nd)}, "
             f"{num(scale * r['ci95'][1], nd)}]")
        return s.replace("+", "{+}").replace("-", "{-}")

    def base(r: dict | None) -> str:
        return "---" if not r else f"{100 * r['mean_group']:.0f}/{100 * r['mean_other']:.0f}"

    rows_a, rows_b = [], []
    for c, lab in CLASS_LABELS.items():
        for scale, slab in (("cells", "all cells"), ("settlements", "at settlement$^{a}$")):
            r = pc[scale].get(c) or {}
            n = (r.get("low_hand") or {}).get("n_aoi", "---")
            name = lab if scale == "cells" else ""
            rows_a.append(f"{name} & {slab} & {f(r.get('low_hand'))} & {base(r.get('low_hand'))} & "
                          f"{f(r.get('low_half'))} & {n} \\\\")
            rows_b.append(f"{name} & {slab} & {f(r.get('log_river'), 1.0, 2)} & {f(r.get('hand_mean'), 1.0, 1)} & "
                          f"{f(r.get('hand_pct'))} \\\\")
    for c in ("wetside", "high"):
        r = pc.get("settlements_not_rounded", {}).get(c) or {}
        slab = "at settlement, coordinates not rounded$^{b}$"
        rows_a.append(f"{CLASS_LABELS[c]} & {slab} & {f(r.get('low_hand'))} & {base(r.get('low_hand'))} & "
                      f"{f(r.get('low_half'))} & {(r.get('low_hand') or {}).get('n_aoi', '---')} \\\\")
        rows_b.append(f"{CLASS_LABELS[c]} & {slab} & {f(r.get('log_river'), 1.0, 2)} & "
                      f"{f(r.get('hand_mean'), 1.0, 1)} & {f(r.get('hand_pct'))} \\\\")
    notes = ("The all-cell rows for water-side names as a group and for highland names are the analysis registered "
             "in amendment B-2; the other all-cell rows (the split by class) were computed with it and are not "
             "registered. $^{a}$Not registered: cells within "
             "177\\,m of the settlement point (the one to four cells around it), averaged per settlement "
             f"({_n(pc['n_settlement_points'])} settlements; AOIs need at least two settlements of each kind). "
             "$^{b}$Excluding settlements whose GeoNames coordinates fall on whole arc-minutes (rounding of up "
             f"to about 0.9\\,km in each coordinate; {pct(pc['share_settlement_points_minute_rounded'])} of "
             "settlement points).")
    a = _table("llcccr", "Nearest name & Scale & HAND $\\leq2$\\,m in cell (pp) & "
               "Class/other (\\%) & Lower half (pp) & AOIs", rows_a,
               "Positive control without flood labels (binary measures): within-AOI differences between cells (or "
               "settlements) whose nearest settlement name is in the class and other cells (or settlements) of the "
               "same AOI, with 95\\% cluster-bootstrap intervals over location groups. ``HAND $\\leq2$\\,m'' is "
               "whether any 25\\,m pixel of the 250\\,m cell lies within 2\\,m of the nearest drainage (the "
               "cell's minimum HAND); it holds for most cells, so the shares in each group are given. Lower half: the "
               "cell lies in the lower half of its AOI by mean elevation. Continuous measures are in Supplementary "
               "Table~\\ref{tab:si-positive-control-b}. " + notes,
               "tab:si-positive-control")
    b = _table("llccc", "Nearest name & Scale & log(1 + river distance) & Mean HAND (m) & HAND percentile (pp)", rows_b,
               "Positive control without flood labels (continuous measures), as in Supplementary "
               "Table~\\ref{tab:si-positive-control}: differences in log(1 + distance to the nearest mapped river), "
               "in the cell's mean height above nearest drainage and in its HAND percentile within the AOI. " + notes,
               "tab:si-positive-control-b")
    return a + "\n" + b


def si_band_contrast(v: dict) -> str:
    """補足表: 距離帯ごとの浸水の対比（v37）。対比そのものは登録外、重み 2 通りと AOI・場所の数での報告は B-2 追記 2。"""
    fc = v["flood_contrast"]["by_band"]
    labels = {"le250m": "$\\leq250$\\,m", "250_500m": "250--500\\,m", "500m_1km": "0.5--1\\,km",
              "gt1km": "$>1$\\,km"}

    def c(r: dict | None) -> str:
        return "---" if not r else f"{_pp(r['estimate'])} [{_pp(r['ci95'][0])}, {_pp(r['ci95'][1])}]"

    rows = []
    for b, lab in labels.items():
        cap, unc = fc["capped"].get(b) or {}, fc["uncapped"].get(b) or {}
        n = cap.get("y") or {}
        rows.append(f"{lab} & {n.get('n_aoi', '---')} & {n.get('n_units', '---')} & {c(cap.get('y'))} & "
                    f"{c(cap.get('resid_M1c'))} & {c(unc.get('y'))} & {c(unc.get('resid_M1c'))} \\\\")
    head = ("Distance to the nearest settlement & AOIs & Places & \\multicolumn{2}{c}{Weights: smaller group} & "
            "\\multicolumn{2}{c}{Weights: uncapped AOI size} \\\\\n"
            " & & & observed & net of terrain & observed & net of terrain")
    return _table("lrrcccc", head, rows,
                  "Within-AOI flood contrasts by distance from the cell to its nearest settlement (not registered; "
                  "added while the registered analyses were written, with their full reporting registered in "
                  "amendment B-2, addendum 2): the share of cells labelled flooded nearest a water-side name minus "
                  "the share nearest other names, in percentage points, observed and net of the held-out "
                  "predictions of the terrain model with the regional composition of names (M1c), with 95\\% "
                  "cluster-bootstrap intervals over location groups. AOIs enter a band when both groups have at "
                  "least 20 cells in it; they are averaged with weights equal to the size of the smaller group, or "
                  "scaled to the uncapped number of evaluable cells of each AOI, which lets a few large AOIs "
                  "dominate.",
                  "tab:si-band-contrast")


def si_odds(v: dict) -> str:
    """補足表: AOI 固定効果つきロジスティック回帰のオッズ比（v39）。"""
    def cell(t: dict | None) -> str:
        if not t:
            return "---"
        return f"{t['odds_ratio']:.2f} [{t['ci95_odds_ratio'][0]:.2f}, {t['ci95_odds_ratio'][1]:.2f}]"

    rows = []
    for c, lab in CLASS_LABELS.items():
        key = "wetside" if c == "wetside" else "categories"
        cells = [cell(v.get(f"{key}_allcells_{adj}", {}).get("terms", {}).get(f"{c}:all"))
                 for adj in ("raw", "adjusted_M1")]
        for adj in ("raw", "adjusted_M1"):
            for band in ("le500m", "500m_1km", "gt1km"):
                cells.append(cell(v[f"{key}_{adj}"]["terms"][f"{c}:{band}"]))
        rows.append(f"{lab} & " + " & ".join(cells) + " \\\\")
    rows.append("\\midrule")
    band_cells = ["---", "---"]
    for adj in ("raw", "adjusted_M1"):
        band_cells += ["1 (reference)", cell(v[f"wetside_{adj}"]["terms"]["band_500m_1km"]),
                       cell(v[f"wetside_{adj}"]["terms"]["band_gt1km"])]
    rows.append("Distance band (all names) & " + " & ".join(band_cells) + " \\\\")
    if "wetside_clean_adjusted_M1" in v:
        clean = ["---", "---"] + [cell(v[f"wetside_clean_{adj}"]["terms"][f"wetside:{band}"])
                                  for adj in ("raw", "adjusted_M1") for band in ("le500m", "500m_1km", "gt1km")]
        rows.append("water-side, built-up and tree cover removed & " + " & ".join(clean) + " \\\\")
    head = ("Nearest name & \\multicolumn{2}{c}{All cells$^{r}$} & \\multicolumn{3}{c}{Unadjusted, by distance} & "
            "\\multicolumn{3}{c}{Given the terrain model, by distance} \\\\\n"
            " & unadjusted & given terrain & $\\leq500$\\,m & 0.5--1\\,km & $>1$\\,km & $\\leq500$\\,m & 0.5--1\\,km & $>1$\\,km")
    bc = v.get("bootstrap_check_wetside")
    check = ""
    if bc:
        b = bc["adjusted_M1"]["wetside:le500m"]["odds_ratio_ci95_bootstrap"]
        check = (f" A cluster bootstrap over location groups ({bc['reps']} replicates, AOIs re-entered with "
                 "their own intercepts) gives similar intervals, e.g. "
                 f"{b[0]:.2f}--{b[1]:.2f} for water-side names within 500\\,m given the terrain model.")
    return _table("lcccccccc", head, rows,
                  "Within-AOI odds ratios of flooding for cells whose nearest settlement name is in the class, "
                  "from logistic regressions with AOI fixed effects, unadjusted or with the logit of the held-out "
                  "terrain-model prediction (M1) as a covariate; 95\\% intervals from cluster-robust "
                  f"(location-group) standard errors; {v['wetside_adjusted_M1']['n_aoi']} AOIs (AOIs whose sampled "
                  "cells are all dry carry no information and are omitted). $^{r}$All cells, one indicator per "
                  "class: the analysis registered in amendment B-2. The distance-banded models (class $\\times$ "
                  "band, with band main effects) were not registered. The distance-band row gives the odds of a "
                  "flood label 0.5--1\\,km and more than 1\\,km from the nearest settlement relative to within "
                  "500\\,m. Built-up and tree cover removed: cells with ESA CCI built-up cover above 10\\% or tree "
                  "cover above 50\\% excluded." + check,
                  "tab:si-odds", landscape=True)


def si_composition(v: dict) -> str:
    """補足表: 名前の地域構成の増分（再較正だけの基準から、v36）。時間外しは物差しが定義されない。"""
    sp = v["splits"]
    rows = []
    for key, lab in (("primary", "Location groups"), ("primary_actgrp", "Flood events grouped$^{e}$"),
                     ("country", "Countries"), ("temporal", "Time$^{u}$"), ("llmcut", "LM cut-off")):
        r = sp.get(key)
        if not r:
            continue
        comp = r["cluster"]["composition"]
        y0 = r["cluster"]["yardsticks"]["Y0"]
        cells = []
        for ck in ("M0->M0c", "M0s->M0c", "T0s->T0c", "M1s->M1c", "M1s->M1rc"):
            c = comp.get(ck)
            cells.append("---" if not c else f"{_r(c['relative_to_Y0'])} {_ci([x / y0 for x in c['ci95_delta']])}")
        rows.append(f"{lab} & " + " & ".join(cells) + " \\\\")
    r = sp.get("llmcut")
    if r:  # B-2 追記 2 の 3: 言語モデルの読みの地域構成を全ての物差しで（クラスタと二方向）
        rows += ["\\midrule",
                 "\\multicolumn{6}{l}{\\emph{Composition of the language-model readings (M1s$\\to$M1rc, LM "
                 "cut-off), by yardstick}} \\\\",
                 " & Y0 & Y1 & Y2 & Y3 & Y4 \\\\"]
        for kind, kname in (("cluster", "cluster"), ("two_way", "two-way")):
            c = r[kind]["composition"].get("M1s->M1rc")
            yv = r[kind]["yardsticks"]
            cells = []
            for y in ("Y0", "Y1", "Y2", "Y3", "Y4"):
                yy = yv.get(y)
                if not c or yy is None:
                    cells.append("---")
                elif yy <= 0:
                    cells.append("undefined")
                else:
                    cells.append(f"{_r(c['delta_bs'] / yy)} {_ci([x / yy for x in c['ci95_delta']])}")
            rows.append(f"\\quad {kname} & " + " & ".join(cells) + " \\\\")
    return _table("lccccc", "Split & M0$\\to$M0c & M0s$\\to$M0c & T0s$\\to$T0c & M1s$\\to$M1c & "
                  "M1s$\\to$M1rc (LM)", rows,
                  "Increment from the regional composition of names (AOI means of the name features; for the "
                  "last column, of the language-model readings), as a fraction of the pre-specified terrain "
                  "increment, with 95\\% cluster-bootstrap intervals (descriptive). M0$\\to$M0c compares the "
                  "flexible position model with a stacked model and therefore includes the gain from "
                  "recalibration; M0s and T0s are the same base models recalibrated without the composition "
                  "(post hoc, amendment B-2), so M0s$\\to$M0c and T0s$\\to$T0c isolate the composition. "
                  "$^{e}$All areas of each flood event on the same side of every split (no T0 arm in this run). "
                  "$^{u}$The pre-specified rule leaves the terrain increment undefined on this split (its 95\\% "
                  "interval includes zero), so these fractions are only indicative. The lower block gives the "
                  "composition of the language-model readings against each yardstick of Supplementary "
                  "Table~\\ref{tab:si-yardsticks} (amendment B-2, addendum 2): point estimate and 95\\% interval "
                  "under the bootstrap over places and under the two-way bootstrap over activations and countries; "
                  "undefined where the yardstick is not positive.", "tab:si-composition")


def real_margin_rate(runs: Path, k: str, split: str) -> tuple[int, int, int] | None:
    """δ = 1.2 の各回の ΔBS の 90% 区間（R の区間 × その回の合成の物差し）が、実データの余白の内側に入った回の数・回数・
    合成の物差しが 0 以下だった回の数（B-2 追記 4）。余白が未定義の分割は None。"""
    ses = _json(runs / f"power_v3{k}_sesoi" / "power_summary.json")
    met = (_json(runs / "world_v2" / f"metrics_{split}.json") or {}).get("primary")
    if not ses or not met or met["decision"] == "undefined-denominator" or not ses.get("reps"):
        return None
    mg = met["margin_delta_bs"]
    inside = nonpos = 0
    for r in ses["reps"]:
        y0 = r["delta_bs_M0_M1"]
        lo, hi = sorted([r["ci90"][0] * y0, r["ci90"][1] * y0])
        inside += -mg < lo and hi < mg
        nonpos += y0 <= 0
    return inside, len(ses["reps"]), nonpos


def si_power_splits(runs: Path) -> str:
    """補足表: 各分割の検出力（主分割は凍結前、他は事後 B-2）。δ = 0 と δ = 1.2 の判定の分布・Wilson 区間・平均 R、
    余白と観測の ΔBS の 95% 上限を δ = 1.2 の平均 ΔBS と比べた比（B-2 追記 3）。"""
    rows = []
    for k, split, lab in (("", "primary", "Location groups (pre-specified)"), ("c", "country", "Countries"),
                          ("t", "temporal", "Time"), ("l", "llmcut", "Language-model cut-off, lexicon effect$^{a}$")):
        nul = _json(runs / f"power_v3{k}_null" / "power_summary.json")
        ses = _json(runs / f"power_v3{k}_sesoi" / "power_summary.json")
        if not nul or not ses:
            rows.append(f"{lab} & \\missing{{pending}} & & & & & & & & \\\\")
            continue

        def rate(x: float, n: int) -> str:
            lo, hi = _wilson(x, n)
            return f"{100 * x:.0f} [{100 * lo:.0f}, {100 * hi:.0f}]"

        rm = real_margin_rate(runs, k, split)
        real = rate(rm[0] / rm[1], rm[1]) if rm else "---"
        if rm and rm[2]:  # 合成の物差しが負の回（R の区間が裏返る）を数えた注
            real += f"$^{{n}}$"
        sr = split_ratios(runs, k, split)
        ratios = ("---" if not sr or sr["margin"] is None else f"{sr['margin']:.1f}",
                  "---" if not sr else f"{sr['upper']:.2f}")
        for d, s, first in (("0", nul, True), ("1.2", ses, False)):
            labels = ", ".join(f"{kk} {vv}" for kk, vv in sorted(s["decision_counts"].items()))
            tail = " & --- & --- & ---" if first else f" & {real} & {ratios[0]} & {ratios[1]}"
            rows.append(f"{lab if first else ''} & {d} & {s['n_reps']} & {rate(s['rate_equivalent_TOST'], s['n_reps'])} & "
                        f"{rate(s['rate_positive'], s['n_reps'])} & {_r(s['mean_R'])} & {labels}{tail} \\\\")
    return _table("p{3.4cm}cccccp{4.6cm}ccc", "Split & $\\delta$ & Replicates & Equivalent (\\%) & Positive (\\%) & "
                  "Mean $R$ & Labels & Equivalent, real margin (\\%) & Margin ratio & Upper ratio", rows,
                  "Operating characteristics of the decision rule on synthetic labels for each split (the "
                  "primary split before the freeze; the others post hoc, amendment B-2), with 95\\% Wilson "
                  "intervals for the Monte-Carlo error. The simulated effect is $\\delta$ times the indicator that "
                  "the nearest name is water-side in the lexicon, calibrated so that $R\\approx0.05$ on the primary "
                  "split ($\\delta=1.2$, an odds ratio of 3.3); on the other splits the same $\\delta$ gives the mean "
                  "$R$ shown, not the smallest effect of interest. Equivalent counts the replicates whose 90\\% "
                  "interval of $R$ lies within $\\pm0.05$, whatever their label. The last column gives, for "
                  "$\\delta=1.2$, the share of replicates whose 90\\% interval of $\\Delta$BS lies within the "
                  "equivalence margin of the real labels on the same split (post hoc, amendment B-2): how often an "
                  "effect of this size would have been called equivalent in our data (---: margin undefined on the "
                  "time split). Margin ratio and upper ratio (amendment B-2, addendum 3): the real-label "
                  "equivalence margin and the upper end of the observed 95\\% interval of $\\Delta$BS on the same "
                  "split, each divided by the mean $\\Delta$BS that $\\delta=1.2$ produced there; these ratios "
                  "ignore that intervals widen under the alternative, so the column before them is the better guide. "
                  "$^{n}$Includes replicates whose simulated terrain increment was negative (the interval of $R$ is "
                  "then reversed; Supplementary Note~10). $^{a}$The effect is injected through the lexicon, so this "
                  "row is the lexicon test on the language-model split; there is no power analysis for the "
                  "language-model reading (H4).",
                  "tab:si-power-splits", landscape=True)


def _wilson(p: float, n: int) -> list[float]:
    z = 1.96
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return [max(0.0, c - h), min(1.0, c + h)]


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
    if 0 < abs(x) < 5e-5:  # 4 桁でも 0.0000 に丸まる値は指数表記で
        e = int(np.floor(np.log10(abs(x))))
        return f"${x / 10**e:+.0f}\\times10^{{{e}}}$".replace("+", "{+}").replace("-", "{-}")
    s = f"{x:+.4f}" if 0 < abs(x) < 0.0005 else f"{x:+.3f}"  # 非零が 0.000 に丸まらないように
    return f"${s}$".replace("+", "{+}").replace("-", "{-}")


def _ci(c: list | None) -> str:
    return "---" if not c else f"[{_r(c[0])}, {_r(c[1])}]"


def _n(x: int) -> str:
    return f"{x:,}".replace(",", "{,}")


def _table(cols: str, head: str, rows: list[str], caption: str, label: str, *, landscape: bool = False) -> str:
    """表の LaTeX。landscape=True は横置き（rotating の sidewaystable、列の多い表を縮めすぎないように）。"""
    body = "\n".join(rows)
    env = "sidewaystable" if landscape else "table"
    pos = "" if landscape else "[h!]"
    return (f"\\begin{{{env}}}{pos}\n\\centering\\footnotesize\\setlength{{\\tabcolsep}}{{4pt}}\n"
            f"\\caption{{{caption}}}\\label{{{label}}}\n"
            f"\\begin{{adjustbox}}{{max width=\\linewidth}}\n"
            f"\\begin{{tabular}}{{{cols}}}\n\\toprule\n{head} \\\\\n\\midrule\n{body}\n"
            f"\\bottomrule\n\\end{{tabular}}\n\\end{{adjustbox}}\n\\end{{{env}}}\n")


def _tost(blk: dict) -> str:
    """§6: TOST の成否を判定 label と独立に併記する。"""
    if "equivalent_TOST" not in blk or blk.get("decision") == "undefined-denominator":
        return "---"
    return "yes" if blk["equivalent_TOST"] else "no"


def _cmp_row(name: str, blk: dict | None, r: str = "R", c95: str = "ci95_R",
             c90: str | None = "ci90_R") -> str:
    if not blk:
        return f"{name} & \\missing{{{name}}} & & & & \\\\"
    v90 = _ci(blk.get(c90)) if c90 else "---"
    return (f"{name} & {_r(blk.get(r))} & {_ci(blk.get(c95))} & {v90} & "
            f"{blk.get('decision', '---')} & {_tost(blk)} \\\\")


def _row7(name: str, r: float | None, c95: list | None, c90: list | None, ratio: list | None,
          verdict: str = "---", tost: str = "---") -> str:
    return (f"{name} & {_r(r)} & {_ci(c95)} & {_ci(c90)} & {_ci(ratio)} & {verdict} & {tost} \\\\")


def si_comparisons(m: dict, v36: dict | None = None) -> str:
    """補足表 1: 主分割の全比較（副次・チャネル・梯子を含む）。区間は固定マージン（物差し = 登録どおりの
    点推定）と比の bootstrap（物差しの揺らぎも含む、記述）の両方。チャネル・H1 の固定マージンは v36。"""
    d01 = m["primary"]["delta_bs_M0_M1"]
    chf = m.get("channels_flexible", {})
    pl, tpl, fpl = m.get("placebo", {}), m.get("total_placebo", {}), m.get("flexible_placebo", {})

    def verdict_row(name: str, blk: dict | None) -> str:
        if not blk:
            return f"{name} & \\missing{{{name}}} & & & & & \\\\"
        return _row7(name, blk["R"], blk["ci95_R"], blk["ci90_R"], blk.get("ci95_R_ratio_bootstrap"),
                     blk["decision"], _tost(blk))

    def fm_row(name: str, comp: str) -> str:
        f = fixed_margin(v36, "primary", comp)
        if not f:
            return f"{name} & \\missing{{v36 {comp}}} & & & & & \\\\"
        return _row7(name, f["R"], f["ci95_R"], f["ci90_R"], f["ratio95"])

    def plac_row(name: str, blk: dict, comp: str | None = None) -> str:
        f = fixed_margin(v36, "primary", comp) if comp else None
        return _row7(name, blk.get("R_placebo"), blk.get("ci95_R_placebo"), f["ci90_R"] if f else None,
                     f["ratio95"] if f else None)

    rows = [
        verdict_row("H2: names within areas, M1c$\\to$M3n", m["primary"]),
        plac_row("\\quad placebo, M1c$\\to$M3np", pl, "placebo"),
        _row7("\\quad real minus placebo", pl.get("R_minus_Rp"), pl.get("ci95_R_minus_Rp"), None, None),
        fm_row("\\quad lexicon only, M1c$\\to$M3nx", "channel_lexicon"),
        fm_row("\\quad learned score only, M1c$\\to$M3nl", "channel_learned"),
        fm_row("H1: position only, M0c$\\to$M0n", "H1"),
        fm_row("\\quad flexible, M0$\\to$M0N", "H1_flexible"),
        verdict_row("H2t: total increment, M1s$\\to$M3s", m.get("total")),
        plac_row("\\quad placebo, M1s$\\to$M3ps", tpl),
        verdict_row("Flexible refit, M1$\\to$M3", m.get("flexible")),
        plac_row("\\quad placebo, M1$\\to$M3p", fpl),
        _row7("\\quad lexicon only, M1$\\to$M3x", (chf.get("M3x") or {}).get("R"), None, None,
              (chf.get("M3x") or {}).get("ci95_R")),
        _row7("\\quad learned score only, M1$\\to$M3l", (chf.get("M3l") or {}).get("R"), None, None,
              (chf.get("M3l") or {}).get("ci95_R")),
    ]
    for key, label in (("ladder", "H3 ladder, within areas"), ("ladder_flexible",
                                                                "H3 ladder, flexible")):
        lad = m.get(key)
        if not lad:
            continue
        for k, (rung, dlt, c) in enumerate(zip(lad["rungs"], lad["delta_names"],
                                               lad["ci95_delta_names"], strict=True)):
            name = f"{label}, T{k}" if k < 4 else f"{label}, T4 = M1"
            arrow = rung.replace("->", "$" + TO + "$")
            rows.append(_row7(f"{name} ({arrow})", dlt / d01, [x / d01 for x in c], None, None))
    return _table("p{5.4cm}rcccll", "Comparison & $R$ & 95\\% (fixed margin) & 90\\% (fixed margin) & "
                  "95\\% (ratio bootstrap) & Verdict & TOST",
                  rows, "All pre-specified comparisons on the primary split (location groups held "
                  "out). $R$ is the change in held-out Brier score as a fraction of the terrain "
                  "increment $\\Delta\\mathrm{BS}(\\mathrm{M0}\\to\\mathrm{M1})$. Fixed-margin intervals "
                  "divide the cluster-bootstrap interval of the change in Brier score by the point estimate "
                  "of the terrain increment (the form of the decision rule); ratio-bootstrap intervals "
                  "resample the ratio and so also carry the uncertainty of the terrain increment "
                  "(descriptive). Location groups are resampled 4,000 times. Verdicts are given for the "
                  "comparisons to which the protocol attaches a decision rule; the same comparison under "
                  "other yardsticks and a two-way bootstrap is in Supplementary "
                  "Table~\\ref{tab:si-yardsticks}.",
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


def si_splits_secondary(primary: dict, extra: dict[str, dict], v36: dict | None = None) -> str:
    """補足表 3: 分割ごとの副次比較（全体の増分・柔軟モデル、llmcut は N-read 版も）と placebo との差。"""
    spec = [("Location groups", primary, "primary"), ("Countries", extra.get("country"), "country"),
            ("Time", extra.get("temporal"), "temporal"), ("Language-model cut-off", extra.get("llmcut"), "llmcut")]
    comps = [("total", "total_placebo", "total, M1s$\\to$M3s"),
             ("flexible", "flexible_placebo", "flexible, M1$\\to$M3"),
             ("H4_total", "H4_total_placebo", "language model, total, M1s$\\to$M3rs"),
             ("H4_flexible", "H4_flexible_placebo", "language model, flexible, M1$\\to$M3r")]
    rows = []
    for split, met, skey in spec:
        if not met:
            rows.append(f"{split} & \\missing{{pending}} & & & & & \\\\")
            continue
        for key, pkey, lab in comps:
            blk, pl = met.get(key), met.get(pkey) or {}
            if not blk:
                continue
            ratio = (((v36 or {}).get("splits", {}).get(skey, {}).get("cluster", {}).get(key) or {})
                     .get("ci95_R_ratio_bootstrap"))
            rows.append(f"{split}: {lab} & {_r(blk['R'])} & {_ci(blk['ci95_R'])} & {_ci(ratio)} & "
                        f"{_ci(pl.get('ci95_R_minus_Rp'))} & {blk['decision']} & {_tost(blk)} \\\\")
    return _table("p{5.6cm}rccclc", "Split: comparison & $R$ & 95\\% interval & Ratio bootstrap, 95\\% & "
                  "real $-$ placebo (95\\%) & Verdict & TOST", rows,
                  "Secondary comparisons under each split: the total increment of names over the "
                  "recalibrated terrain model (including the regional composition of names) and the "
                  "flexible model refitted with the name features, for the lexicon and learned score "
                  "and, on the language-model split, for the language-model reading. The 95\\% interval is the "
                  "fixed-margin interval of the decision table (Methods); the fourth "
                  "column gives the interval from resampling the ratio of the two increments, and the fifth the "
                  "interval of $R$ minus the same comparison with names permuted among the settlements of each "
                  "AOI.", "tab:si-splits-secondary")


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


def si_strata(m: dict, v36: dict | None = None) -> str:
    """補足表: H5 の層（記述、登録）と、建物群・大陸・言語の層（v36、事後、90% 区間つき）。"""
    h5 = m.get("H5_strata", {})
    names = {"nearest_name_in_lexicon": "Nearest name in the lexicon",
             "not_nearest_name_in_lexicon": "Nearest name not in the lexicon",
             "nearest_name_wetside": "Nearest name water-side",
             "not_nearest_name_wetside": "Nearest name not water-side",
             "low_relief": "Slope at or below the median",
             "not_low_relief": "Slope above the median"}
    rows = [f"{lab} & {_n(h5[k]['n'])} & {100 * h5[k]['wet_rate']:.1f} & {_r(h5[k]['R_within'])} & --- \\\\"
            for k, lab in names.items() if k in h5]
    if v36:
        st = v36["splits"]["primary"]["strata"]
        rows.append("\\midrule")
        for k, lab in (("built_up_gt_10pct", "Built-up cover above 10\\%"),
                       ("not_built_up_gt_10pct", "Built-up cover at most 10\\%")):
            if k in st:
                s = st[k]
                rows.append(f"{lab} & {_n(s['n_cells'])} & {100 * s['wet_rate']:.1f} & {_r(s['R'])} & "
                            f"{_ci(s['ci90_R'])} \\\\")
        rows.append("\\midrule")
        for k in sorted(k for k in st if k.startswith("continent_")):
            s = st[k]
            rows.append(f"{CONTINENTS.get(k[10:], k[10:])} & {_n(s['n_cells'])} & {100 * s['wet_rate']:.1f} & "
                        f"{_r(s['R'])} & {_ci(s['ci90_R'])} \\\\")
        rows.append("\\midrule")
        langs = sorted((k for k in st if k.startswith("language_")), key=lambda k: -st[k]["n_cells"])
        for k in langs:
            s = st[k]
            if s["n_cells"] < 20000:
                continue
            rows.append(f"Language: {k[9:]} & {_n(s['n_cells'])} & {100 * s['wet_rate']:.1f} & {_r(s['R'])} & "
                        f"{_ci(s['ci90_R'])} \\\\")
    return _table("lrrrc", "Stratum & Cells & Flooded (\\%) & $R$ within stratum & 90\\% interval", rows,
                  "Descriptive strata for the primary comparison: the H5 strata of the protocol (first block, no "
                  "intervals as pre-specified; built-up cover and continent are also strata of the protocol and are "
                  "shown with post hoc intervals) and the post hoc strata by the first language of the country of "
                  "the nearest settlement (languages with at least 20,000 cells), with 90\\% intervals from the "
                  "cluster bootstrap over location groups within the stratum. ---: the terrain increment within the "
                  "stratum is not positive, so $R$ is undefined. $R$ within a "
                  "stratum is the names increment in the stratum divided by the terrain increment in the same "
                  "stratum; small strata have small and unstable terrain increments, hence wide intervals.",
                  "tab:si-strata")


def si_metrics(m: dict) -> str:
    """補足表 5: 副次評価量。"""
    sm = m.get("secondary_metrics", {})
    names = {"M1": "M1 (terrain model)", "M3": "M3 (refitted with names)",
             "M1s": "M1s (recalibrated M1)", "M3s": "M3s (+ names)",
             "M1c": "M1c (+ regional composition)", "M3n": "M3n (+ names within areas)"}
    rows = [f"{lab} & {sm[k]['brier']:.5f} & {sm[k]['log_loss']:.5f} & {sm[k]['pr_auc']:.3f} & "
            f"{100 * sm[k]['capture10']:.1f} \\\\" for k, lab in names.items() if k in sm]
    return _table("lrrrr", "Model & Brier & Log loss & PR-AUC & Top-10\\% capture (\\%)", rows,
                  "Secondary metrics on the primary split (pooled over held-out folds).",
                  "tab:si-metrics")


def sens_macros(runs: Path) -> dict[str, str]:
    """§7 の感度の判定の要約（本文用）。全て metrics_primary.json から。"""
    got = {tag: (m or {}).get("primary") for tag, _ in SENS
           if (m := _json(runs / f"world_v2_{tag}" / "metrics_primary.json"))}
    out = {"NSens": str(len(SENS)), "NSensDone": str(len(got)), "NSensWord": _word(len(SENS))}
    if len(got) < len(SENS):
        MISSING.append(f"sensitivity runs {len(got)}/{len(SENS)}")
    out["NSensEquiv"] = str(sum(1 for b in got.values() if b and b["decision"] == "equivalent-null"))
    out["NSensEquivWord"] = _word(out["NSensEquiv"])
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
        nb = (_json(runs / "world_v2_s2a_nocap" / "metrics_primary.json") or {}).get("brier", {})
        if "M1c" in nb and "M1" in nb:  # 地域構成で積み上げた基準が、自分の柔軟モデルより悪い（B-2 追記 4、記述）
            out |= {"BrierNoCapMOnec": f"{nb['M1c']:.5f}", "BrierNoCapMOne": f"{nb['M1']:.5f}",
                    "NoCapLossShare": pct((nb["M1c"] - nb["M1"]) / nocap["delta_bs_M0_M1"])}
        cells: Counter = Counter()
        for f in sorted((runs / "world_v2_s2a_nocap" / "predictions" / "primary").glob("fold*.npz")):
            u, c = np.unique(np.load(f)["unit"], return_counts=True)
            cells.update(dict(zip(u.tolist(), c.tolist(), strict=True)))
        if cells:
            n = np.array(list(cells.values()), dtype=float)
            out |= {"KishNoCap": f"{n.sum() ** 2 / (n ** 2).sum():.0f}",
                    "LargestShareNoCap": pct(n.max() / n.sum())}
    return out


def write_si_tables(primary: dict | None, extra: dict[str, dict], post: dict[str, dict | None]) -> None:
    """補足表を表ごとの別ファイルに（si.tex が対応する Note の直後で読む）。post = 事後の v36–v39。"""
    v36, v37, v38, v39 = (post.get(k) for k in ("v36", "v37", "v38", "v39"))
    tables = {"comparisons": lambda: si_comparisons(primary, v36),
              "splits": lambda: si_splits(primary, extra),
              "splits_secondary": lambda: si_splits_secondary(primary, extra, v36),
              "sensitivity": lambda: si_sensitivity(DATA / "runs"),
              "strata": lambda: si_strata(primary, v36), "metrics": lambda: si_metrics(primary),
              "yardsticks": lambda: si_yardsticks(v36) if v36 else "\\missing{v36}",
              "posthoc": lambda: si_posthoc(v36, v37, v38),
              "class": lambda: si_class(primary, v37) if v37 else "\\missing{v37}",
              "positive_control": lambda: si_positive_control(v37) if v37 else "\\missing{v37}",
              "odds": lambda: si_odds(v39) if v39 else "\\missing{v39}",
              "band_contrast": lambda: si_band_contrast(v37) if v37 else "\\missing{v37}",
              "composition": lambda: si_composition(v36) if v36 else "\\missing{v36}",
              "power_splits": lambda: si_power_splits(DATA / "runs"),
              "registered": lambda: si_registered(v36) if v36 else "\\missing{v36}",
              "ratio": lambda: si_ratio(v36) if v36 else "\\missing{v36}"}
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
        ver = REPO / "docs" / "verification"
        post = {k: _json(ver / f) for k, f in (("v36", "v36_yardstick_robustness.json"),
                                               ("v37", "v37_positive_control.json"),
                                               ("v38", "v38_ems_sensors.json"),
                                               ("v39", "v39_within_aoi_odds.json"))}
        m |= detail_macros(primary, post["v36"])
        m |= split_macros(extra)
        m["NSecondaryVerdicts"] = str(int(m["NSecondaryVerdicts"])
                                      + sum(1 for k in ("total", "flexible") if primary.get(k)))
        fig_main(_forest_rows(primary, extra, post["v36"]), primary["primary"]["decision"])
        if primary.get("ladder"):
            fig_ladder(primary["ladder"], primary["primary"]["delta_bs_M0_M1"])
        write_si_tables(primary, extra, post)
        m |= robustness_macros(post["v36"])
        m |= positive_control_macros(post["v37"])
        m |= sensor_macros(post["v38"])
        m |= odds_macros(post["v39"])
        m |= power_split_macros(DATA / "runs")
        m |= hgb_macros(_json(ver / "v40_hgb_iterations.json"))
        if post["v37"]:
            fig_names(post["v37"], post["v39"])
        m |= sens_macros(DATA / "runs")
        diag = run / "diagnostics_primary.png"  # 補足図 1（集計 world_aggregate が描く診断図）
        if diag.is_file():
            shutil.copyfile(diag, OUTDIR / "si_fig_diagnostics.png")
        v35 = _json(REPO / "docs" / "verification" / "v35_name_class_contrast.json")
        m |= contrast_macros(v35)
    write_macros(m)
    write_source_data()
    digest = hashlib.sha256((OUTDIR / "summary_numbers.tex").read_bytes()).hexdigest()[:12]
    print(f"macros: {len(m)} (summary_numbers.tex {digest}); missing: {len(MISSING)}")
    for w in MISSING:
        print(f"  missing: {w}")


if __name__ == "__main__":
    main()
