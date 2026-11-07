"""投稿用の 1 ファイルを作る: \\input を展開し、figures/numbers.tex のマクロを値に置き換え、文献（main.bbl）を中に入れ、
図を本文で出る順に Fig1.pdf, Fig2.pdf, ... として submission/ に写す（Springer の指示: LaTeX は 1 つの .tex で、図の
ファイル名は Fig と番号）。展開後に \\missing や未定義のマクロ、AUTHOR-GATED の印が残っていれば止まる。

使い方（docs/paper/paper_2_natural_hazards/ から）: python submission/flatten.py  → submission/manuscript.tex
（先に make で main.bbl を作っておく。改訂 B-19 で NHESS の版から写した）
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "submission" / "manuscript.tex"
NEWCMD = re.compile(r"^\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$")
INPUT = re.compile(r"\\input\{([^}]+)\}")


def macros() -> dict[str, str]:
    out: dict[str, str] = {}
    for line in (HERE / "figures" / "numbers.tex").read_text(encoding="utf-8").splitlines():
        m = NEWCMD.match(line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


COMMENT = re.compile(r"(?<!\\)%.*$")


def strip_comment(line: str) -> str:
    """行の中の注釈（著者向けの覚え書きなど）を消し、% だけを残す（行末の空白の扱いは変えない）。"""
    return COMMENT.sub("%", line)


def expand_inputs(text: str, depth: int = 0) -> str:
    if depth > 5:
        msg = "\\input nesting too deep"
        raise SystemExit(msg)

    def rep(m: re.Match) -> str:
        name = m.group(1)
        if name == "figures/numbers":
            return ""
        path = HERE / (name if name.endswith(".tex") else f"{name}.tex")
        body = path.read_text(encoding="utf-8")
        body = "\n".join(strip_comment(line) for line in body.splitlines() if not line.lstrip().startswith("%"))
        # 行末の % が続く閉じ括弧を食わないよう、改行で終える（例: \authorcontribution{\input{...}}）
        return expand_inputs(body, depth + 1) + "\n"

    return INPUT.sub(rep, text)


def substitute(text: str, m: dict[str, str]) -> str:
    # 長い名前から順に。\Name{} と \Name（後ろが英字でない）の両方を値に
    for k in sorted(m, key=len, reverse=True):
        v = m[k].replace("\\", "\\\\")
        text = re.sub(rf"\\{k}\{{\}}", v, text)
        text = re.sub(rf"\\{k}(?![A-Za-z])", v, text)
    return text


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=["doi", "nodoi"], default="doi",
                    help="doi: main.tex から submission/with_doi/ に。nodoi: main_nodoi.tex（DOI を待たずに出す版）から submission/no_doi/ に")
    variant = ap.parse_args().variant
    stem = "main" if variant == "doi" else "main_nodoi"
    outdir = HERE / "submission" / ("with_doi" if variant == "doi" else "no_doi")  # 2026-11-07: DOI 入りの版も専用の場所に
    outdir.mkdir(parents=True, exist_ok=True)
    m = macros()
    main_tex = "\n".join(strip_comment(line) for line in (HERE / f"{stem}.tex").read_text(encoding="utf-8").splitlines())
    text = expand_inputs(main_tex)
    for _ in range(5):  # マクロの値が別のマクロを含むことがある（例: 規則の文の中の数値）ので、変わらなくなるまで
        new = substitute(text, m)
        if new == text:
            break
        text = new
    bbl = HERE / f"{stem}.bbl"
    if not bbl.is_file():
        print(f"{stem}.bbl missing: build {stem}.tex first")
        return 1
    text = text.replace("\\bibliography{refs}", bbl.read_text(encoding="utf-8"))
    figs = re.findall(r"\\includegraphics(?:\[[^]]*\])?\{(figures/[^}]+\.pdf)\}", text)
    for k, f in enumerate(dict.fromkeys(figs), start=1):
        shutil.copyfile(HERE / f, outdir / f"Fig{k}.pdf")
        text = text.replace("{" + f + "}", "{" + f"Fig{k}.pdf" + "}")
    left = sorted(set(re.findall(r"\\missing\{[^}]*\}", text)))
    unknown = sorted({k for k in re.findall(r"\\([A-Z][A-Za-z]+)", text) if k in m})
    notes = re.findall(r"AUTHOR-GATED|TODO|FIXME|to be inserted", text)
    if left or unknown or notes:
        print("notes left:", len(notes)) if notes else None
        print("unresolved:", left[:10], unknown[:10])
        return 1
    out = outdir / "manuscript.tex"
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out.relative_to(HERE)} ({len(text):,} chars, {len(m)} macros)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
