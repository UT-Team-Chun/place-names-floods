"""投稿用の 1 ファイルを作る: \\input を展開し、figures/numbers.tex のマクロを値に置き換える。

Copernicus の原稿は独自コマンドを定義できない（\\newcommand 不可）。開発中は数値を全てマクロで書き
（成果物由来・手打ち禁止）、投稿の直前にこのスクリプトで展開する。展開後に \\missing や未定義の
マクロが残っていれば止まる。

使い方（docs/paper/paper_2_nhess/ から）: python submission/flatten.py  → submission/manuscript.tex
"""

from __future__ import annotations

import re
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
    m = macros()
    main_tex = "\n".join(strip_comment(line) for line in (HERE / "main.tex").read_text(encoding="utf-8").splitlines())
    text = expand_inputs(main_tex)
    for _ in range(5):  # マクロの値が別のマクロを含むことがある（例: 規則の文の中の数値）ので、変わらなくなるまで
        new = substitute(text, m)
        if new == text:
            break
        text = new
    left = sorted(set(re.findall(r"\\missing\{[^}]*\}", text)))
    unknown = sorted({k for k in re.findall(r"\\([A-Z][A-Za-z]+)", text) if k in m})
    notes = re.findall(r"AUTHOR-GATED|TODO|FIXME", text)
    if left or unknown or notes:
        print("notes left:", len(notes)) if notes else None
        print("unresolved:", left[:10], unknown[:10])
        return 1
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT.relative_to(HERE)} ({len(text):,} chars, {len(m)} macros)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
