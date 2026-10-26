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
        body = "\n".join(line for line in body.splitlines() if not line.lstrip().startswith("%"))
        return expand_inputs(body, depth + 1)

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
    text = expand_inputs((HERE / "main.tex").read_text(encoding="utf-8"))
    text = substitute(text, m)
    left = sorted(set(re.findall(r"\\missing\{[^}]*\}", text)))
    unknown = sorted({k for k in re.findall(r"\\([A-Z][A-Za-z]+)", text) if k in m})
    if left or unknown:
        print("unresolved:", left[:10], unknown[:10])
        return 1
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT.relative_to(HERE)} ({len(text):,} chars, {len(m)} macros)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
