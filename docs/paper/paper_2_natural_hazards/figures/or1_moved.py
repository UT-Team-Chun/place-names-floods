"""Online Resource 1 の Text S4・S5 が読む表の写し（2026-11-07、本文の組み直し）。

本文から Online Resource 1 へ移した表（分割ごとの比較、感度分析、事後の区間の上端）は make_figures.py が本文用に
生成する（app_table_splits.tex・app_table_sensitivity.tex・tab_upper.tex）。その写しを or1_*.tex に作り、
Online Resource 1 の中で通るように参照を直す: 本文のラベルは xr-hyper の接頭辞 P- を付けて「of the paper」を足し、
自分自身を指す「of Online Resource~1」を外す。数値には触れない。make supplement が組版の前に走らせる。
"""

from __future__ import annotations

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
MOVED = {"app_table_splits.tex": "or1_app_table_splits.tex",
         "app_table_sensitivity.tex": "or1_app_table_sensitivity.tex",
         "tab_upper.tex": "or1_tab_upper.tex"}
# Online Resource 1 の中で定義されるラベル（本文には無い）
OR1_LABELS = {"tab:si-splits", "tab:si-sensitivity", "tab:upper", "tab:deviations"}
_REF = re.compile(r"(Sects?\.|Tables?|Figs?\.)~\\ref\{((?:sec|tab|fig|app):[^}]*)\}([a-z]*)")
_BARE = re.compile(r"\\ref\{((?:sec|tab|fig|app):[^}]*)\}")  # 種類の語の無い参照（and~\ref{...} など）


def to_or1(tex: str) -> str:
    """本文向けの参照を Online Resource 1 向けに直す（パネルの文字 \\ref{fig:main}b はその後に of the paper）。"""

    def ref(m: re.Match) -> str:
        kind, label, panel = m.group(1), m.group(2), m.group(3)
        if label in OR1_LABELS:
            return m.group(0)
        return f"{kind}~\\ref{{P-{label}}}{panel} of the paper"

    tex = _REF.sub(ref, tex)
    def bare(m: re.Match) -> str:
        return m.group(0) if m.group(1) in OR1_LABELS else f"\\ref{{P-{m.group(1)}}}"

    tex = _BARE.sub(bare, tex)
    tex = re.sub(r"Eq\.~\(\\ref\{(eq:[^}]*)\}\)", r"Eq.~(\\ref{P-\1}) of the paper", tex)
    return re.sub(r"\s+of\s+Online\s+Resource~1", "", tex)  # 改行をまたぐ「of Online\nResource~1」も


def main() -> None:
    """移した表の写しを書く。"""
    for src, dst in MOVED.items():
        s = (HERE / src).read_text(encoding="utf-8")
        head = f"% 生成: figures/or1_moved.py（{src} の写し、Online Resource 1 用）\n"
        (HERE / dst).write_text(head + to_or1(s), encoding="utf-8")
        print(f"wrote figures/{dst}")


if __name__ == "__main__":
    main()
