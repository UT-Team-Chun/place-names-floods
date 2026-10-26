"""行政 GIS 配布 ZIP の展開ヘルパー。

日本の行政データ ZIP はエントリ名が Shift-JIS（flag_bits の UTF-8 ビット無し）で
入っており、そのまま展開すると macOS/Linux で文字化け・展開不能になる。
cp437→shift_jis 復号で安全に展開する（失敗時は原名のまま）。
"""

from __future__ import annotations

import pathlib
import zipfile


def extract_sjis_zip(zip_path: pathlib.Path, dest_dir: pathlib.Path) -> list[str]:
    """ZIP を dest_dir に展開し、展開後相対パス名の一覧を返す。"""
    out: list[str] = []
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            name = _fix_name(info)
            target = dest_dir / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(zf.read(info))
            out.append(name)
    return out


def _fix_name(info: zipfile.ZipInfo) -> str:
    if info.flag_bits & 0x800:  # UTF-8 指定済み
        return info.filename
    try:
        return info.filename.encode("cp437").decode("shift_jis")
    except (UnicodeEncodeError, UnicodeDecodeError):
        return info.filename
