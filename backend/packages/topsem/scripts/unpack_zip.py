"""行政 GIS ZIP（Shift-JIS エントリ名）を安全に展開する CLI。

使い方::

    uv run python packages/topsem/scripts/unpack_zip.py <zip> <dest_dir>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.util.archives import extract_sjis_zip


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zip_path", type=Path)
    parser.add_argument("dest_dir", type=Path)
    args = parser.parse_args()
    names = extract_sjis_zip(args.zip_path, args.dest_dir)
    for n in names:
        print(n)
    print(f"-- {len(names)} files -> {args.dest_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
