# %% KSJ W05 河川（Geoshape）をラベル水系ぶんだけ集める一回作業。
# 使い方: uv run python packages/topsem/scripts/fetch_rivers.py [--force]
# 生成: data/covariates/ksj_rivers/{10桁}.json + {6桁}_codes.json
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.covariates.rivers import BASIN_CODE, fetch_basin


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    total = 0
    for label, code in sorted(BASIN_CODE.items()):
        paths = fetch_basin(code, force=args.force)
        total += len(paths)
        print(f"{code} {label:30s} {len(paths):3d} river files")
    print(f"→ {total} files under data/covariates/ksj_rivers/")


if __name__ == "__main__":
    main()
