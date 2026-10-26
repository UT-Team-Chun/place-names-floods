"""N-read の入力: AOI のセルに割り当てられた P 名（最近傍 K=3）の一意な表。

出力 parquet の列: gid, asciiname, name, cc, country（countryInfo の国名）。
ラベルは読まない（features.parquet の割当 gid と GeoNames 部分集合だけを使う）。

    uv run python packages/topsem/scripts/nread_names.py --out ../data/global/nread/names.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from topsem.config import DATA_DIR

ROOT = DATA_DIR / "global"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "nread" / "names.parquet"))
    args = ap.parse_args()
    gids: set[int] = set()
    for f in sorted((ROOT / "labels").glob("EMSR*/features.parquet")):
        t = pq.read_table(f, columns=["P_gid1", "P_gid2", "P_gid3"]).to_pandas()
        g = np.unique(t.to_numpy().ravel())
        gids.update(int(x) for x in g if x >= 0)
    gn = pd.read_parquet(ROOT / "geonames" / "aoi_subset.parquet")
    gn = gn[gn["gid"].isin(gids)][["gid", "asciiname", "name", "cc"]].copy()
    names = {}
    for line in (ROOT / "geonames" / "countryInfo.txt").read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            f = line.split("\t")
            names[f[0]] = f[4]
    gn["country"] = gn["cc"].map(names).fillna(gn["cc"])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    gn.sort_values("gid").to_parquet(out, index=False)
    print(f"{len(gn):,} unique P names → {out}")


if __name__ == "__main__":
    main()
