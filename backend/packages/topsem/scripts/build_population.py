"""解析母集団（一般地名を含む）を構築して data/population/ に出す。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/build_population.py \
        --input ../disaster-toponyms/data/raw/geolonia/latest.csv [--name-col name]

入力 CSV がない場合は app リポジトリ側で Geolonia 住所テーブルを取得しておく
（docs/data/data_manifest.md 参照）。
"""

from __future__ import annotations

import argparse
import gzip
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from topsem.config import DATA_DIR
from topsem.corpus.population import (
    NameIndex,
    infer_name_column,
    population_stats,
    tag_rows,
)
from topsem.corpus.snapshot import load_snapshot
from topsem.runs import write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="全国地名 CSV（Geolonia 等）")
    parser.add_argument("--name-col", default=None, help="地名列名（省略時推定）")
    parser.add_argument("--tag", default=None)
    parser.add_argument("--app-repo", default=None)
    args = parser.parse_args()

    tag = args.tag or datetime.now(UTC).strftime("pop_%Y%m%dT%H%M%SZ")
    snapshot = load_snapshot(Path(args.app_repo) if args.app_repo else None)
    df = pd.read_csv(args.input)
    name_col = args.name_col or infer_name_column(df)
    tagged = tag_rows(df, name_col, NameIndex.from_snapshot(snapshot))
    stats = population_stats(tagged)

    out_dir = DATA_DIR / "population" / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_dir / "population.csv.gz", "wt", encoding="utf-8") as f:
        tagged.to_csv(f, index=False)
    write_json(
        out_dir / "summary.json",
        {
            "tag": tag,
            "created_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "input": str(args.input),
            "name_col": name_col,
            "corpus": snapshot.manifest(),
            "stats": stats,
        },
    )
    print(f"== 母集団構築 ({out_dir}) ==")
    for k, v in stats.items():
        print(f"{k:24s}: {v}")


if __name__ == "__main__":
    main()
