"""Gate 1 監査: コーパスの成立性を数えて data/runs/ に残す。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/audit_corpus.py [--tag audit_YYYYMMDD]

モデルより必ず先にこれを回す（CLAUDE.md）。
"""

from __future__ import annotations

import argparse
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.corpus.audit import audit_built_geojson, audit_corpus
from topsem.corpus.snapshot import load_snapshot
from topsem.runs import start_run, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=None, help="run 名（省略時: audit_<UTC 日時>）")
    parser.add_argument("--app-repo", default=None, help="disaster-toponyms のパス")
    args = parser.parse_args()

    tag = args.tag or datetime.now(UTC).strftime("corpus_audit_%Y%m%dT%H%M%SZ")
    snapshot = load_snapshot(Path(args.app_repo) if args.app_repo else None)
    result = audit_corpus(snapshot)
    geojson = audit_built_geojson(snapshot)

    run = start_run(tag)
    summary = {
        **run.provenance(),
        "corpus": snapshot.manifest(),
        "audit": result,
        "built_geojson": geojson,
    }
    write_json(run.path("summary.json"), summary)

    print(f"== Gate 1 監査 ({run.run_dir}) ==")
    print(f"app SHA        : {snapshot.app_sha[:12]}")
    print(f"地名レコード   : {result['n_toponyms']}")
    print(f"  洪水系を含む : {result['n_toponyms_flood_any']}")
    print(f"  反証あり     : {result['n_toponyms_with_dispute']}")
    print(f"  最大 level   : {result['toponym_evidence_max_level_hist']}")
    print(f"  位置あり     : {result['toponyms_with_location']} {result['location_precision']}")
    print(f"  別名・状態   : {result['toponyms_by_hazard_type']} {result['toponyms_by_status']}")
    print(
        f"辞書要素       : {result['n_elements']}（洪水系 {len(result['element_ids_flood_any'])}, 替字 {result['n_elements_substitution']}）"
    )
    if geojson:
        print(
            f"ビルド版       : {geojson['n_features']} features / "
            f"洪水系 {geojson['n_flood_any']} / 洪水型なし {geojson['n_without_hazard_type']} / "
            f"反証フラグ {geojson['n_disputed']} / level 分布 {geojson['evidence_level_hist']}"
        )
    else:
        print("ビルド版       : geojson 無し（app 側 make build-data 後に再実行可能）")


if __name__ == "__main__":
    main()
