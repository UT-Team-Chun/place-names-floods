"""data/labels/ の全 unit から 3 種の分割 manifest を凍結する。

使い方（backend/ から）::

    uv run python packages/topsem/scripts/freeze_splits.py

生成: data/splits/{watershed-loro-1, event-loo-1, temporal-1}/manifest.json
既存バージョンがあっても上書きしない（version 名を上げて再実行する）。
"""

from __future__ import annotations

import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topsem.analysis.cell_features import discover_label_units
from topsem.config import SPLITS_DIR
from topsem.splits.manifest import (
    FoldDef,
    SplitSpec,
    SplitUnit,
    make_watershed_holdout,
    save_split,
)

_GEO_EXCL_FIRST_SUFFIX = 5  # A-7: v5 から watershed 系は幾何排他を伴う


def event_start_dates(all_units: list[dict]) -> dict[str, str]:
    """イベントごとの版日付最小（= 発生開始）。dict 反復順に依存しない。"""
    dates: dict[str, str] = {}
    for u in all_units:
        if u["event_date"]:
            cur = dates.get(u["event"])
            if not cur or u["event_date"] < cur:
                dates[u["event"]] = u["event_date"]
    return dates


def temporal_folds(units: list[tuple[str, str]], dates: dict[str, str]) -> list[FoldDef]:
    """評価イベントより**前**のイベントのみで訓練する fold（前が無い評価は作らない）。"""
    folds = []
    for e in sorted(dates):
        earlier = {e2 for e2, d in dates.items() if d < dates[e]}
        hold = [SplitUnit(event=e, watershed=w) for (e2, w) in sorted(units) if e2 == e]
        if earlier and hold:
            folds.append(FoldDef(name=f"time_{e}", holdout=hold))
    return folds


def main() -> None:
    all_units = discover_label_units()
    units = [(u["event"], u["watershed"]) for u in all_units]
    if not units:
        msg = "data/labels/*/derived/*/metadata.json が見つからない"
        raise FileNotFoundError(msg)
    dates = event_start_dates(all_units)
    now = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    suffix = int(sys.argv[1]) if len(sys.argv) > 1 else 1  # 1 なら既存 -1 名（上書きしない）
    wsv, evv, tpv = (
        f"watershed-loro-{suffix}", f"event-loo-{suffix}", f"temporal-{suffix}")
    for version in (wsv, evv, tpv):
        if (SPLITS_DIR / version).exists():
            print(f"skip (exists): {version} — 変更は version を上げて")
    ws = make_watershed_holdout(units, wsv, "watershed_holdout", dates)
    if suffix >= _GEO_EXCL_FIRST_SUFFIX:
        # A-7 (r4' F3): 隣接水系窓との物理地面共有を評価時マスクで除去。
        # 250 m = 1 cell 半径。感度 500/1000 m は実行時 spec 読替で扱う
        ws.geo_exclusion_m = 250.0
        ws.notes += " | v5: geo_exclusion 250 m (A-7)"
    ev = make_watershed_holdout(units, evv, "event_holdout", dates)
    tp = SplitSpec(
        version=tpv,
        kind="temporal_holdout",
        created_utc=now,
        event_dates=dates,
        folds=temporal_folds(units, dates),
        notes="time-held-out（評価より前のイベントのみで訓練）",
    )
    ws.created_utc = ev.created_utc = now
    for spec in (ws, ev, tp):
        if (SPLITS_DIR / spec.version).exists():
            continue
        p = save_split(spec)
        print(f"{spec.version:18s} {len(spec.folds)} folds -> {p}")


if __name__ == "__main__":
    main()
