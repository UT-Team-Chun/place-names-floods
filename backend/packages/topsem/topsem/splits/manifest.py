"""分割 manifest（data/splits/<version>/manifest.json）の読み込み・検証。

凍結された分割のみ解析に使える（CLAUDE.md 漏洩防止 #1）。
manifest は「単位（event, watershed）→ fold」の対応そのものを持ち、
行のその場分割は行わない。
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field

from topsem.config import SPLITS_DIR


class SplitUnit(BaseModel):
    event: str
    watershed: str


class FoldDef(BaseModel):
    name: str
    holdout: list[SplitUnit] = Field(default_factory=list)


class SplitSpec(BaseModel):
    version: str
    kind: Literal["watershed_holdout", "event_holdout", "temporal_holdout", "strict"]
    created_utc: str = ""
    event_dates: dict[str, str] = Field(default_factory=dict)  # event -> YYYY-MM-DD
    folds: list[FoldDef]
    notes: str = ""
    geo_exclusion_m: float = 0.0  # A-7: watershed系のみ有効。隣接窓の幾何排他半径

    @property
    def units(self) -> list[SplitUnit]:
        seen: dict[tuple[str, str], None] = {}
        for f in self.folds:
            for u in f.holdout:
                seen.setdefault((u.event, u.watershed), None)
        return [SplitUnit(event=e, watershed=w) for (e, w) in seen]

    def fold_mask(self, fold: FoldDef, events: list[str], watersheds: list[str]) -> list[bool]:
        hold = {(u.event, u.watershed) for u in fold.holdout}
        return [(e, w) in hold for e, w in zip(events, watersheds, strict=True)]

    def train_mask(self, fold: FoldDef, events: list[str], watersheds: list[str]) -> list[bool]:
        """fold の訓練側マスク（他の分割系では holdout の補集合）。

        temporal では**評価イベントより後のイベントは訓練からも除外**する
        （両側除外 — ~holdout の素朴補集合は未来漏洩。2026-09-23 レビュー A-1）。
        """
        hold = {(u.event, u.watershed) for u in fold.holdout}
        if self.kind != "temporal_holdout":
            return [(e, w) not in hold for e, w in zip(events, watersheds, strict=True)]
        eval_dates = {self.event_dates.get(e, "") for (e, _) in hold if e}
        d = min((x for x in eval_dates if x), default="")
        return [
            ((e, w) not in hold) and bool(self.event_dates.get(e, "")) and self.event_dates[e] < d
            for e, w in zip(events, watersheds, strict=True)
        ]

    def _evaluable(self, units: set[tuple[str, str]]) -> set[tuple[str, str]]:
        """temporal では『より前のイベントが 1 つもない』 unit は評価不能で正当。

        他の分割系では全 unit が評価対象。
        """
        if self.kind != "temporal_holdout":
            return units
        return {
            (e, w)
            for (e, w) in units
            if self.event_dates.get(e)
            and any(d and d < self.event_dates[e] for d in self.event_dates.values())
        }

    def validate(self, events: list[str], watersheds: list[str]) -> dict[str, object]:
        """漏洩チェック。全 unit がちょうど 1 fold の holdout に現れること等。"""
        units = set(zip(events, watersheds, strict=True))
        all_hold = {(u.event, u.watershed) for f in self.folds for u in f.holdout}
        # fold 内に同一キーの unit が複数あるレイヤ（同流域を束ねる sub-basin unit）は
        # 正常。**別々の fold** に現れる才是真正的違反なので、fold 単位 set で数える。
        counts: dict[tuple[str, str], int] = {}
        for f in self.folds:
            for k in {(u.event, u.watershed) for u in f.holdout}:
                counts[k] = counts.get(k, 0) + 1
        problems = []
        missing = self._evaluable(units) - all_hold
        dup = {u for u, c in counts.items() if c > 1}
        extra = all_hold - units
        if missing:
            problems.append(f"{len(missing)} units never held out")
        if dup:
            problems.append(f"{len(dup)} units held out in >1 fold")
        if extra:
            problems.append(f"{len(extra)} holdout units absent from data")
        temporal = self.holds_dates_check() if self.event_dates else []
        problems += temporal
        return {
            "ok": not problems,
            "problems": problems,
            "n_units": len(units),
            "n_units_evaluable": len(self._evaluable(units)),
            "n_units_held_out": len(all_hold & units),
            "n_folds": len(self.folds),
            "temporal_violations": temporal,
        }

    def holds_dates_check(self) -> list[str]:
        """temporal fold の時間漏洩検査（2026-09-23 書き換え・レビュー A-1）。

        訓練集合は train_mask の定義そのもの（より前のイベントのみ）なので、
        違反し得るのは fold 定義が壊れている場合のみ。非 temporal 系は日程が
        軸なので空返す（旧実装は event-loo にも素朴補集合で警告し、意味が違った）。
        """
        if self.kind != "temporal_holdout":
            return []
        out = []
        for f in self.folds:
            held_events = {u.event for u in f.holdout}
            if len(held_events) > 1:
                out.append(f"fold {f.name}: holds >1 event {sorted(held_events)}")
                continue
            out.extend(
                f"fold {f.name}: eval event {e} 版日付不明（訓練側が決定不能）"
                for e in held_events
                if not self.event_dates.get(e)
            )
        return out


def save_split(spec: SplitSpec) -> Path:
    d = SPLITS_DIR / spec.version
    d.mkdir(parents=True, exist_ok=True)
    path = d / "manifest.json"
    path.write_text(spec.model_dump_json(indent=2) + "\n", encoding="utf-8")
    return path


def load_split(version: str) -> SplitSpec:
    path = SPLITS_DIR / version / "manifest.json"
    if not path.is_file():
        msg = f"分割 {version} の manifest 無し（{path}）。先に freeze_split を実行"
        raise FileNotFoundError(msg)
    return SplitSpec.model_validate_json(path.read_text(encoding="utf-8"))


def make_watershed_holdout(
    units: list[tuple[str, str]],
    version: str,
    kind: Literal["watershed_holdout", "event_holdout"] = "watershed_holdout",
    event_dates: dict[str, str] | None = None,
    notes: str = "",
) -> SplitSpec:
    """節ごとの leave-one-out（この研究の既定: watershed = 水系/流域）。"""
    if kind == "watershed_holdout":
        keys = sorted({w for _, w in units})
        folds = [
            FoldDef(
                name=f"hold_{w}",
                holdout=[SplitUnit(event=e, watershed=w2) for (e, w2) in sorted(units) if w2 == w],
            )
            for w in keys
        ]
    else:
        keys = sorted({e for e, _ in units})
        folds = [
            FoldDef(
                name=f"hold_{e}",
                holdout=[SplitUnit(event=e2, watershed=w) for (e2, w) in sorted(units) if e2 == e],
            )
            for e in keys
        ]
    return SplitSpec(
        version=version,
        kind=kind,
        event_dates=event_dates or {},
        folds=folds,
        notes=notes,
    )
