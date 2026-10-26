# Project rules (analysis)

This file marks the repository root for the scripts (`next(p for p in parents if (p / "CLAUDE.md").is_file())`)
and records the rules the analysis follows. Notes on the development computing infrastructure are omitted.

- **Pre-specification.** Hypotheses, estimand, splits, decision rules, the smallest effect size of
  interest and the sensitivity analyses are fixed in `docs/protocol/preregistration_v2.md` before the
  confirmatory runs. Any later change is an amendment in `docs/protocol/amendments/` with its reason, date
  and the analyses it affects. Thresholds, primary metrics and target hazards are not changed after
  results are seen.
- **Leakage prevention.** Splits come from the protocol, never from random splitting. Anything fitted
  (name scores, second-stage models, standardisation) is fitted on the training side only. No information
  dated after a flood enters the features for that flood (see `docs/data/point_in_time.md`). The same place,
  aliases of a place, and floods of the same place at other times never sit on both sides of a split.
- **Run contract.** Every run writes `data/runs/<run>/` with a summary (provenance: input digests,
  protocol version, models, seeds, git revision), metrics with per-split breakdowns, cell-level held-out
  predictions, and diagnostics.
- **Traceability.** Every number in the paper is generated from run artefacts by
  `docs/paper/paper_2_global/figures/make_figures.py`; none is typed by hand.
