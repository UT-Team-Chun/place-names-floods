# Amendments B-1 and B-4 — corrections of implementation errors (English translation)

> English translation (2026-10-27) of `docs/protocol/amendments/B-1_s6_group_buffer_fix.md` and
> `docs/protocol/amendments/B-4_lexicon_norway_yaml_fix.md`. The Japanese files are the record; where the two differ,
> the Japanese files govern. Translator's notes are in [square brackets].

---

## B-1 — Correction and rerun of an implementation defect in sensitivity analysis §7-6 (places merged with a 10 km buffer) (2026-10-26)

Status: **done, rerun, approved (2026-10-26; the author reports that the co-author has confirmed it)**. A fix that
aligns the implementation with the registered definition, not a change of conditions.
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924).

### Before (implementation at the freeze, commit e625924 / image 1792961983094088000)

§7-6 was registered as "location groups with boundary width 10 km → 20 km (buffer 5 → 10 km) —
`--group-buffer-km 10`". The implementation (`world_eval.py`)

- computed the cell cap (10,000 per location group, `quotas_for`) with **the groups of the 10 km buffer**, but
- computed the fold assignment and the bootstrap units (`folds_for`) with **the groups of the 5 km buffer (195)**.

Cause: `build_table` passed the buffer width in `DataFrame.attrs["buffer_km"]`, and `pd.concat(blocks, axis=1)` in
`add_channels` dropped the attrs (in pandas 3.0.6, when the attrs of the inputs differ), so `folds_for` fell back
to the default 5 km. With the 10 km buffer there are 166 location groups (195 with 5 km).

### After

- The buffer width is passed to `folds_for` as an argument (`buffer_km=sens.buffer_km`); `DataFrame.attrs` is not
  relied on. Regression test `test_folds_for_uses_the_given_group_buffer` (the argument of 10 km applies even if 5 km
  remains in attrs)
- §7-6 is rerun with the corrected code (run tag `world_v2_s6_buffer_b1`, 10 folds, 4 CPU, 6 GiB). No other analysis
  is rerun

### Scope

- **§7-6 only.** The primary, country, time and llmcut splits and the other sensitivity analyses (§7-1 to 5 and 7)
  all use the 5 km buffer (the default), so even when the attrs are lost they return to the same value and the
  results do not change
- `DataFrame.attrs["split"]` (the branch that selects folds 0/1 in the time splits) is set after the same concat and
  is not lost. Even if it were lost, the time splits evaluate fold 0 only and select the same rows as the default path
  (training = fold ≥ 0 and ≠ 0 = fold 1), so the results do not change (checked)

### Is this a change made after seeing results?

**Yes.** The defect was found after the aggregate of s6 had been seen (equivalent-null, R = −0.004 [−0.011, 0.002]),
from the number of bootstrap units being 195, the same as the primary split. The fix aligns the implementation with
the registered definition and does not change the decision rules, thresholds or targets. The original run
(`world_v2_s6_buffer`) is kept as "the run with an implementation defect"; the paper's supplementary table reports the
corrected run as the result of §7-6 and notes the original run.

### Result of the rerun (2026-10-26)

`world_v2_s6_buffer_b1` (commit 356ecaf, image 1792969492729423000): 166 location groups (as intended), 1,292,256
cells, terrain increment 8.73e-03, R = +0.0003 [−0.0044, +0.0053] (90% [−0.0037, +0.0046]) → **equivalent-null**
(the defective original run was also equivalent-null, R = −0.0042). The verdict does not change.

---

## B-4 — Correction and rerun of an implementation defect: Norway (NO) was not read in the lexicon's country → language mapping (2026-10-26)

Status: **registered (2026-10-26 evening, before the rerun with the corrected code)**. The frozen rules are not
changed. The code is corrected so that the frozen lexicon (`docs/protocol/v2/lexicon_hydro_v1.yaml`, sha256 as at the
freeze) is read as intended. Treated like B-1 (the implementation defect of sensitivity §7-6): both the defective and
the corrected run are reported.
Protocol: `preregistration_v2.md` (FROZEN 2026-10-25T20:58:37Z, commit e625924).

### Defect

- The country → language mapping of the frozen lexicon contains the unquoted key `NO: ["no"]`. PyYAML's (YAML 1.1)
  `safe_load` reads unquoted `NO` as the Boolean false, and `topsem.world.lexicon.Lexicon` turns the key into
  `"FALSE"` with `str(cc).upper()`. The languages of country code `NO` became empty, and the lexicon never matched
  a Norwegian place name
- Scanning the whole lexicon, this is the only key or value that is not a string (the language block `no` itself is
  read correctly)
- Impact: in the sample of the confirmatory analysis's primary split, 33,813 cells (11 AOIs) have a Norwegian
  nearest settlement. With the intended mapping, the nearest name of 5,197 of them (15.4%, 0.34% of all cells)
  matches the lexicon. The lexicon features (and their AOI means, the regional composition) were 0 in every run.
  The relation with the flood labels was not looked at in this check
- How it was found: the blind review (blind round 1, referee 1) pointed it out after reading the lexicon and code in
  the public repository. The SI statement that "Norway is the only country whose language has no lexicon entries"
  misread this defect as a gap in the lexicon

### Fix

- The YAML loading of `Lexicon` is changed to a reader that, like YAML 1.2, treats only `true` / `false` as Booleans
  (`yes` / `no` / `on` / `off` are read as strings). The frozen lexicon file is not changed
- A test checks that every country key is a two-letter string and that the analysed countries map to the intended
  languages

### Rerun (with the corrected code; otherwise identical to the confirmatory analysis)

| run tag | split | arms | purpose |
|---|---|---|---|
| `world_v2_lexfix` | primary (10 folds), country (10), temporal (1), llmcut (1, with N-read) | default (same as the confirmatory analysis) | the confirmatory analysis with the corrected reading of the lexicon |

- Aggregation with `world_aggregate` (the decision table of §6). The yardstick robustness of v36 is also computed for
  the primary split
- C2 of B-3 (GFM, refitting) is done with the corrected code, and C1 (rescoring) with the held-out predictions of
  both world_v2 and world_v2_lexfix (an addition to the rules of B-3; both decided before the GFM results are seen)
- The sensitivity analyses of §7 and the post hoc runs of B-2 are rerun with the corrected code if the primary
  verdict changes. If it does not, the original runs are reported, noting that the defect affects 0.34% of all cells
- The power analysis (before the freeze, synthetic labels) put the synthetic effect on the lexicon's water-side
  indicator, so the Norwegian cells received no effect. As a design-stage analysis it is reported as it is, with this
  point noted

### Reporting rules

- The defect and the fix are described in the text and SI, and the values of the confirmatory analysis are given for
  both the defective run (the output of the frozen code) and the corrected run. Which is written as primary: if the
  verdicts are the same, the run with the frozen code is primary (the result of the procedure at the freeze) and the
  corrected run is given alongside. If the verdicts differ, both are written with equal weight and the difference is
  stated in the abstract
- After seeing the results, the choice of which run is primary is not revisited
