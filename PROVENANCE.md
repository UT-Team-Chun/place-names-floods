# Provenance of the frozen protocol

This repository is a curated copy of the project's development repository, prepared for publication.
Its history starts at the copy; the development repository is not public. The freeze is documented by
the commit identifiers below, by the content checksums that let anyone verify that the files here are
the frozen files, and by the job records of the computing cluster described in the paper.

## Freeze

- Protocol `docs/protocol/preregistration_v2.md` frozen at **2026-10-25T20:58:37Z**: freeze-target commit
  `e62592444170907104d8375473d40b41cc053097`, freeze commit `1c1faf6841eef2b732df552ee31622a81f422fbc`,
  tag `prereg-v2-frozen` (development repository).
- Before the freeze, no model had been fitted to the observed flood labels. The labels were used only as a
  label-only incidence table and as AOI-level flooded fractions to which the synthetic labels of the pilots
  and the power analysis were calibrated. Through the latter, the pilots showed that the regional
  composition of names (AOI means of the name features) tracks AOI-level flood rates; this association
  between real flood rates and real names, seen before the freeze, is why the primary estimand holds the
  composition fixed (see the paper's Methods and Supplementary Note 3).
- The freeze was recorded in version control and in the job records of the computing cluster, which ran
  the confirmatory analysis from a container image built after the freeze. **The protocol was not deposited
  with a third-party registry before the analysis**; the paper therefore describes it as pre-specified.

## Verifying the frozen content

The files in this repository are byte-identical to their frozen versions:

| File | sha256 |
|---|---|
| `docs/protocol/preregistration_v2.md` (as of tag `prereg-v2-frozen`, unchanged since) | `4d9870ad888ebf88457cb0c8c1f2283264168fbab59ebad3db9599008bfd6dde` |
| `docs/protocol/v2/lexicon_hydro_v1.yaml` | `6c9999921956ba13f929b18ef921afcd0e05334293ef0989d9357c6a8d59c58d` |
| `docs/data/ems_inventory_20261025.json` | `f5c55eff3e470bce1fee7ae228641d1d3b21361f30b6d5ac9b87eef438028699` |

```bash
shasum -a 256 docs/protocol/preregistration_v2.md docs/protocol/v2/lexicon_hydro_v1.yaml \
  docs/data/ems_inventory_20261025.json
```

The digests of the other frozen inputs (GeoNames dump, coastline, incidence table, language-model
outputs, power-analysis runs) are listed in §10 of the protocol and in
`docs/verification/v33_freeze_record.json`.

## Changes after the freeze

- `docs/protocol/amendments/B-1_s6_group_buffer_fix.md`: an implementation defect in one sensitivity analysis
  (the location-group buffer did not reach the fold assignment) was found after its first result, corrected
  and re-run; the verdict did not change. Both runs are reported in the paper.
- `docs/protocol/amendments/B-2_posthoc_review_analyses.md` and its four addenda: post hoc analyses added after
  the confirmatory results in response to five rounds of internal review (alternative yardsticks and a two-way
  bootstrap, recalibration-only reference models, grouping by flood event, reassigned folds, a label-free positive
  control, distance bands, results by class of name, within-area odds ratios, the mapping sensor, effect-size
  translation, power for the secondary splits). Most were registered, with their computations and reporting rules,
  before they were run; the paper marks those that were not. All are labelled post hoc and change no pre-specified
  verdict. Scripts: `backend/verification/r_okauchi/v36`–`v40`; outputs: `docs/verification/`.

The confirmatory run records label the protocol "DRAFT until frozen" because the code wrote a fixed string
there; the status string is now read from the protocol's status line. The power-analysis runs record the
lexicon under its draft version label (sha256 prefix `4f9203535d64`); its entries are identical to the frozen
lexicon, which differs only in a status comment and the version label.
