# topsem — 研究の計算ロジック

地名の環境記憶と洪水ハザード評価（[README](../../../README.md)）の解析本体。
アプリ（disaster-toponyms）は**読取専用・git SHA 固定**で参照し、コピーしない。

## モジュール対応表（研究の概念 → コード）

| 研究の概念（計画 § / 事前登録 §） | モジュール |
|---|---|
| コーパス来歴・固定（§4.3） | `topsem.corpus.snapshot` / `audit` / `population` |
| 観測ラベルとマスク（§6.3） | `topsem.labels.masks` |
| ベースライン共変量・時点管理（§5.1） | `topsem.covariates` + `docs/data/point_in_time.md` |
| 特徴経路 A / B / C（§5.2） | `topsem.features` |
| 空間単位・凍結分割（§5.3） | `topsem.splits` + `data/splits/` |
| M0〜M5（§5.4） | `topsem.models` |
| 指標（Brier 主評価量・capture 等） | `topsem.evaluation.metrics` |
| 実行ディレクトリ契約 | `topsem.runs` |

## スクリプト（`backend/` から）

```bash
uv run python packages/topsem/scripts/audit_corpus.py       # Gate 1 監査（最初に必ず）
uv run python packages/topsem/scripts/build_population.py \
    --input ../disaster-toponyms/data/raw/geolonia/latest.csv   # 解析母集団
```

`TOPSEM_APP_REPO` でアプリ側リポジトリの位置を変更可（デフォルト: 姉妹ディレクトリ）。

## sandbox

検証用スクラッチは `sandbox/{user}/{n|s}_{YYYYMMDD}_{topic}.py`。
有益なスクラッチは必ず `scripts/` か `topsem/` に昇格させる（CLAUDE.md）。
