# Gateways Package Tests

このディレクトリには、gatewaysパッケージのユニットテストが含まれています。

## ディレクトリ構成

```
tests/
├── README.md                                    # このファイル
└── llm/
    └── operations/
        └── test_llm_operations.py               # LLM操作のテスト
```

## テストファイル

### [llm/operations/test_llm_operations.py](./llm/operations/test_llm_operations.py)

LLM操作クラスのユニットテスト。OpenAI APIクライアントをモックして動作を検証します。

**テスト対象**: `gateways.llm.operations.llm_operations.LLMOperations`

- `test_llm_operations_initialization()` - クラスの初期化テスト
- `test_llm_operations_generate()` - テキスト生成のテスト
- `test_llm_operations_generate_with_custom_params()` - カスタムパラメータでの生成テスト

## テストの実行

### すべてのテストを実行
```bash
cd backend/packages/gateways
uv run pytest tests/ -v
```

### 特定のテストファイルのみ実行
```bash
uv run pytest tests/llm/operations/test_llm_operations.py -v
```

### 特定のテスト関数のみ実行
```bash
uv run pytest tests/llm/operations/test_llm_operations.py::test_llm_operations_initialization -v
```

## テスト結果の確認

すべてのテストが成功すると、以下のような出力が表示されます:

```
tests/llm/operations/test_llm_operations.py::test_llm_operations_initialization PASSED
tests/llm/operations/test_llm_operations.py::test_llm_operations_generate PASSED
tests/llm/operations/test_llm_operations.py::test_llm_operations_generate_with_custom_params PASSED

============================== 3 passed in 0.23s ===============================
```

## テストの特徴

### モックの使用

これらのテストは、実際のOpenAI APIを呼び出すことなく動作を検証するため、`unittest.mock`を使用しています。

```python
with patch("gateways.llm.operations.llm_operations.get_openai_client") as mock_client:
    # モックされたクライアントで動作を検証
```

これにより:
- API keyなしでテストが実行できる
- テストが高速
- 外部サービスに依存しない

## 今後の拡張

以下のテストを追加する予定:

- `tests/rdb/` - データベース操作のテスト
- `tests/aws/` - AWS S3操作のテスト
- `tests/llm/connections/` - OpenAI接続のテスト

## 関連ファイル

- [gateways/llm/operations/llm_operations.py](../gateways/llm/operations/llm_operations.py) - テスト対象のLLM操作クラス
- [gateways/llm/connections.py](../gateways/llm/connections.py) - OpenAI接続管理
- [backend/api/tests/](../../../api/tests/) - API側の統合テスト
- [backend/packages/alg/tests/](../alg/tests/) - algパッケージのテスト
