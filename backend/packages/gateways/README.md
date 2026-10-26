# Gateways Package

外部システム接続パッケージ。
バックエンドとジョブアプリケーション間で共有される外部システムとの接続・操作を提供します。

> **📖 このREADMEについて**
>
> このファイルには、**backend/packages/gateways/ のディレクトリ構造**、**開発コマンド**、**実装パターン**、**コード規約**が記載されています。
>
> - **SWEの作業手順やDB操作の追加方法** → [docs/roles/software_engineer.md](../../../docs/roles/software_engineer.md) を参照してください
> - **アーキテクチャ設計の詳細** → [backend/docs/architecture.md](../../docs/architecture.md) を参照してください
> - **依存関係図** → [backend/docs/dependencies.mmd](../../docs/dependencies.mmd) を参照してください

---

## 📂 ディレクトリ構造

```
backend/packages/gateways/
└── gateways/                # パッケージ本体
    ├── rdb/                 # RDB接続
    │   ├── connections.py   # DB接続管理
    │   ├── models/          # SQLAlchemyモデル (ORM)
    │   └── operations/
    │       └── repositories/  # Repository実装 (CRUD操作)
    │
    ├── llm/                 # LLM接続 (OpenAI等)
    │   ├── connections.py   # クライアント生成
    │   ├── models/          # データモデル
    │   └── operations/      # LLM操作
    │
    └── aws/                 # AWS接続 (S3等)
        ├── connections.py   # クライアント生成
        ├── models/          # データモデル
        └── operations/      # AWS操作
```

---

## 🎯 各サービスの役割

### `rdb/` - RDB接続・データベース操作

**Repository Pattern** を採用し、SQLAlchemyを使用します。

#### connections.py
- `get_sync_session()`: 同期セッションを取得 (FastAPIの依存性注入で使用)
- `get_async_session()`: 非同期セッションを取得 (将来の拡張用)
- Alembicによるマイグレーション管理

#### models/
- SQLAlchemyモデル (ORM)
- `Base`: すべてのORMモデルの基底クラス
- テーブル定義 (例: `User`, `Project` など)

#### operations/repositories/
- CRUD操作の実装
- **重要**: ORMオブジェクトをそのまま返さず、必ずPydanticモデルに変換して返す
- トランザクション管理は呼び出し側 (usecase) の責任

---

### `llm/` - LLM接続 (OpenAI等)

#### connections.py
- `get_openai_client()`: OpenAIクライアントの取得

#### operations/
- LLM操作の実装 (テキスト生成、埋め込み生成など)

---

### `aws/` - AWS接続 (S3等)

#### connections.py
- `get_s3_client()`: S3クライアントの取得

#### operations/
- AWS操作の実装 (S3アップロード、ダウンロードなど)

---

## ⚠️ 重要なルール

### 1. ORMオブジェクトの流出禁止

RepositoryはORMオブジェクトをそのまま返さず、必ずPydanticモデルに変換して返してください。

### 2. トランザクション管理

Repositoryはトランザクション管理を行いません（`commit`, `rollback`しない）。

- Repositoryでは必要に応じて`flush()`のみ使用 (IDを取得する場合など)
- トランザクション管理は呼び出し側 (usecase) の責任

### 3. セッション管理

FastAPIの`Depends(get_sync_session)`で依存性注入します。

---

## 🔍 よくある質問

### Q. 新しいRepositoryを追加したい

**A.** `gateways/rdb/operations/repositories/` にファイルを作成してください。

### Q. マイグレーションを実行したい

**A.** Alembicを使用してスキーマ変更を管理します:

```bash
# マイグレーションファイル生成
alembic revision --autogenerate -m "Add users table"

# マイグレーション実行
alembic upgrade head
```

### Q. LLMを使いたい

**A.** `gateways.llm`を使用してください。

### Q. AWSを使いたい

**A.** `gateways.aws`を使用してください。

---

## 📚 参考ドキュメント

- **プロジェクト全体の仕様**: [.claude/CLAUDE.md](../../.claude/CLAUDE.md)
- **アーキテクチャ詳細**: [docs/architecture.md](../../docs/architecture.md)
- **依存関係図**: [docs/dependencies.mmd](../../docs/dependencies.mmd)
