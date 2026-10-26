# gateways/

ゲートウェイレイヤー。外部システム（DB, AWS, OpenAI等）との接続と操作を担当します。

## 役割

- 外部システムとの接続管理
- データの永続化・取得
- 外部APIの呼び出し
- ORMモデルとドメインモデルの変換

## ディレクトリ構成

```
gateways/
├── connections/            # クライアント生成
│   ├── aws_client.py       # AWS (boto3) クライアント
│   ├── openai_client.py    # OpenAI クライアント
│   └── rdb.py              # RDB (SQLAlchemy) セッション
├── db_models/              # ORM定義
│   ├── base.py             # SQLAlchemy Base
│   ├── audit_log.py        # 監査ログテーブル
│   └── ...                 # 各テーブルのORMモデル
└── operations/             # 実際の操作実装
    ├── repositories/       # CRUD操作（Write/Read）
    ├── queries/            # 複雑な読み取り（JOIN, 集計）
    └── llm_operations.py   # LLM関連の操作
```

## 各ディレクトリの役割

### `connections/` (クライアント生成)
- 外部システムへの接続クライアントの生成
- セッション管理
- 接続プールの管理

### `db_models/` (ORM定義)
- SQLAlchemyによるテーブル定義
- データベーススキーマの表現

### `operations/` (操作実装)
- **`repositories/`**: 単純なCRUD操作
  - 例: `get()`, `create()`, `update()`, `delete()`
  - **Commitはしない**（Usecaseが責任を持つ）
- **`queries/`**: 複雑な読み取り操作
  - 例: JOIN、集計、複雑な検索条件
  - 戻り値は`domain_models`または`schema`に変換
- **`llm_operations.py`**: LLM関連の操作
  - OpenAI APIの呼び出しなど

## 実装ルール

### ✅ やるべきこと

- **ドメインモデルへの変換**: ORMオブジェクトをドメインモデルに変換して返す
- **エラーハンドリング**: 外部システムのエラーを適切にハンドリング
- **接続管理**: 接続プールやセッションの適切な管理

### ❌ やってはいけないこと

- **Repositoryでのcommit**: `operations/repositories/`では`commit()`しない
- **ORMオブジェクトの流出**: ORMオブジェクトを直接返さない
- **ビジネスロジック**: ゲートウェイ層にビジネスロジックを書かない

## 使用例

### connections/rdb.py
```python
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from config.config import Config

config = Config()
engine = create_engine(config.database.url)
SessionLocal = sessionmaker(bind=engine)

def get_db() -> Session:
    """DBセッションの取得"""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
```

### db_models/user.py
```python
from sqlalchemy import Column, String, Integer
from gateways.db_models.base import Base

class UserModel(Base):
    """ユーザーテーブル"""
    __tablename__ = "users"

    id = Column(String, primary_key=True)
    name = Column(String, nullable=False)
    age = Column(Integer, nullable=False)
```

### operations/repositories/user_repository.py
```python
from sqlalchemy.orm import Session
from gateways.db_models.user import UserModel
from domain_models.user_domain import UserDomain

class UserRepository:
    def __init__(self, db: Session):
        self.db = db

    def get(self, user_id: str) -> UserDomain:
        """ユーザー取得"""
        user = self.db.query(UserModel).filter_by(id=user_id).first()
        if not user:
            raise ValueError(f"User not found: {user_id}")

        # ORMモデルをドメインモデルに変換
        return UserDomain(
            user_id=user.id,
            name=user.name,
            age=user.age,
        )

    def create(self, domain: UserDomain) -> None:
        """ユーザー作成"""
        user = UserModel(
            id=domain.user_id,
            name=domain.name,
            age=domain.age,
        )
        self.db.add(user)
        # ※ commitはしない（Usecaseが責任を持つ）

    def update_score(self, user_id: str, score: float) -> None:
        """スコア更新"""
        user = self.db.query(UserModel).filter_by(id=user_id).first()
        if user:
            user.score = score
        # ※ commitはしない
```

### operations/queries/user_query.py
```python
from sqlalchemy.orm import Session
from gateways.db_models.user import UserModel
from domain_models.user_domain import UserDomain

class UserQuery:
    def __init__(self, db: Session):
        self.db = db

    def find_adult_users(self) -> list[UserDomain]:
        """成人ユーザーの検索"""
        users = self.db.query(UserModel).filter(UserModel.age >= 20).all()

        # ORMモデルをドメインモデルのリストに変換
        return [
            UserDomain(
                user_id=u.id,
                name=u.name,
                age=u.age,
            )
            for u in users
        ]
```

### operations/llm_operations.py
```python
from openai import OpenAI
from config.config import Config

class LLMOperations:
    def __init__(self):
        config = Config()
        self.client = OpenAI(api_key=config.llm.openai_api_key)

    def generate_text(self, prompt: str) -> str:
        """テキスト生成"""
        response = self.client.chat.completions.create(
            model="gpt-4",
            messages=[{"role": "user", "content": prompt}]
        )
        return response.choices[0].message.content
```

## 依存関係

```
features → gateways
gateways → domain_models (ドメインモデルへの変換のため)
gateways → config (接続情報の取得のため)

❌ gateways → features
❌ gateways → alg
```
