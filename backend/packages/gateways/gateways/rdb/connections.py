from collections.abc import Generator

from sqlmodel import Session, create_engine

from gateways.rdb.config import DatabaseConfig

# グローバルなengineとsessionmakerを保持
_engine = None


def initialize_engine(config: DatabaseConfig | None = None) -> None:
    """データベースエンジンを初期化する.

    Args:
        config: Database configuration. If None, loads from environment.

    """
    global _engine  # noqa: PLW0603

    config = config or DatabaseConfig()

    _engine = create_engine(
        config.database_url,
        pool_size=config.pool_size,
        max_overflow=config.max_overflow,
        pool_recycle=config.pool_recycle,
        pool_pre_ping=True,
        pool_use_lifo=True,
        echo=config.database_echo,
    )


def get_sync_session() -> Generator[Session, None, None]:
    """FastAPIの依存性注入で使用するセッション取得関数.

    Yields:
        Session: データベースセッション

    Example:
        @router.get("/users/{user_id}")
        def get_user(user_id: int, session: Session = Depends(get_sync_session)):
            ...

    """
    if _engine is None:
        msg = "Engine is not initialized. Call initialize_engine() first."
        raise RuntimeError(msg)

    with Session(_engine) as session:
        yield session
