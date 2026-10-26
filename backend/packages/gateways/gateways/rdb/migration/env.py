import sys
from logging.config import fileConfig
from pathlib import Path
from typing import Literal

from alembic import context
from sqlalchemy import engine_from_config, pool
from sqlmodel import SQLModel

# Add the site-packages directory to sys.path to import gateways package
# This is necessary when running Alembic in Docker or other isolated environments
site_packages = Path("/usr/local/lib/python3.12/site-packages")
if site_packages.exists() and str(site_packages) not in sys.path:
    sys.path.insert(0, str(site_packages))

from gateways.rdb.config import DatabaseConfig  # noqa: E402
from gateways.rdb.models.schema import *  # noqa: F403, E402

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

# Interpret the config file for Python logging.
# This line sets up loggers basically.
if config.config_file_name is not None:
    fileConfig(config.config_file_name)


config.set_main_option("sqlalchemy.url", DatabaseConfig().database_url)
target_metadata = SQLModel.metadata


# sqlmodelを使ったautostringを生成すると、以下のエラーに引っ掛かるため
# gateways/rdb/migration/versions/d772bf336b6a_init.py:27:
# error: Module has no attribute "sql"  [attr-defined]
def render_item(_type_, obj, _autogen_context) -> str | Literal[False]:  # noqa: ARG001, ANN001
    """Render SQLModel types as SQLAlchemy types."""
    from sqlmodel.sql.sqltypes import AutoString  # noqa: PLC0415

    if isinstance(obj, AutoString):
        # Convert SQLModel AutoString to SQLAlchemy String
        return f"sa.String(length={obj.length})"
    return False


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    url = config.get_main_option("sqlalchemy.url")
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        render_item=render_item,
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """
    connectable = engine_from_config(
        config.get_section(config.config_ini_section, {}),
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )

    with connectable.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
            render_item=render_item,
        )

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
