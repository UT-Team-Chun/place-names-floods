"""Base model for all SQLModel tables."""

from datetime import UTC, datetime

from sqlalchemy import Column
from sqlalchemy.types import TIMESTAMP
from sqlmodel import Field


def utc_now() -> datetime:
    """Get current UTC time with timezone info."""
    return datetime.now(UTC)


class TimestampMixin:
    """Mixin to add created_at and updated_at timestamps."""

    created_at: datetime = Field(
        sa_column=Column(
            TIMESTAMP(timezone=True),
            nullable=False,
        ),
        default_factory=utc_now,
        description="作成時間",
    )
    updated_at: datetime = Field(
        sa_column=Column(
            TIMESTAMP(timezone=True),
            nullable=False,
        ),
        default_factory=utc_now,
        description="更新時間",
    )
