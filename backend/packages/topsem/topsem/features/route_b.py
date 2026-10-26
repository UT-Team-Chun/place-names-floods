"""route_b: 地名文字列の埋め込み。

パイロット実装 = 文字 n-gram ハッシュベクトル（次元固定 512、L2 正規化）。
汎用事前学習埋め込み（計画 §5.2 の本命）は版・学習期間の記録と暗記チェックが
揃い次第そちらに置き換える。置き換えたら protocol §7 の該当 TBD を解消する。
"""

from __future__ import annotations

import hashlib
import math

DIM = 64


def _ngrams(s: str, n_lo: int = 1, n_hi: int = 3) -> list[str]:
    s = f"^{s}$"
    out: list[str] = []
    for n in range(n_lo, n_hi + 1):
        out.extend(s[i : i + n] for i in range(len(s) - n + 1))
    return out


def embed_name(name: str, dim: int = DIM) -> list[float]:
    """文字 n-gram をハッシュで dim 次元に散らし、L2 正規化。"""
    v = [0.0] * dim
    for g in _ngrams(str(name)):
        h = int.from_bytes(hashlib.blake2b(g.encode("utf-8"), digest_size=4).digest(), "big")
        v[h % dim] += 1.0
    norm = math.sqrt(sum(x * x for x in v))
    if norm == 0.0:
        return v
    return [x / norm for x in v]
