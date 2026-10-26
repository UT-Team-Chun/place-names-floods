import math

from topsem.features.route_a import ElementIndex, name_features
from topsem.features.route_b import DIM, embed_name

ELEMS = (
    {
        "id": "januke",
        "surfaces": ["蛇抜", "蛇退"],
        "readings": ["ジャヌケ"],
        "hazard_types": ["debris_flow"],
        "specificity": 0.95,
    },
    {
        "id": "awa",
        "surfaces": ["阿"],
        "readings": ["アワラ"],
        "hazard_types": ["flood", "inundation"],
        "specificity": 0.7,
        "substitution_of": "原",
    },
)


def test_route_a_surface_and_reading_hits() -> None:
    idx = ElementIndex.from_elements(ELEMS)
    f = name_features("蛇抜沢", None, idx)
    assert f["a_n_elements"] == 1.0
    assert f["a_flood_hit"] == 0.0
    assert f["a_spec_max"] == 0.95
    f2 = name_features("原", "アワラ", idx)  # 読みヒット（換字元の表記）
    assert f2["a_flood_hit"] == 1.0
    assert f2["a_n_substituted"] == 1.0
    f3 = name_features("新宿", None, idx)
    assert f3["a_any_hit"] == 0.0


def test_route_b_normalized_and_deterministic() -> None:
    v = embed_name("蛇落地悪谷")
    assert len(v) == DIM
    assert math.isclose(math.sqrt(sum(x * x for x in v)), 1.0, rel_tol=1e-6)
    assert embed_name("蛇落地悪谷") == v


def test_route_b_similarity_ordering() -> None:
    a = embed_name("蛇落地悪谷")
    b = embed_name("蛇落地")
    c = embed_name("新宿駅東口")
    dot = lambda x, y: sum(p * q for p, q in zip(x, y, strict=True))  # noqa: E731
    assert dot(a, b) > dot(a, c)
