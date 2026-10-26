import pandas as pd

from topsem.corpus.population import NameIndex, infer_name_column, population_stats, tag_rows


def _index() -> NameIndex:
    return NameIndex(
        curated={"蛇抜け": "蛇抜け", "ジャヌケ": "蛇抜け"},
        hazard={"蛇抜": "januke", "沢": "sawa", "阿": "awa"},
        flood=frozenset({"awa"}),
    )


def test_tag_rows_flags_hazard_and_flood() -> None:
    df = pd.DataFrame({"name": ["蛇抜沢", "阿原", "新宿", None]})
    tagged = tag_rows(df, "name", _index())
    assert tagged["is_flood_word"].tolist() == [False, True, False, False]
    assert tagged["hazard_elements"].tolist() == ["januke,sawa", "awa", "", ""]
    assert tagged["is_curated"].tolist() == [False, False, False, False]


def test_population_stats_shapes() -> None:
    df = pd.DataFrame({"name": ["蛇抜", "阿原"]})
    stats = population_stats(tag_rows(df, "name", _index()))
    assert stats["n_rows"] == 2
    assert stats["frac_flood_word"] == 0.5
    assert stats["frac_any_hazard_word"] == 1.0


def test_infer_name_column() -> None:
    assert infer_name_column(pd.DataFrame({"pref": [1], "name": ["a"]})) == "name"
    assert infer_name_column(pd.DataFrame({"大字": ["a"]})) == "大字"
