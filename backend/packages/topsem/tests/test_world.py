"""v2 全球検定の部品（topsem.world）の単体テスト。ネットワーク不要な合成データで。"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import shapely
from pyproj import Transformer

from topsem.world.ems import (
    CELL_M,
    build_aoi_labels,
    flood_class,
    laea,
    select_products,
)
from topsem.world.lexicon import CATEGORIES, Lexicon, cell_lexicon_features
from topsem.world.names import assign_names
from topsem.world.terrain import _block, _box_mean, _neigh_min

REPO = next(p for p in Path(__file__).resolve().parents if (p / "CLAUDE.md").is_file())
LEXICON = REPO / "docs" / "protocol" / "v2" / "lexicon_hydro_v1.yaml"


# ---- EMS ラベル -----------------------------------------------------------------


def test_flood_class_from_obj_desc() -> None:
    assert flood_class("Riverine flood") == "riverine"
    assert flood_class("Flash flood") == "flash"
    assert flood_class("Coastal flood") == "coastal"
    assert flood_class("Flooded area") == "other"
    assert flood_class("") == "other"


def test_select_products_one_per_aoi_prefers_overview_latest() -> None:
    prods = [
        {"aoi": "01", "kind": "PRODUCT", "scale": "DETAIL", "version": 3,
         "url": "x_DETAIL01_v3"},
        {"aoi": "01", "kind": "PRODUCT", "scale": "OVERVIEW", "version": 1, "url": "x_OV_v1"},
        {"aoi": "01", "kind": "PRODUCT", "scale": "OVERVIEW", "version": 2, "url": "x_OV_v2"},
        {"aoi": "01", "kind": "MONIT", "scale": "OVERVIEW", "version": 9, "url": "x_MON_v9"},
        {"aoi": "02", "kind": "PRODUCT", "version": 1, "url": "y_v1"},
    ]
    sel = select_products(prods)
    assert [p["aoi"] for p in sel] == ["01", "02"]
    assert sel[0]["url"] == "x_OV_v2"  # MONIT は除外、DETAIL より OVERVIEW、版は最新


def _square(lon0: float, lat0: float, dlon: float, dlat: float) -> dict:
    ring = [[lon0, lat0], [lon0 + dlon, lat0], [lon0 + dlon, lat0 + dlat],
            [lon0, lat0 + dlat], [lon0, lat0]]
    return {"type": "Polygon", "coordinates": [ring]}


def _fc(features: list[tuple[dict, dict]]) -> str:
    return json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "geometry": g, "properties": p} for g, p in features]})


def _make_product(tmp: Path, *, footprint: bool) -> Path:
    """経度 0.1° × 緯度 0.1° の AOI、西半分が河川浸水・北東角は Not Analysed。"""
    z = tmp / "EMSR999_AOI01_DEL_PRODUCT_v1.zip"
    aoi = _square(10.0, 45.0, 0.1, 0.1)
    flood = _square(10.0, 45.0, 0.05, 0.1)
    hole = _square(10.08, 45.08, 0.02, 0.02)
    fp = [(_square(9.9, 44.9, 0.3, 0.3), {"obj_type": "Image Footprint"})] if footprint else []
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("EMSR999_AOI01_DEL_PRODUCT_areaOfInterestA_v1.json",
                    _fc([(aoi, {"emsr_id": "EMSR999"})]))
        zf.writestr("EMSR999_AOI01_DEL_PRODUCT_observedEventA_v1.json", _fc([
            (flood, {"event_type": "5-Flood", "obj_desc": "Riverine flood",
                     "notation": "Flooded area"})]))
        zf.writestr("EMSR999_AOI01_DEL_PRODUCT_imageFootprintA_v1.json",
                    _fc([*fp, (hole, {"obj_type": "Not Analysed"})]))
    return z


@pytest.mark.parametrize("footprint", [True, False])
def test_build_aoi_labels_fractions(tmp_path: Path, *, footprint: bool) -> None:
    cols, meta = build_aoi_labels(_make_product(tmp_path, footprint=footprint))
    inside = cols["f_aoi"] >= 0.99
    assert inside.sum() > 100
    west = cols["lon"] < 10.045
    east = cols["lon"] > 10.055
    # 西半分の AOI 内セルはほぼ全て浸水、東半分は浸水なし
    assert np.mean(cols["f_riverine"][inside & west] > 0.9) > 0.9
    assert np.all(cols["f_riverine"][inside & east] < 0.05)
    assert np.all(cols["f_flash"] == 0.0)
    # Not Analysed の北東角は観測域から外れる（footprint が無い製品でも同じ）
    ne = (cols["lon"] > 10.085) & (cols["lat"] > 45.085) & inside
    assert ne.any()
    assert np.all(cols["f_observed"][ne] < 0.1)
    assert meta.get("footprint_implicit") is (not footprint)
    assert meta["grid"]["cell_m"] == CELL_M
    assert meta["obj_desc"] == {"Riverine flood": 1}


@pytest.mark.parametrize("typed", [True, False])
def test_build_aoi_labels_generation0(tmp_path: Path, *, typed: bool) -> None:
    """第 0 世代（2012–2013）: AreaOfInterest / Crisis_event_A / Hydrology_A_LAF。"""
    z = tmp_path / "EMSR044_04TORGAU_DELINEATION_OVERVIEW_v1_vector.zip"
    pre = "EMSR-044_04Torgau_20130606_DEU_FL_v1_"
    flood_props = {"SBTYPdes": "EM009 - Flood", "INTRP": "1"} if typed else {"INTRP": "998"}
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr(pre + "AreaOfInterest.json", _fc([(_square(10.0, 45.0, 0.1, 0.1), {})]))
        zf.writestr(pre + "Crisis_event_A.json",
                    _fc([(_square(10.0, 45.0, 0.05, 0.1), flood_props)]))
        zf.writestr(pre + "Hydrology_A_LAF.json",
                    _fc([(_square(10.06, 45.0, 0.02, 0.02), {"txt": "water"})]))
        # 島（陸地）が入る任意層は恒常水域として読まない
        zf.writestr(pre + "Hydrology_A_Optional_LAF.json",
                    _fc([(_square(10.06, 45.06, 0.03, 0.03), {"txt": "island"})]))
    cols, meta = build_aoi_labels(z)
    inside = cols["f_aoi"] >= 0.99
    west = inside & (cols["lon"] < 10.045)
    assert meta["schema_generation"] == 0
    assert meta["layers"]["hydro"].endswith("Hydrology_A_LAF.json")
    assert meta.get("footprint_missing") is True
    # 種別の記載が無いので "other"（§2.4 の内陸規則で河川とみなすかを決める）
    assert np.mean(cols["f_other"][west] > 0.9) > 0.9
    assert np.all(cols["f_riverine"] == 0.0)
    island = inside & (cols["lon"] > 10.065) & (cols["lat"] > 45.065)
    lake = inside & (cols["lon"] > 10.062) & (cols["lon"] < 10.078) & (cols["lat"] < 45.018)
    assert np.all(cols["f_perm_water"][island] == 0.0)
    assert np.mean(cols["f_perm_water"][lake] > 0.9) > 0.9


def test_laea_roundtrip_is_metric() -> None:
    crs = laea(10.0, 45.0)
    tr = Transformer.from_crs(4326, crs, always_xy=True)
    x0, y0 = tr.transform(10.0, 45.0)
    x1, _ = tr.transform(10.0 + 1 / 111_320 / np.cos(np.radians(45.0)), 45.0)
    assert abs(x0) < 1e-6
    assert abs(y0) < 1e-6
    assert abs((x1 - x0) - 1.0) < 0.01  # 経度 1 m 相当 ≈ 1 m


# ---- 辞書照合 -------------------------------------------------------------------


@pytest.fixture(scope="module")
def lex() -> Lexicon:
    return Lexicon(LEXICON)


@pytest.mark.parametrize(
    ("name", "cc", "expected"),
    [
        ("Frankfurt am Main", "DE", {"crossing"}),
        ("Donauworth", "DE", {"riverine"}),
        ("Oberhausen", "DE", {"high"}),
        ("Oxford", "GB", {"crossing"}),
        ("Vega de San Mateo", "ES", {"riverine"}),
        ("Slavonski Brod", "HR", {"crossing"}),
        ("Lunca Mures", "RO", {"riverine"}),
        ("Char Fasson", "BD", {"riverine"}),
        # 誤爆しやすい形を入れていないことの回帰テスト
        ("Eindhoven", "NL", set()),  # -ven（湿地）を suffix にしない
        ("Bexley", "GB", set()),  # -ey（島）を入れない
    ],
)
def test_lexicon_examples(lex: Lexicon, name: str, cc: str, expected: set[str]) -> None:
    hits = {k for k, v in lex.match(name, cc).items() if v}
    assert hits == expected


def test_lexicon_is_country_scoped(lex: Lexicon) -> None:
    # ドイツ語の -furt は英国の地名には適用しない（言語は地名の国で決まる）
    assert lex.match("Frankfurt", "DE")["crossing"] == 1
    assert lex.match("Frankfurt", "GB")["crossing"] == 0
    assert lex.coverage("DE")
    assert not lex.coverage("ZZ")


def test_lexicon_has_no_non_ascii_forms(lex: Lexicon) -> None:
    forms = [e.form for ents in lex.by_lang.values() for e in ents]
    assert all(f.isascii() for f in forms)


def test_cell_lexicon_features_weights_and_missing(lex: Lexicon) -> None:
    gids = np.array([[1, 2, -1], [2, 1, 1]])
    dists = np.array([[0.0, 1000.0, np.inf], [500.0, 0.0, 0.0]])
    lookup = {1: ("Oxford", "GB"), 2: ("Hill Top", "GB")}
    f = cell_lexicon_features(gids, dists, lookup, lex)
    assert f["lex_crossing_n1"].tolist() == [1.0, 0.0]
    assert f["lex_high_n1"].tolist() == [0.0, 1.0]
    # 2 行目: Oxford が距離 0 で 2 回 → 重み 1 + 1
    assert f["lex_crossing_k"][1] == pytest.approx(2.0)
    # 1 行目: Hill Top は 1 km → 重み 0.5、欠損（-1）は寄与しない
    assert f["lex_high_k"][0] == pytest.approx(0.5)
    assert set(CATEGORIES) <= {k.split("_")[1] for k in f if k.startswith("lex_")}


# ---- 地名割当 -------------------------------------------------------------------


def test_assign_names_nearest_and_counts() -> None:
    crs = laea(10.0, 45.0)
    tr = Transformer.from_crs(crs, 4326, always_xy=True)
    pts_xy = np.array([[0.0, 0.0], [500.0, 0.0], [5000.0, 0.0], [0.0, 400.0]])
    lon, lat = tr.transform(pts_xy[:, 0], pts_xy[:, 1])
    gn = pd.DataFrame({
        "gid": [11, 12, 13, 14], "lon": lon, "lat": lat,
        "fclass": ["P", "P", "P", "H"],  # H（水域名）は割り当てない
        "asciiname": ["a", "b", "c", "river"], "cc": ["IT"] * 4,
    })
    out = assign_names((np.array([100.0]), np.array([0.0])), crs, gn, (9.9, 44.9, 10.2, 45.1))
    assert out["P_gid1"].tolist() == [11]
    assert out["P_gid2"].tolist() == [12]
    assert out["P_d1"][0] == pytest.approx(100.0, abs=0.5)
    assert out["P_n1km"].tolist() == [2]
    assert out["P_n3km"].tolist() == [2]
    assert 14 not in out["PLTS_gid1"].tolist() + out["PLTS_gid2"].tolist()


# ---- 地形集約の部品 ----------------------------------------------------------------


def test_block_and_neighbour_min() -> None:
    a = np.arange(16, dtype=float).reshape(4, 4)
    b = _block(a, 2, 2, 2)
    assert b.shape == (2, 2, 4)
    assert sorted(b[0, 0].tolist()) == [0.0, 1.0, 4.0, 5.0]
    m = _neigh_min(np.array([[3.0, 1.0], [2.0, 5.0]]))
    assert m.tolist() == [[1.0, 1.0], [1.0, 1.0]]


def test_box_mean_ignores_nan() -> None:
    a = np.array([[1.0, np.nan, 3.0], [np.nan, 5.0, np.nan], [7.0, np.nan, 9.0]])
    m = _box_mean(a, 1)
    assert m[1, 1] == pytest.approx(5.0)  # (1+3+5+7+9)/5
    assert m[0, 0] == pytest.approx(3.0)  # (1+5)/2


def test_shapely_make_valid_available() -> None:
    # ems.read_layer が依存する（shapely 2 系の配列 API）
    bow = shapely.from_wkt("POLYGON((0 0, 1 1, 1 0, 0 1, 0 0))")
    assert shapely.make_valid(bow).is_valid


def test_landcover_year_is_always_before_the_flood() -> None:
    from topsem.world.landcover import CCI_GROUPS, group_lut, lc_year  # noqa: PLC0415

    assert lc_year("2013-06-02T00:00:00") == 2012
    assert lc_year("2021-07-14T10:00:00") == 2020
    assert lc_year("2026-10-19T03:10:00") == 2020  # 最新版で頭打ち（なお洪水より前）
    lut = group_lut()
    assert lut[0] == 0  # 欠測
    assert lut[190] == list(CCI_GROUPS).index("built") + 1
    assert lut[180] == list(CCI_GROUPS).index("wetland") + 1
    codes = [c for v in CCI_GROUPS.values() for c in v]
    assert len(codes) == len(set(codes))  # 1 区分は 1 群にだけ属する
