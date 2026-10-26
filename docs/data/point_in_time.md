# 時点管理台帳（point-in-time register）

「評価イベントより**前**の情報だけを入れる」ための版台帳。各レイヤーについて、
**いつの版をどこに置いたか**を 1 行 = 1 バージョンで記録する。ここで版が引けない
レイヤーは M1 に入れてはいけない（CLAUDE.md「漏洩禁止」）。

| レイヤー | 版 (vintage) | 取得日 | 保存先 / SHA | 使える評価イベント | 備考 |
|---|---|---|---|---|---|
| 地名コーパス（app） | app `sha` ごとに自動記録 | -- | `topsem.corpus.snapshot` | -- | YAML 来歴は app リポの git 履歴 |
| 洪水浸水想定区域 | 千曲川（国管理）: 現行版のみ（公表年不記載・取得 2026-09-20）。長野県管理: 更新 2026-04-15 明示 | 取得日 2026-09-20 | `data/covariates/flood_hazard/<river>/<muni>/` | 2019 評価には**使用禁止**（災害後に改訂された可能性）。2024 評価も改訂履歴不明なら感度扱い | 遡及版が不要な「現況のみ」共変量として strict 分割でのみ使用可 |
| DEM（基盤地図情報 5m） | TBD | -- | `data/covariates/dem5/<v>/` | -- | 更新頻度低・版の目印は公開日 |
| 河川線形（KSJ W05 / Geoshape） | 2006–2009 取得版（W05 vintage） | 2026-09-20 ダンプ | `data/covariates/ksj_rivers/`（10 桁 code JSON） | 2010 年以降のイベント | 河道改修は反映されない = 保守側誤差。M1 接続時に summary へ版記録 |
| 土地利用（MRLC） | 20xx 年版 | -- | `data/covariates/landuse/<v>/` | -- | 5 年刻みで版が明確 |
| 自然災害伝承碑 | 国土地理院 更新日 | -- | `data/covariates/monuments/<v>/` | -- | ダンプ日の記録必須 |
| 観測浸水ラベル | イベント発生日 + 資料公開日 | -- | `data/labels/<event>/` | そのイベントのみ | ラベル自体は「未来」なので特徴には厳禁 |

## v2（全球検定・Study 2）

評価事象は EMS activation の事象時刻（2012–2026）。各層の版と、事象との前後関係:

| レイヤー | 版 (vintage) | 取得日 | 保存先 / SHA | 事象との前後 | 備考 |
|---|---|---|---|---|---|
| 観測浸水ラベル（EMS DEL） | AOI ごとに最初の delineation の最新版 | 2026-10-25 | NAS `data/global/labels/<tag>/`（meta.json に zip sha256・cells_sha256） | ラベル自体 | 特徴には厳禁 |
| 地形 GLO-30 DSM | 2021 公開版（TanDEM-X 撮像 2011–2015） | AWS COG（実行時） | `copernicus-dem-30m` | 2012–2015 の事象と重なり得る | 事象で不変とみなす。protocol §7-8 で 2015 年以前を除外する感度 |
| HAND | ASF v1（GLO-30 2021 版から算出） | AWS COG | `glo-30-hand/v1/2021` | 同上 | 同上 |
| 河道 HydroRIVERS | v1.0（SRTM 2000 由来） | 2026-10-25 | NAS `data/global/hydrorivers/` | 全事象の前 | |
| 土地被覆 ESA CCI / C3S | **事象の前年**（2021 年以降の事象は 2020 年版）、300 m 年次 | Planetary Computer COG（実行時） | `labels/<tag>/landcover_meta.json` に年とタイル URL | 全事象の前 | WorldCover 2021 は 2012–2020 年の事象に事後情報になるので使わない |
| 地名 GeoNames | 2026-10-25 ダンプ | 2026-10-25 | NAS `data/global/geonames/allCountries.zip`（sha256 は .sha256） | 全事象の後 | 集落名は長く続くとみなす（限界として報告）。過去版のダンプは無い |
| 海岸線 Natural Earth | 10m coastline 5.0.0-pre9 | 2026-10-25 | sha256 bfa04cdb… | 不変とみなす | 種別未記載製品の内陸判定にのみ使用（特徴ではない） |

## 運用ルール

1. 共変量を読み込むコードは必ず `(layer, vintage)` を引数に取り、
   `summary.json` の `covariate_versions` に全レイヤーの版を書く
2. 新イベントを追加したら、そのイベント日より後の版しかないレイヤーを
   洗い出してこの台帳に印をつける
3. 「版が不明なファイル」を `data/covariates/` に置かない
