# B-13 初見の査読（第 9 回）に応じた組み直しと事後の解析

状態: **計画（2026-10-29T12:19Z、commit c096169、下の B・C の量を計算する前に固定）**。確証の判定（equivalent-null）は変えない。
プロトコルの判定はそのまま記録として報告し、その誤り率の問題を直す事後の検定（B1）を並べる。全て事後（探索）。
対象プロトコル: [`preregistration_v2.md`](../preregistration_v2.md)（FROZEN 2026-10-25T20:58:37Z、commit e625924）。

書いた時点で見ていた値: 確証の解析と B-2〜B-12 の結果。ただし B-9 の C1〜C4（クラスタで計算中）、B-12 の B11（名前の
係数だけに罰則、手元で計算中）と C1（1 km の注入、クラスタで計算中）は見ていない。下の B・C の量はどれも計算していない。

## 背景

初見の査読 第 9 回（NHESS 版、本文 33 頁・補足 32 頁、開発リポジトリ 40b0a3d の版、公開リポジトリ 5622df7。報告は
[reviews/](../../research/reviews/) の `2026-10-29_paper2_blind9_*.md`）は、編集者が「審査の前に差し戻し」、査読者 3 人が
major revision。4 人に共通の要点:

1. **判定の規則は誤り率を制御しない**: 固定したモデルの区間を点の物差しで割る規則は、真の R が余白のとき 32〜33% の
   標本で同等と判定する（B-12 の B1）。当て直しで R の標準偏差は約 2.3 倍。分子と分母の共分散を運び、二段目を当て
   直し、上位の単位で再抽出する有効な同等性の検定（比の余白に対する Fieller 型の TOST）を解釈の主に（査読者 2 の M1・
   M4、編集者の M3）
2. **数と結びつきの誤り**: 討論の OR 1.83 は当て直した場所グループの区間（上限 0.013）のもので、二方向の区間（上限
   0.036）なら名目で約 2.1（推定量の尺度で約 2.6）。p0 の 18% は辞書を作った符号化者 1 の値で、符号化者 2 なら約 23%。
   希釈の例は符号化者 1 の率（査読者 2 の M6、査読者 3 の M5、編集者の M6）
3. **正の向きの結果の報告不足**: H4 の全体の増分（positive-small）、GFM の地域内の R 0.011（placebo を上回る）が要旨・
   結論に無い（編集者の M5、査読者 2 の M9）
4. **対照の漏れ**: 地表水（1984–2021）と再起の記録は当の事象と後の事象を含む。事象の前に分かる対照で希釈を示す
   （査読者 1 の M1・M7、編集者の M1）
5. **集落の近くの推定**: 識別しない推定に本文の多くを割き、上限なしの推定は選んだもの、少ないクラスタの区間、裏付けの
   規則、ラベルの弱さを結論に（査読者 2 の M7、査読者 1 の M2・M3・M6、査読者 3 の M7）
6. **辞書**: 同形異義（-ach と -bach、-ore と Store、mont-、Ried、kacha、taung）、国に付くブロック（イタリア全土を
   ドイツ語のブロックで、ベルギーの -eau を Au）、言語で違う照合の規則、相対の修飾語、照合が符号化者の一致だけ
   （査読者 3 の M2〜M6、編集者の M2）
7. **長さと造語**: water-side・not-meaning・M の記号が多い。題が長い（34 語）。本文を 4 部に（確証の結果と全ての
   不確かさ、指定した少数の確認、設計の天井と希釈、読みの妥当性）
8. **確かめられないもの**: 凍結の記録、未確認の英訳、データの寄託、生成 AI の方針（著者の判断）

## A. 書き直し（計算なし。論文と補足で直す）

1. **題**（約 20 語、事前固定は第三者が確かめられないので protocol-based に）: "Water-related terms in settlement names,
   assigned to the nearest grid cells, did not improve within-area predictions of mapped river flood extent: a
   protocol-based test in 500 Copernicus EMS areas"
2. **要旨**: 設計の限界（集落の近くだけの効果は余白に届かない、集落を通じた割り当ては事象の前に分かる情報も薄める）を
   数の前に。R・活性化・基礎の利得を一文で定義。結果のラベル（一回の写像、浸水の多い地域を選んだ）を書く。プロトコルの
   判定は記録とし、B1 の妥当な検定の結果を書く。二方向の区間は余白の内に入らず、上限としての読みは結果の後に決めたと
   書く。H4 の全体の増分と GFM の地域内の R を書く。OR の換算は「名の最寄りの全ての地面に一様に広がる効果」と明記するか
   外す。「Brier score がわずかに増えた」と書く
3. **第 3 節**: プロトコルの判定（記録）と B1 の妥当な検定を並べ、活性化をまとめた分割を両方のコードで主の隣に。H4 の
   全体の増分（positive-small）と柔軟な H1（正）を本文と表 C1 に。「H1 は外れた」は入れ子の形の予測についてだと書く
4. **第 4 節**: 主の不確かさを B1 にまとめ、二方向・当て直しの区間は構成要素として短く。分母の代替（Y3 と ΔBS、Brier
   skill の増分）、南アジアの影響（B6）、基準のモデルの漏れ（B5）。緩めた二段目は版ごとに判定を書く（大陸の交互作用は
   negative、JRC の文の矛盾を直す）
5. **第 5 節**: 漏れのない対照（C1）を主に、当の事象を含む対照は「設計が拾えるもの」として残す。上限は「対応する」と
   書き、どの区間の上限かを結ぶ（B2）。無効果の基準ごとの範囲を同じ重みで。モンテカルロの範囲は信頼区間でないと書く。
   希釈の例は両方の符号化者の率で（B3）。「過去の EMS の洪水が記録のある AOI で拾われた」は数と判定で書き直す
6. **第 6 節**: (a) 読み: 表 5 は「言語モデルの符号との一致」と書く。符号化者の取り違えを直す。辞書の系統的な誤りの
   表（B10）を補足に、本文は要点。N-read は場所の知識を含む読みと書く。(b) 集落の近く: 一段落に。計画した GFM の推定を
   計画した区間で、上限なしの EMS の推定は「地図のラベルで正の関連は見えない」探索の結果として、識別しないことを添えて。
   類の話は補足に。裏付けの規則を外す。ラベルの弱さは確かめていない仮説として書く。GFM の推定の扱いを両方の向きで同じに
7. **討論・結論**: practical → simple。NHESS の読者が得るもの（集落を通じた割り当ては事象の前に分かる情報も薄める一般の
   教訓、集落の規模の検定に要るもの: GFM は市街地を評価しないので干渉性・二重反射の方法、Global Flood Database）
8. **用語**: water-side → water-related（渡河・谷・「下」の修飾語を含むと定義に書く）。not-meaning は「placebo と区別
   できない」と言い換え。表 1 を用語集に（M0〜M3s、T0〜T4、Y0〜Y4）。activation・location group・ground を初出で定義
9. **補足・付録の直し**: 年表に改訂 12 と第 8・9 回の批評。表 S5（計画の状態）の節の番号を直し、改訂 12・13 の行を
   足す。表 S6 の題の矛盾と欠けた注 b。表 S17 の B2 の活性化の行に国を。B4 の行の文字の崩れ。「Sentinel-1 以外の AOI」
   の数（68・111・93・108）の数え方を題に。表 3 の題と見出し。図 1b の「T4 = M1 + HAND」。図 4 の印と色の説明。図 S6 の
   行の名を平易に。表に無い数（距離の中央値 1.5 km、高地の OR、校正の切片と傾き、集落の近くの R 0.060、H4 の値、相関
   0.04、豪州の 9.7%、氾濫原の 1.45、較正の利得、付録 C の量）を表に。文献（ESA WorldCover、Natural Earth、GSW の版、
   空間の交差検証と浸水の起こりやすさの地図の文献、Owen と Nadeau–Bengio を本文で、補足の未引用の文献）。sha256 を
   全桁で。日本の先行の検定の概要を補足に
10. **数と言い回しの直し**: 33% → 32〜33%、「more than 98%」→「at least 98%」、片側 p の言い回し、0.017 の文の矛盾、
    土地被覆は 2022 年以降は 2020 年の図、地名辞典の取得は全ての事象より後、種別の 154 の内訳、付録 C の「3 分の 2」の
    尺度、641 の推定は改訂 2〜10、アフリカの 0.5% と表 5 の 1%、表 2 の柔軟な T0・T1 の上限が同じ値、表 S13 の JRC の
    90%・95% が同じ値（どれも既存の値の確認と直し）
11. **再現の限界**: arm64 と amd64 で模擬が最大 1.4 MC SE 違う、GFM のラベルの作り直しが 192 のうち 189 で一致、を
    本文か補足に

## B. 既存の成果物からの計算（手元、スクリプト v70。乱数の種 20261102）

1. **妥当な同等性の検定**（Fieller 型 TOST）: L± = ΔBS ∓ 0.05·Δ_base（ΔBS = BS(M1c) − BS(M3n)、Δ_base = BS(M0) −
   BS(M1)）を同じ再標本で計算し、L+ の片側 95% 上限 < 0 かつ L− の片側 95% 下限 > 0 なら同等（R の 90% の
   Fieller 区間が ±0.05 の内）。再抽出: (a) 場所グループ、モデル固定（4,000 回）; (b) 場所グループ、二段目を各回で
   当て直す（500 回、B-11 の B7 と同じ当て直し）; (c) 活性化を単位（国の中の活性化を二段で再抽出、二段目を当て直す、
   500 回）; (d) (c) に、基礎のモデルを当て直す分散を 13 通りの fold の割り当て（v47）から Nadeau–Bengio の補正で足す
   正規近似。主比較（確証の実行と直した実行）、H4 の入れ子と全体の増分、国の分割。動作特性: B-12 の B1 の二重の再抽出
   （モデル固定）で、真の R を ±0.05 と 0 に置いたときに L 規則が同等と判定する割合（一様のずらしと比例のずらし）
2. **上限の換算**: 各区間の上限（固定 0.003、当て直し、二方向 0.036、B1 の各上限）に対応する注入の OR（名目と推定量の
   尺度）を v64 の曲線で。無効果の 4 つの基準ごと。モンテカルロの範囲つき
3. **符号化者 2 の p0**: 読まなかった名の層の、符号化者 2 の符号で水に関わる由来の割合（区間つき、直接数える）。希釈の
   表と誤読の補正（ten_misclass の確率的な偏りの分析）を両方の符号化者の率で
4. **集落の近くの規模の層**: v55 の規模の三分位（上限の表）を本文の数に（新しい計算なし、表を作るだけ）
5. **基準のモデルの漏れ**: 主分割（直したコードの二段目）で、文字モデルの得点の AOI 平均を除いた M1c と M3n の R
6. **影響**: 主の R（確証の実行）と B1 の (b) について、影響の大きい場所グループ 10 を一つずつ除いた値と、影響に基づく
   実効のクラスタ数（Σ影響² に対する）
7. **分母の代替と絶対の量**: 主比較・H4 の入れ子と全体の増分で、Y3（定数の予測からの利得）での R、ΔBS（Brier の単位）、
   気候値に対する Brier skill の増分
8. **名の係数の同時検定**: 二段目の M3n を全てのセルで当て、名の 14 の係数が全て 0 の同時 Wald 検定（場所グループで
   クラスタ頑健）。B-7 の追記 3 の「held-out の係数の検定」（v55 の heldout_coef_test・wald_tests）を探して報告
9. **集落の近くの点検**: (a) 座標が整数の分でない集落だけで、上限なしの EMS の OR（pairs）; (b) 活性化を一つずつ除いた
   上限なしの EMS の OR（点）と、上位 2 つ（EMSR052・EMSR130）を除いた pairs 区間; (c) 世代 0・1 の産物を除いた OR;
   (d) 計画した GFM の推定と上限なしの EMS の推定を、名を AOI の中で入れ替えた置換（200 回）の分布と比べた p; (e) EMS と
   GFM の 2×2 の一致を距離帯と市街地の割合で; (f) 水系側 × 500 m × センサー（光学とレーダー）の交互作用の検定（上限の
   表）; (g) GeoNames の地物の符号が PPLX（集落の一部）・PPLH・PPLQ・PPLW（歴史上・廃れた・壊れた集落）の集落を除いた
   上限なしの EMS の OR と計画した GFM の OR
10. **辞書の点検**: 同形異義と言語の割り当ての数え上げ（-ach で -bach でない名、-ore・-ora と Store、mont- の修道院
    の名、Ried、kacha、taung、ベルギーの -eau の Au、イタリアでドイツ語のブロックが読んだ名、イランの ab の全体の語と
    接尾の数え方、ala-/yli-）。ブロックと国の表、照合の位置（全体・頭・尾）の言語ごとの表（補足に）。手元の GeoNames
    の name 列（Unicode）を、ä→a・ö→o・ü→u・å→a・ø→o・æ→ae・ß→ss と結合の印を外す規則で読んだときに、asciiname と
    比べて得る・失う照合の数を言語ごとに
11. **N-read**: 区分と unknown の割合をブロックと国ごとに
12. **二つの符号化者**: 一致（κ）と、GeoNames の ID・国・層・二人の符号を並べたファイル（公開、照合の標本とは別）
13. **AOI で釣り合わせた指標**: AOI ごとの ΔAUC・ΔPR-AUC の中央値（M1c → M3n）、AOI で重みを揃えた Brier skill の増分
14. **GFM の当て直しの崩れの診断**（手元の成果物がある範囲で）: AOI ごとの GFM と EMS の浸水割合、GFM での AOI ごとの
    基礎の利得
15. **相対の修飾語の対**（査読者 3 の M7。補足）: 同じ国の集落で、名から C2 (d) の修飾語を除いた残り（修飾語が頭の
    位置で付く言語では、頭を除いた残り）が同じで 3 文字以上、一方が下の側、他方が上の側の修飾語を持ち、互いに 10 km
    以内の対。両方の座標が整数の分でないものに限る。対ごとに (a) 地点の HAND の差（下 − 上）、(b) 500 m 以内の評価
    できるセルの EMS の浸水割合の差（両方にセルがあるとき）、(c) 計画した GFM の領域で同じ量。対の数、差の中央値、
    下の側が低い（浸水が多い）対の割合と符号検定の p、対を場所グループで再抽出した 90% 区間

## C. クラスタの計算

1. **漏れのない対照**（二段目に入れる。基礎のモデルには入れない）: (a) 地表水の発生頻度（GSW 1984–2021）を、2022 年
   以降の事象の AOI だけで（セルに直接と、集落の 500 m 以内で要約して集落を通じて）; (b) EMS の他の事象の再起の記録を、
   その事象より前の事象だけから（セルに直接と集落を通じて）; (c) JRC の全球の河川洪水ハザード図（100 年確率の浸水の
   有無）を集落の 500 m 以内で要約して集落を通じて。主分割、直したコード。判定表の判定も
2. **誤りを直した辞書**（派生のファイル `docs/verification/v70_lexicon_corrected.yaml` を、凍結した辞書から
   `backend/verification/r_okauchi/v70_lexicon_corrected.py` で機械的に作る。凍結したファイルは変えない）。直すのは
   査読者 3 が名指した項目と、分けるよう求めた三つの類だけ。それ以外の項目の当否は地名学者の照合（D5）に残す。
   照合の規則（英字以外で分ける、接辞の後に 2 文字以上）と六つの区分・14 の名の特徴は変えない（外した形は読まない）。
   形は全ての位置（語・頭・尾・中）で外す。ただし (b) の -ach と mont- は書いた位置だけ。
   - (a) v57 の掃除（2 文字以下の頭・尾の項目と、英語の尾の -burn・-bourne・-borne・-pool・-well・-wells を外す）
   - (b) 同形異義（査読者 3 の M3 (iii)・(ii)）: de の尾の ach（尾の bach は残すので -ach は -bach の形だけ読む）、
     de の ried、fr の頭の mont（語の mont・monts は残す）、no の ore・ora、fa の kacha、my の taung、pt と es の
     monte、en の bank・banks
   - (c) 地名学の上で洪水の意味を持たないと名指された川沿いの項目（M3 (ii)）: sv の ang・angen・angar・vall、
     no の eng・enga・engene・voll・vollen、da の eng・enge（no の eng と同じ語）、no の oy・oya、fr の ile・iles・
     isle・ilot、sv と no の sand・sanden。島の語は他の言語では残す（川の中の島にも使い、扱う範囲が川の洪水の地域なので）
   - (d) 相対の修飾語（上・下の対の一方を示す語。M3 (ii)、M7）。下の側: en low・lower・nether、de unter・nieder、
     nl laag・neder・beneden、fr bas・basse・sous、it basso・bassa・bassi・sotto・inferiore、es bajo・baja・bajos・
     bajas・abajo、pt baixa・baixo、el kato、sh donji・donja・donje、sl spodnji・spodnja・spodnje、bg dolna・dolni・
     dolno・dolen、sq poshte・poshtme、ro jos、hu also、cs dolni・dolny・dolne・dolna、pl dolny・dolna・dolne、
     uk nyzhnii・nyzhnia・nyzhnie・nizhniy・nizhnyaya・nizhneye・nizhne・nizhni、lt zemutinis・zemutine・zemieji、
     et alam、fi ala、sv・no・da nedre、tr asagi、ar tahta・tahtani・tahtaniya・sufla・sufli、fa nichla・nichli・
     nichle・payan・payin・paeen・pain・payeen、hi nichla、id bawah・hilir・ilir、tl ibaba・baba、sw chini、
     mg ambany、ka kvemo。上の側: en high・upper、de ober、nl hoog・boven、fr haut・haute・hauts、it alto・alta・
     alti・sopra・superiore、es alto・alta・altos・altas・arriba、pt alto・cima、el ano・epano、sh gornji・gornja・
     gornje、sl gorenji・gorenja・gorenje・zgornji・zgornja・zgornje・gornji、bg gorna・gorni・gorno・goren、
     sq siper・siperme、ro sus、hu felso、cs horni・horny・horne・horna、pl gorny・gorna・gorne、uk verkhnii・
     verkhnia・verkhnie・verkhniy・verkhnyaya・verkhneye・verkhne、lt aukstutinis・aukstoji、et ulem、fi yla、
     sv・no・da ovre、tr yukari、ar fawqani・foqani・fouqani・fawqa・oulia・ulya・ulia、fa bala、hi upar、
     id atas・hulu・ulu、tl itaas、sw juu、mg ambony、ka zemo。高さの形容詞（de hoch・tief、cs・pl・uk の vysok-
     など）と谷・低地の名詞は残す
   - (e) 泉・井戸（M3 (ii)、水系側から外すよう求めた）: en well・wells・spring・springs・soak、de born・brunn・
     brunnen・quell・quelle、fr fontaine・fontaines・source・sources・puits、it fonte・fontana・fontane・sorgente・
     pozzo・pozzi、es fuente・fuentes・pozo・pozos・manantial・ojo、pt fonte・fontes・poco・nascente・cacimba、
     el pigi・pigadi・pigadia・vrysi・vrisi・vryses・vrises・krini・kefalovryso・vrys・vris、sh vrelo・vrela・izvor・
     izvori・bunar・studenac・vrutak、sl studenec・izvir・studenci、bg izvor・izvori・kladenets・kladenci・bunar・
     cheshma、sq burim・burimi・pus・puse、ro izvor・izvorul・izvoare・fantana・fantanele、hu kut・forras、cs studnice・
     studna・studnicka・pramen・prameny、pl zrodlo・zrodla・studnia・studzianka、uk krynytsia・krinitsa・kolodez・
     kolodiazi・dzherelo・istok、lt saltinis・saltiniai、lv avots・avoti、et allikas、fi lahde、sv kalla・kallan、
     no・da kilde、tr pinar・kuyu・bulak・kaynak、ar ain・ayn・aioun・ayoun・bir・biar・abar・hassi、fa chashma・
     chashmeh・chah・chahi、vi gieng、id sumur、tl bukal、sw chemchemi・kisima、bantu chitsime、af fontein・put、
     mg loharano、qu puquio・pukyu、nah ameyal、ka tskaro・tskarosi、kk bulak・bulaq、so ceel、am minch・minchi
   - (f) 早瀬・滝（M3 (ii)）: en falls、fi koski、sv fors・forsen、no foss・fossen
   - (g) 名指された欠け（M3 (iv)）を足す: fr の語 sur（riverine）、cs と pl の語 nad（riverine）、de の語と尾の hamm・
     damm・deich（riverine）、es の語 rambla（water）、hu の語と尾の fok（water）、km の語 trapeang（water）
   - (h) 国とブロック（M4 (i)）: de のブロックをベルギーとイタリアから外す（ドイツ・オーストリア・スイス・
     リヒテンシュタイン・ルクセンブルクだけに）

   この辞書で、主比較（world_v2_b13lexfix、腕 M0、積み重ね M1c・M3n・M3np）と、集落の近くの OR（上限ありとなし、
   EMS と計画した GFM）。凍結した辞書・v57 の辞書との照合の数の差を言語ごとに。規則はこの計画のとおりで、結果を見て
   変えない
3. **早期終了を場所グループで**: M0・M1 の反復の回数を、訓練側の場所グループを分けた内側の検証で選び（無作為の 1 割の
   早期終了をやめて）、主分割を当て直す。R、分母、定数の予測との比較
4. **1 km の注入**（B-12 の C1）: 走行中の結果を入れる

## D. 著者の判断・別の研究（この改訂では行わない）

1. 生成 AI の方針（B-12 の D1）。記載の位置と、どの部分を言語モデルが下書きし著者が何を確かめたか
2. 英訳の確認（プロトコルと全ての改訂）
3. 確かめられる凍結の記録（開発の履歴とクラスタの記録を第三者に）。開発リポジトリは共有しないという指示があるので
   著者が決める
4. 派生データの寄託（Zenodo、査読者が使える形で）
5. 地名学者の照合（300 名と手引き）。結果が出たら、照合した読みで主比較と集落の近くの OR をやり直す（規則は照合の前に
   固定する）
6. ラベルと基礎の作り直し（詳細図、監視図による最大の範囲、画像の時刻、世代、GFM の場面ごとの作り直しと助言の印、
   FABDEM と流路に沿う HAND、JRC の欧州の図、GFPLAIN250m、防御の代理）は集落の規模の研究で。この論文では、それらに
   依る主張（集落の近くでラベルが弱い、など）を確かめていない仮説として書く
7. 判型（Brief communication）または誌の選択
8. 集落の集合と読みの作り直し（PPLX・PPLH・PPLQ・PPLW と人名の記録を除いて最寄りの割り当てをやり直す、Unicode の名と
   別名で読む、名ごとに言語と層を判定する、少数の名の登録簿）は集落の規模の研究で

## 報告の規則

- 全て事後（探索）で、確証の判定は変えない。プロトコルの判定は記録として報告し、B1 の妥当な検定を解釈の主に置く
  （付録 A の逸脱に加える）
- 結果の向きにかかわらず報告する。B1 で同等が言えない形があれば、そのとおり書く
- 計算した後に規則・単位・回数を変えない。実装の誤りは追記で直す

## 追記 1（実装の記録。2026-10-29T12:37Z、commit a7b0f33、B・C の量を計算する前）

1. **C1(a)**: `--min-event-year 2022`（2022 年以降の事象の AOI だけの表で、分割・基礎のモデル・二段目を作り直す）。
   集落を通じた地表水は既存の `settlement_gsw.parquet`（集落の点から 500 m 以内の occurrence・seasonality の平均。
   B-6 の 3(b)）で M1cG → M3oG、恒常水域を除いた平均で M1cGn → M3oGn。セルに直接は既存の `cell_gsw.parquet`
   （B-8 の C4）で M1cW → M3oW。名の主比較（M1c → M3n、placebo M3np）も同じ表で。実行 world_v2_b13gsw22
2. **C1(b)**: `--recurrence --recurrence-prior`。他の事象は、範囲が重なり activation が違い、その事象より 30 日を
   超えて前のもの（B-7 の 4 は前後どちらも）。集落を通じて M1cE → M3oE、セルに直接 M1cEc → M3oEc。実行
   world_v2_b13recprior
3. **C1(c)**: `--jrc-set`。集落の点から 500 m 以内のセル（JRC のタイルの画素が半分以上あるもの）の 100 年の浸水の
   割合の平均、その中に浸水するセルがあるか（計画の「有無」）、値があるかの印の三つを、最寄りの集落から割り当てる
   （M1cJ → M3oJ）。比べる相手として、同じ実行でセルに直接の 100 年の浸水の割合（M1cJc → M3oJc）も出す（計画に無い
   比べる相手を、結果の前に足す）。実行 world_v2_b13jrcset
4. **C2**: 実行 world_v2_b13lexfix（`--lexicon docs/verification/v70_lexicon_corrected.yaml`、sha256 c1ed4f18300a…）。
   派生の辞書は 2,538 → 2,148 項目（項目は形と位置の組。外したのは v57 の掃除 15、同形異義 15、川沿いの項目 37、
   相対の修飾語 181、泉・井戸 145、早瀬・滝 10。足したのは 13）。記録は `docs/verification/v70_lexicon_corrected_record.json`
5. **C3**: `--grouped-es`。基礎のモデル（M0 と M1。内側の交差適合と訓練側の全体の当て直しの全て）の各当て直しで、
   訓練側の場所グループの 1 割（四捨五入、1 以上）を当て直しの種で選んで内側の検証にし、残りで早期終了なしに
   200 回まで当て、検証の log loss が最小の回数で訓練側の全体を当て直す（従来の既定は無作為の 1 割で、10 回改善
   しなければ止める）。実行 world_v2_b13gres
6. C の実行は全て主分割（10 fold、1 Job 1 fold、4 CPU・14 Gi）、直したコード、既定の 5 つの種
7. **B1**: (a) は確証の実行（world_v2）と直した実行（world_v2_lexfix）の held-out の予測から、場所グループの再抽出
   4,000 回。(b)(c) は B-11 の B7 と同じく、直したコードの二段目の設計（world_v2_b8stage2_local）で当て直す（確証の
   実行の二段目の設計は保存していない）。(c) の二段の再抽出は、国を等確率で国の数だけ引き、引いた国ごとにその国の
   activation をその数だけ引き直す。(d) は v47 の 13 通りの fold の割り当て（world_v2_lexfix と world_v2_fs1b・fs2b・
   fs3〜fs12）の held-out の予測から L± を計算し、その分散を Nadeau–Bengio で補正して (c) の分散に足す（正規近似、
   片側 95%）。L± はセルあたりの Brier の尺度（重みの和で割る）。Fieller の区間も (a)〜(c) の再標本の分散と共分散から
   正規近似で出す
8. **B15**: 地点の HAND は v37 と同じ（集落の点から 177 m 以内のセルの hand_min の平均）、座標が整数の分かも v37 の
   規則。集落は解析の表（v41 の表、上限あり）でいずれかのセルの最寄りの集落。対は、同じ国・同じ基の名の下の側と上の
   側の集落を距離の近い順に一つずつとる（各集落は一回だけ）。浸水の割合は同じ表の 500 m 以内のセルで、GFM は v41 の
   GFM のラベルのあるセルで。区間は下の側の集落の場所グループで再抽出（2,000 回）

## 追記 2（C1(a) の実装の直し。2026-10-29T12:40Z、commit 39f2082、C の量を計算する前）

追記 1 の 1 の `--min-event-year 2022` は、AOI を絞ってから場所グループの上限の中でセルを抜き直すので、主の表に無い
セルが入り、セルごとの地表水の表（`cell_gsw.parquet`、主の表の 1,507,486 セルだけ）に値が無いセルが出る。そこで
`--post-min-event-year 2022` に直す: 主の表（全ての AOI で抜いたセル）を作った後で、事象が 2022 年より前の AOI の行を
除き、場所グループと fold はその後で作り直す（2022 年以降の事象は 210 AOI、581,512 セル、69 activation）。他は
追記 1 のとおり

## 追記 3（B1〜B5・B7・B13 の実装の記録。2026-10-29T12:47Z、commit 03c120e、B の量を計算する前）

1. **B1(c)**: 二つ以上の国にまたがる activation は、セルが最も多い国に入れる（入れ子の再抽出のため）
2. **B1 の区間**: 片側の限界は再標本の L± の 5%・95% 点。Fieller の区間は再標本の ΔBS と Δ_base の分散・共分散
   からの正規近似（90%）。比の百分位の区間も並べる。(b)(c) の点は重みの無い当て直し（= 直した実行）
3. **B1(d)**: L± の分散 = (c) の再標本の Var(ΔBS)・Var(Δ_base)・Cov から作る分散 + 13 の割り当ての L± の
   Nadeau–Bengio の分散（n_test/n_train = 1/9）。Fieller の区間も同じく足した分散で
4. **B2**: 上限は R の尺度で、直した実行の物差しを掛けて Brier の尺度にし、無効果の基準 ΔBS_null を引いた値に平均の
   利得が届く OR（v59 の bound の読み）。形 A と B。B1 の上限は (a) 直した実行、(b)、(c)、(d) の Fieller の上限。
   推定量の尺度は v68 の B3 の規則（B-10 の B11 の曲線）
5. **B3**: 符号化者 2 の coded の行で、水に関わる = 水系側の五つの区分のどれか（v57b と同じ）。Wilson の 95% 区間。
   確率的な偏りの分析は v57 の pairs の再標本（保存したもの）に Beta(k+1, n−k+1) の率を両方の符号化者について
   （20,000 回）
6. **B5**: 基準（M1c）と名のモデル（M3n）の両方から文字モデルの得点の AOI 平均を除く（二段目は B-8 の C2 の設計で
   当て直す）。比べる相手として同じ当て直しの標準の R も
7. **B7**: Y3 = 定数の予測の Brier − M1 の Brier。定数は、主分割では fold ごとの訓練側の浸水の割合、一つの分割の
   llmcut では評価側の割合（B-12 の B17）。Brier skill の増分 = ΔBS ÷ その定数の Brier
8. **B13**: 浸水と乾きの両方がある AOI ごとに、M3n と M1c の AUC・PR-AUC の差。Brier skill の増分は AOI 自身の浸水の
   割合 p で p(1 − p) を気候値として AOI ごとに計算し、AOI に同じ重みで平均と中央値

## 追記 4（B6・B8〜B12・B14・B15 の実装の記録。2026-10-29T12:57Z、commit c394b67、これらの量を計算する前）

スクリプトは `backend/verification/r_okauchi/v70b_b13_part2.py`（出力 `docs/verification/v70_b13_local2.json`）。

1. **B6**: 影響は場所グループの jackknife（モデル固定）の ψ = R − R(−g)。実効のクラスタ数 = (Σψ²)² ÷ Σψ⁴。確証の
   実行と直した実行のそれぞれで影響の大きい 10 グループを一つずつ除いた R。B1(b) の側は、直した実行の影響の順の
   上位 10 を、二段目を当て直して（B-8 の C2 の設計）除く
2. **B8**: 各 fold の評価側のセル（全てのセル）で、[基礎の logit・地域構成・名の特徴] を標準化した説明変数のロジスティック
   回帰（切片一つ、v39.fit）。場所グループでクラスタ頑健。値が一定の列は外す。名の係数が全て 0 の同時 Wald 検定
3. **B9**: 上限なしの EMS の推定（B-11 の B2 と同じ表・直した読み・pairs の再抽出 2,000 回）と計画した GFM の推定
   （v41 の表、クラスタ頑健）。(a) 最寄りの集落の GeoNames の座標に v37 の丸めの規則。(b) 500 m 以内に水系側の
   セルがある activation を一つずつ除いた点（度数の重みの IRLS）と、EMSR052・EMSR130 を除いた pairs の区間。
   (c) 記録の schema_generation が 0・1 の AOI を除く。(d) 各 AOI の中で、（AOI, 最寄りの集落）の組の水系側の印を
   組の間で入れ替え（200 回）、上限なしの EMS の点（IRLS）と計画した GFM の点。p は |log OR| の両側で
   (1 + 数) ÷ (1 + 200)。(e) 計画した GFM の表のセルで、距離帯と 10 m の市街地の割合（WorldCover の lc_built）の
   三分位。(f) 上限の表（v41、EMS のラベル）で、センサーの分かる AOI（Sentinel-1 と他の SAR = レーダー、光学）、
   水系側 × 500 m 以内 × 光学の項。(g) 最寄りの集落の地物の符号が PPLX・PPLH・PPLQ・PPLW のセルを除く（両方の推定）
4. **B10**: 解析する集落 = 上限なしの表でいずれかのセルの最寄りの集落。凍結した辞書を直した読み込みで。位置は照合した
   項目ごとに数える。Unicode の名は ä・ö・ü・å・ø を一文字に、æ → ae、ß → ss、結合の印を外して読む
5. **B11**: B10 と同じ集落の N-read の区分の割合を、国に当てたブロック（+ でつなぐ）ごとと国ごとに
6. **B12**: 両方の符号化者の符号を GeoNames の ID・名・国・層と並べた公開のファイル `docs/verification/v70_two_coders.csv`。
   どれかの照合の標本（v51・v68・v69）にある名は除く（照合の答えを漏らさないため）。一致は v57b の値
7. **B14**: 計画した GFM の表の AOI ごとの EMS と GFM の浸水の割合。GFM のラベルで学習し採点した実行（world_v2_gfm）の
   AOI ごとの地形の利得（M0 の Brier − M1 の Brier）
8. **B15**: 追記 1 の 8 のとおり。集落の値は上限の表（v37 の表、v41 と同じセル）で、最寄りの集落がその集落のセルから
   作る

## 追記 5（C の集計と、直した辞書での集落の近くのオッズ比の実装の記録。2026-10-29T13:14Z、commit 44fbbfc、これらの量を計算する前）

スクリプトは `backend/verification/r_okauchi/v71_b13_cluster.py`（出力 `docs/verification/v71_b13_cluster.json`）。

1. **C1〜C3 の集計**: 各実行の held-out の予測から、確証の解析と同じ場所グループの bootstrap（4,000 回、v56 の
   `_verdicts`）で R・区間・判定表の判定。C1 は名の主比較（M1c → M3n、placebo M3np）と各対照（セルに直接と集落を
   通じて）。C3 は加えて、fold ごとの訓練側の浸水の割合を予測とする定数の Brier と、選ばれた反復の回数の分布
2. **C2 の集落の近くのオッズ比**: 最寄りの集落の名を直した辞書で読み直した水系側の印で、上限ありの EMS（v37 の表、
   第 5 節の推定量、クラスタ頑健）、計画した GFM の推定（v41 の表）、上限なしの EMS（B-11 の B2 と同じ表と pairs の
   再抽出 2,000 回）
3. **C4（B-12 の C1）**: 形 BK の利得の曲線から、余白に届く注入の OR（v59 の reach）と、固定と当て直しの区間の上限が
   届く OR（無効果の 4 つの基準ごと）。v64 の C3 と同じ読み

## 追記 6（B1(c) の実装の直しと、C の集計の記録。2026-10-29T14:20Z、commit 07c9788、B1(b)(c)(d)・B2 を計算する前）

1. **B1(c) の国の引き方**: 国の符号の表（v68 の B7 の符号）には 70 の国があるが、activation をセルの最も多い国に入れる
   （追記 3 の 1）と、activation を持つ国は 55 で、残る 15 は他の国に入れた activation の AOI にセルがあるだけの国。
   最初の実装は 70 の国から引き、activation の無い国を引いた回で 0 で割って止まった（約 95 分の計算の後、結果は保存
   されていない）。追記 1 の 7 の「国を等確率で国の数だけ引き」は、activation を持つ 55 の国を等確率で 55 回引くことと
   して直す（`topsem.world.b13.nested_weights`）。引いた国の中の activation の引き直しは追記 3 のまま
2. **途中の保存**: (b) と (c) の再標本を方式ごとに保存し（`data/verification/v70_b1_refit_reps_{groups,nested}.npy`）、
   保存があれば読む。回数・種・当て直しの規則は変えない（(b) の 500 回の種は SEED + b、(c) は SEED + 100,000 + b）
3. **C1〜C4 の集計の記録**: v71 の C3 の集計は予測のファイルに fold の番号が無く止まったので、load_split と同じ順
   （ファイル名の順）にファイルの長さから fold の番号を作った（定数の予測の規則は追記 5 のまま）。C1(b) の予測は NAS
   からの引き上げで fold8 が途中で切れたので、その二つのファイルを引き上げ直し、全ての fold が読めることを確かめてから
   集計した

## 追記 7（B9 の pairs の再抽出の一部をクラスタで。2026-10-29T16:40Z、commit f25fffe、B9(b)(c)(g) を計算する前）

1. **場所**: B9 の (b) EMSR052・EMSR130 を除いた表、(c) 世代 0・1 を除いた表、(g) PPLX・PPLH・PPLQ・PPLW を除いた表の
   pairs の再抽出（v68.pairs_or、2,000 回、種は v68 の SEED）を、手元では 2 並列で 1 回あたり 2 時間を超えるので、
   クラスタで計算する。表と部分集合の作り方・種・回数は v70b と同じ（`backend/verification/r_okauchi/v72_b9_pairs_cluster.py`
   が v70b の関数で同じ部分集合を作り、v68.pairs_or を同じ鍵で呼ぶ）。上限なしの表は NAS の
   `data/verification/v68_uncapped.parquet`（手元と同じファイル。sha256 を記録する）
2. **集計**: 各回の結果（`data/verification/v68_pairs_v70_b9{b_no052_130,c_gen2plus,g_fcodes}.npy`）を NAS から引き上げ、
   手元の v70b の B9 がそれを読む（同じ鍵の保存があれば再計算しない作り）。(a) と (b) の活性化の一つ抜き、(d) は手元で計算する
3. **計算機の違い**: クラスタは amd64、手元は arm64。IRLS の浮動小数点の結果は機械でわずかに違いうるので、手元で計算した
   (a) との違いは区間の端の数字に出ない程度と見込むが、その確認はしない（主比較の検証では両者が R の小数 5 桁まで一致した）
4. **実装の記録**: 最初の投入（イメージ 1793292010673138000）は、v36 がリポジトリの根を CLAUDE.md で探すのにイメージに
   CLAUDE.md が無く、読み込みで止まった（何も計算していない）。Dockerfile で CLAUDE.md を写し、同じ引数で投入し直した
