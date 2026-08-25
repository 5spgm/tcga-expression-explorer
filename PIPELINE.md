# TCGA Expression Explorer

TCGAのRNA-seqデータを、旧・中期・新の3世代パイプラインで並べてboxplot比較するための
静的Webツール一式。GitHub Pagesでの公開を想定している。

## 構成

```
tcga_tool/
  preprocess_tcga.py   前処理スクリプト(手元/サーバーで実行し、data/を生成する)
  site/                GitHub Pagesにそのまま置く静的サイト
    index.html
    style.css
    app.js
    data/              前処理済みJSON(今はCOADのテストデータのみ同梱)
```

## 1. 前処理スクリプトを実行する

各がん種ごとに、旧・中期・新のファイルパスを指定して実行する。

```bash
pip install pandas openpyxl xlrd mygene

python3 scripts/preprocess_tcga.py \
    --cancer-type COAD \
    --new-tpm      /path/to/20221216_COAD_TPM.xls \
    --mid-fpkmuq   /path/to/COAD_FPKM_UQ.xls \
    --old-normcount /path/to/COAD-HiseqV2-20150129.txt \
    --out-dir ./site/data
```

がん種の数だけこれを繰り返す(`--cancer-type` と各ファイルパスを差し替えるだけ)。
`--entrez-cache` は省略時 `./entrez_cache.json` に作られ、全がん種で使い回せる
(同じEntrez IDは何度も問い合わせない)。

> **注意: スクリプトが完了するまで `site/data` へコピー(または閲覧)しないこと。**
> 遺伝子ごとのJSONは処理しながら1つずつ書き出されるが、`_index.json` と
> `manifest.json` は全遺伝子を処理し終えた最後にまとめて書き込まれる。
> 実行中にコピーすると、一部の遺伝子だけ書き出し済みだったり、特定の世代
> (新/中期/旧)のデータが途中までしか揃っていなかったりする不完全な状態
> になり、「特定の世代だけ表示されない」といった紛らわしい症状の原因になる。
> ターミナルに `完了: N遺伝子分のJSONを書き出しました` と出るまで待つこと。

> **既知の問題と対策(2026年8月時点で解消済み)**: 実データ規模(数万遺伝子×
> 数百検体)で「約2.9秒/遺伝子」という極端な遅さと、その後の修正でも
> segmentation faultが発生する問題があった。原因は、遺伝子ごとの値を
> wide(遺伝子×検体)からlong(1行1組)形式に変換して処理していたため、
> 数千万行規模に膨れ上がりメモリを圧迫していたこと。現在のバージョンは
> wide行列のままnumpyの2次元配列として保持し、遺伝子ごとに行を1回
> スライスするだけの実装に変更済み(6万遺伝子×501検体の合成データで
> 約60秒、クラッシュなしで完走することを確認)。

> **もう一つの既知の問題(2026年8月時点で解消済み)**: 中期データ(FPKM-UQ)
> で、期待される約6万遺伝子のうち2万遺伝子ほどしか出力されず、ACTBのような
> ハウスキーピング遺伝子すら見つからないことがあった。原因は、このファイルが
> 拡張子こそ`.xls`だが実体はR の `write.table()` のデフォルト書式(1列目=
> row.namesにヘッダー名を付けない)のTSVであり、ヘッダー行の項目数がデータ行
> より1つ少なかったこと。pandasがこれを「先頭列はインデックスだろう」と
> 自動推測してしまい、遺伝子ID列がDataFrameの列ではなくインデックスに
> 吸い込まれ、代わりに最初のサンプル列が誤って遺伝子ID列として扱われて
> いた(結果、大半の遺伝子が消失・破損していた)。現在のバージョンは
> `read_matrix_file()` 内でこのパターンを検出し、自動的に列を復元する。

## 1.5. Tumorを分子サブタイプで分割表示する(任意・複数分類を切り替え可能)

Tumor検体を、選んだ分類方式(iCluster、MSI、CIMPなど)ごとに複数の箱に分けて
表示できる。画面上部の「Tumor分類」プルダウンで方式を切り替えられる(「なし」
を選べば従来通りTumor/Normalの2群)。

### 対応している分類

| 出典 | 分類 | 内容 |
|---|---|---|
| Hoadley et al. 2018 (Cell 173(2):291-304.e6) | iCluster | 全33がん種横断の統合クラスタ(28クラスタ) |
| Liu et al. 2018 (Cell 173(4):963-985.e16) | Molecular Subtype | 消化管腺癌のメイン分類(CIN/MSI/GS/HM-SNV/EBV) |
| 同上 | MSI Status | マイクロサテライト不安定性(MSS/MSI-H/MSI-L) |
| 同上 | CIMP | Hypermethylation category。CIMP-H/CRC・GEA CIMP-L/Non-CIMP/CIMP EBV。この表現型はIssa先生らが提唱した概念に基づく分類 |
| 同上 | Colorectal CMS | Guinney et al. 2015のCMS1-4(大腸がんのみ) |
| BCR clinical patient ファイル | TNBC | 免疫組織化学によるER/PR/HER2から Triple-negative / Non-triple-negative(乳がんのみ) |
| Thorsson et al. 2018 (Immunity 48(4):812-830.e14) | Immune Subtype | 全33がん種横断の免疫サブタイプ C1-C6 |
| 同上 (TCGA Subtype列) | 各がん種の既発表サブタイプ | 乳がんは**PAM50**(LumA/LumB/Basal/Her2/Normal)、消化管はCIN/MSI/GS/EBV など |

### 手順

1. 論文のSupplementary Tableを入手する
   - Hoadley 2018: Table S6 (`mmc6.xlsx`)
   - Liu 2018: Table S1 (`mmc2.xlsx`, "Master Patient Table"シート)
2. 汎用フォーマット(`patient_id,subtype`の2列CSV)に変換する:
   ```bash
   python3 scripts/make_icluster_subtype_table.py --input mmc6.xlsx --output icluster_subtypes.csv
   python3 scripts/make_liu_subtype_tables.py --input mmc2.xlsx --out-dir ./subtypes
   ```
   後者は `subtypes/Molecular_Subtype.csv`, `MSI_Status.csv`, `CIMP.csv`,
   `Colorectal_CMS.csv` の4ファイルを一度に生成する。

   乳がんのTNBC分類は、論文ではなくBCRのclinical patientファイルから作る:
   ```bash
   python3 scripts/make_tnbc_subtype_table.py \
       --input nationwidechildrens_org_clinical_patient_brca.txt \
       --output subtypes/BRCA_TNBC.csv \
       --detail-csv subtypes/BRCA_TNBC_detail.csv
   ```
   判定は ER(`er_status_by_ihc`)、PR(`pr_status_by_ihc`)、HER2
   (`her2_status_by_ihc`。Equivocal/Indeterminate/未評価のときだけ
   `her2_fish_status` で補う)による。3つすべてNegativeならTriple-negative、
   1つでもPositiveならNon-triple-negative(残りが未確定でも確定できる)。
   陽性ゼロだが未確定項目がある患者は判定不能としてCSVに出さない
   (`--emit-unknown` で 'Unknown' として含めることもできる)。
   `--detail-csv` に患者ごとのER/PR/HER2の内訳が出るので、判定の妥当性は
   そちらで確認できる。

   **この分類は2値なので `--min-subtype-n` は `:1` にすること**
   (既定の5でも実害はないが、明示しておくと安全)。

   Thorsson 2018 の Supplementary Table (`mmc2.xlsx`, シート `PanImmune_MS`)
   からは、1ファイルで**2種類**の分類が取り出せる:
   ```bash
   python3 scripts/make_immune_subtype_tables.py --input mmc2.xlsx --out-dir ./subtypes
   ```
   - `subtypes/Immune_Subtype.csv` — 免疫サブタイプ C1-C6。全33がん種横断
     なので、iClusterと同じくどのがん種でも使い回せる(9,126患者)。
     ラベルには論文の呼称を併記する(例: `C2 IFN-gamma Dominant`)。
     素の `C1`〜`C6` にしたい場合は `--raw-immune-labels`。
   - `subtypes/tcga_subtype/<接頭辞>.csv` — がん種ごとの既発表サブタイプ。
     **乳がん(`BRCA.csv`)はPAM50**。接頭辞はTCGAのstudy略号と一致しない
     ことがある(`GI.csv` は COAD/READ/STAD/ESCA をまとめたもの、
     `GBM_LGG.csv` は GBM/LGG、`OVCA.csv` は OV)。どのがん種にどのファイルを
     渡すかは、生成時に表示される対応表で確認すること。

   > PAM50の `Normal`(Normal-like)は、正常組織の混入によるアーティファクト
   > の可能性が指摘されている群。解釈の際は留意すること。
3. `scripts/preprocess_tcga.py` に `--tumor-subtype-table "表示名=path.csv"` を
   **必要な数だけ繰り返し指定**する(表示名がプルダウンの選択肢になる)。
   末尾に `:数字` を付けると、そのスキームだけ専用の `--min-subtype-n` を
   指定できる(省略時は共通の `--min-subtype-n` に従う):
   ```bash
   python3 scripts/preprocess_tcga.py \
       --cancer-type COAD \
       --new-tpm ... --mid-fpkmuq ... --old-normcount ... \
       --tumor-subtype-table "iCluster (Hoadley 2018)=icluster_subtypes.csv" \
       --tumor-subtype-table "Molecular Subtype (Liu 2018)=subtypes/Molecular_Subtype.csv:1" \
       --tumor-subtype-table "MSI Status (Liu 2018)=subtypes/MSI_Status.csv:1" \
       --tumor-subtype-table "CIMP (Liu 2018)=subtypes/CIMP.csv:1" \
       --tumor-subtype-table "Colorectal CMS (Liu 2018)=subtypes/Colorectal_CMS.csv:2" \
       --min-subtype-n 5 \
       --out-dir ./site/data
   ```
   `--min-subtype-n`(共通既定値、既定5)は「このがん種内での出現数が
   これ未満のsubtypeを自動的に'Other'にまとめる」閾値。**カテゴリー数が
   多い分類(iClusterの28クラスタなど)には有効だが、もともと3〜6種類しか
   ないMSI/CIMP/Molecular Subtype/CMSのような分類では、既定の5だと少数派の
   群(CIMP-HやMSI-Hなど、むしろ見たい群であることが多い)まで'Other'に
   吸収されてしまう**。その場合はスキームごとに`:1`や`:2`のように緩めるか、
   実質OFFにするとよい。
4. `icluster_subtypes.csv` と `subtypes/*.csv` はどちらも全がん種共通なので、
   他のがん種を処理するときも同じファイルをそのまま使い回せる(該当する
   患者だけが自動的に拾われる)。COAD/READ以外のがん種ではColorectal CMSは
   該当患者がいないため、その方式はプルダウンに出てこない。

**注意点**
- `mygene.info` への問い合わせにはこのスクリプトを実行するマシンからの
  インターネット接続が必要。社内ネットワークで塞がれている場合は、
  別マシンで動かすか、VPN/プロキシ設定を確認すること。
- 実行後に画面に出る `symbol直接一致: n / N 遺伝子` と
  `未マッピングの遺伝子 n件` は必ず確認する。想定より低い場合は
  遺伝子ID形式が想定と違う可能性がある。
- `sample type code 内訳` も毎回確認する。コードごとに Tumor / TumorExtra /
  Normal のどれに振り分けられたかが表示される(下の「1.6.」参照)。想定外の
  コードがあれば警告が出る。
- `--mid-fpkm` や `--new-fpkm` / `--new-fpkmuq` を追加すれば、サイト側に
  自動で選択ボタンが増える(コード変更不要)。
- **中期データは、がん種によって手元にある値の種類が違う**(COADは
  FPKM-UQ、BRCAはFPKM)。`--mid-fpkmuq` と `--mid-fpkm` は別の引数なので、
  ファイルに合った方を指定すること。指定した値の種類(`FPKM` / `FPKM_UQ`)は
  パネル上部に常に表示され、y軸ラベルにも入るので、混同する心配はない。
  両方ある場合はボタンで切り替えられる。
- **中期データがそもそも無いがん種もある**(2022年のDR32更新以前に
  取得していなかった場合)。`--mid-*` を省略すればそのまま実行でき、
  サイト側では中期パネルに「このがん種では未取得のデータです」と表示される
  (遺伝子IDが対応付かなかっただけの場合は別の文言になるので区別できる)。

## 1.6. 転移巣・再発巣の扱い(既定は除外、チェックボックスで合流)

TCGAバーコードの sample type code を3群に振り分ける。

| 群 | コード | 内容 | 既定の表示 |
|---|---|---|---|
| `Tumor` | 01, 03, 05, 09 | 原発巣(Primary) | Tumor群として表示 |
| `TumorExtra` | 02, 04, 06, 07, 08, 40 | 再発巣・転移巣など | **表示しない** |
| `Normal` | 10〜14 | 非腫瘍部 | Normal群として表示 |

転移巣・再発巣は遺伝子JSONの `tumor_extra`(および分類方式別の
`tumor_extra_by_scheme`)に**別枠で**書き出される。原発巣の `tumor` とは
混ざらないので、既定では原発巣だけの分布が描かれる。

サイト側の**「転移巣・再発巣を含む」チェックボックス**をONにすると、
`tumor_extra` が `tumor` に合流して1つの箱として描画される(分類方式を
選んでいる場合は、subtypeごとの箱それぞれに合流する)。脚注には
`[転移巣・再発巣 7件を含む]` / `[... を除外]` と現在の状態が表示される。

チェックボックスは `manifest.json` の `extra_tumor.n_samples` を見て
自動で有効・無効が切り替わる。該当検体が無いがん種(例: COAD)では
グレーアウトし、ラベル横に件数と内訳がツールチップで出る。

想定外のコードが出た場合は警告を出したうえで、腫瘍系(コード10未満)なら
安全側に倒して `TumorExtra`(既定で除外)に入れる。振り分けを変えたい場合は
スクリプト冒頭の `PRIMARY_TUMOR_CODES` / `EXTRA_TUMOR_CODES` /
`NORMAL_CODES` を編集する。

> **注意**: この3分割は既存の出力JSONとは形式が違う。`tumor_extra` を
> 持たない古いJSON(この機能より前に生成したもの)を読んでもエラーには
> ならず、単にチェックボックスが無効化されるだけだが、その場合の `tumor`
> には転移巣・再発巣が混ざったままになっている。**該当コードを含むがん種は
> 前処理をやり直すこと。**

## 1.7. 同一患者の重複vial(01A / 01B)とFFPE検体

TCGAバーコードの4番目のフィールド末尾の文字(`TCGA-A7-A13D-01**A**`)は
**vial文字**で、同じ検体から採られた何本目のバイアルかを示す。同じ患者・
同じsample typeで `01A` と `01B` の両方がデータに含まれていると、その患者が
boxplotに2回カウントされてしまう。

### `--dedup-vials`(既定 `first`)

同一患者×同一sample typeに複数vialがある場合、1つだけを採用する。
既定はvial文字が若い方(A)。`last` で遅い方、`none` で除去しない。
**ファイルごとに独立して適用される**(どのvialが入っているかは世代によって違う)。

実行時に、どの患者のどのvialを落としたかが一覧表示されるので必ず確認すること。

### `--exclude-vials`(既定: 空)

指定したvial文字を**重複の有無に関わらず**落とす。例: `--exclude-vials B`。

> **強く注意**: 「Bはとりあえず全部消す」は多くの場合まずい。BRCAの実データで
> 確認したところ、`B` vialのうち大半は**その患者にとって唯一の検体**で
> (新世代TPMで28件、旧世代で27件)、重複しているのは6件だけだった。
> 一律除去するとその27〜28人がまるごと解析から消える。まず既定の
> `--dedup-vials first` だけで実行し、表示される重複組の一覧を見てから
> 判断すること。

### FFPE検体を落としたい場合

**vial文字はFFPE由来かどうかを表さない。** `01B` が必ずFFPEということはなく、
逆にFFPE検体が `01A` のこともある。FFPEかどうかはGDCのメタデータの
`is_ffpe` フラグが唯一の正確な情報源で、バーコードからは判別できない。

GDCのsample sheet(またはGDC APIの `samples.is_ffpe`)からFFPE検体の
バーコード一覧を作り、`--exclude-samples` に渡す。

> **重要**: FFPE情報は**検体レベル**の属性なので、`nationwidechildrens_org_
> clinical_patient_*.txt`(患者レベル)には**入っていない**。同じBCRアーカイブの
> **biospecimen sample** ファイル
> (`nationwidechildrens_org_biospecimen_sample_brca.txt`、`is_ffpe` 列)
> か、GDCのsample sheetを使うこと。

同梱の `scripts/make_ffpe_exclude_list.py` が、BCR biospecimen形式・GDC sample sheet
形式のどちらでも列名を自動検出して除外リストを作る(3行ヘッダーのBCR形式にも
対応済み):

```bash
python3 scripts/make_ffpe_exclude_list.py \
    --input nationwidechildrens_org_biospecimen_sample_brca.txt \
    --output ffpe_exclude_brca.txt
# 列名が想定と違う場合は --barcode-column / --ffpe-column で明示指定できる

python3 scripts/preprocess_tcga.py \
    --cancer-type BRCA \
    --new-tpm ... --mid-fpkm ... --old-normcount ... \
    --exclude-samples ffpe_exclude_brca.txt \
    --out-dir ./site/data
```

`--exclude-samples` は1行1件のテキストで、フルバーコード
(`TCGA-XX-XXXX-01A`)でも患者ID(`TCGA-XX-XXXX`)でも受け付ける。
先頭行が `barcode` / `sample_id` / `patient_id` ならヘッダーとして読み飛ばす。

適用した設定は `manifest.json` の `filters` に記録されるので、後から
どの条件で作ったデータか追跡できる。

## 1.8. 新世代データの検体リスト作成を検証する(check_rnaseq_list.py)

新世代データをGDCから落として1つの行列にまとめる際、GDC sample sheetから
「どのファイルを読むか」を選ぶ工程がある。ここでの選び方を誤ると、特定の
検体や患者が**警告なしに**消えるが、出来上がった行列を見ただけでは気づけない。

`scripts/check_rnaseq_list.py` は、GDC sample sheet を読んで
`20250913_RNAseq_list.R` の選別ロジックを再現し、どの検体・どの患者が
落ちるかを一覧表示したうえで、修正版の検体リストを出力する。

```bash
python3 scripts/check_rnaseq_list.py \
    --sample-sheet gdc_sample_sheet.2025-09-13.tsv \
    --output corrected_file_list.tsv
```

検出する3つの落とし穴:

1. **`sheet[,8]=="Tumor"` の完全一致**
   GDC標準のSample Type値は `Primary Tumor` / `Solid Tissue Normal` /
   `Metastatic` などなので、そのままでは一致しない。特に `Metastatic` は
   Tumorにも Normalにも一致せず、必ず落ちる。
2. **`grep("-01A", ...)`**
   コメントは「01Aと01Bが両方ある場合、01Aを落とす」だが、コードは逆に
   01Aだけを残している。さらに、複数ファイルを持つ患者で 01A が無い場合、
   その患者はまるごと消える。ファイルが1つだけの患者はこのgrepを通らない
   ため01Bでも残る、という非対称性があり発見しにくい。
3. **`!duplicated(sheet3[,6])`**
   6列目は Case ID(患者)なので患者単位の重複除去になっている。
   意図はおそらく Sample ID(7列目)単位。

出力された `corrected_file_list.tsv` を R側の `sheet` の代わりに読ませれば、
取りこぼしなく作り直せる。

## 2. サイトをローカルで確認する

ブラウザの `fetch()` はfileプロトコルだとCORSでブロックされるため、
簡易サーバー経由で開く。

```bash
cd site
python3 -m http.server 8000
# ブラウザで http://localhost:8000 を開く
```

同梱の `data/COAD/` はアップロードいただいたサンプル(先頭100行程度)から
生成したテストデータなので、遺伝子数は99件のみ。実データで再生成すれば
全遺伝子が使えるようになる。

## 3. データを自前サーバーに置く場合

`site/app.js` の先頭にある `DATA_BASE_URL` を書き換える。末尾にスラッシュは付けない。

```js
const DATA_BASE_URL = "https://tcga-data.example.jp/data";
```

自前サーバー側で対応が必要な点は2つ。

**(a) CORS**: レスポンスに以下のヘッダが必要(GitHub Pagesのドメインから読み込むため):
```
Access-Control-Allow-Origin: *
```
(あるいは `https://<username>.github.io` を明示的に許可)

**(b) HTTPS化**: GitHub Pagesは常にHTTPSで配信されるため、HTTPの外部リソースは
mixed contentとしてブラウザにブロックされる。自前サーバー側もHTTPS化が必須。

### HTTPS化の方法

- **固定グローバルIP + ポート開放ができる場合**: nginx + certbot (Let's Encrypt)
  で通常通りHTTPS化する。
- **大学・組織内ネットワークなどでポート開放が難しい/許可が下りない場合**:
  [Cloudflare Tunnel](https://developers.cloudflare.com/cloudflare-one/connections/connect-networks/)
  (`cloudflared`)を使うと、ポート開放なし(アウトバウンド接続のみ)でHTTPSの
  公開URLが得られる。まず `cloudflared tunnel --url http://localhost:8080`
  でクイックトンネル(一時的な `*.trycloudflare.com` URL)を試し、恒久運用
  するなら独自ドメインで named tunnel を作るとよい。

nginxの設定例(ローカルで `site/data/` を配信し、CORSヘッダを付与):
```nginx
server {
    listen 8080;  # ポート開放する場合は 80/443 に、certbotと組み合わせる
    server_name _;
    root /path/to/tcga_expression_explorer/site;

    location /data/ {
        add_header Access-Control-Allow-Origin "*" always;
        add_header Access-Control-Allow-Methods "GET, OPTIONS" always;
        autoindex off;
    }
}
```

CORS/HTTPS対応が難しい場合は、`site/data/` にそのままJSON一式をコピーして
リポジトリに含めてしまう方法もある(容量が許せば一番シンプル。ただし
`site/.gitignore` で `data/` を除外しているので、その場合は`.gitignore`から
該当行を削除すること)。

## 4. GitHub Pagesへの公開

`site/` の中身をリポジトリのルート(またはdocs/)に置き、
Settings → Pages でブランチを指定するだけ。`data/` を含めない場合は
上記3.の設定が必要。

## 未対応・今後の課題

- 旧データの `symbol` 不明("?")遺伝子のEntrez解決は、mygene.infoの
  ネット接続が使える環境で再実行して確認すること。
- サンプルタイプコード `02`(Recurrent)など少数派コードの扱い方針は
  データを見ながら要調整。
- 30種類以上のがん種すべてを処理すると `data/` 全体のファイル数は
  数十万に達する見込み(がん種数 × 遺伝子数)。GitHubの1リポジトリ
  あたりのファイル数上限やcloneの遅さが気になる場合は、
  自前サーバー配信(3.)を優先すること。
