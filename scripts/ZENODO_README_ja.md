# TCGA Expression Explorer — 解析に用いたデータ一式

TCGAのRNA-seqデータを、3世代の定量パイプラインで比較するために整えたものです。
可視化ツール **TCGA Expression Explorer** の入力そのものであり、
公開されている図はすべてこのデータから再現できます。

- ツール: https://5spgm.github.io/tcga-expression-explorer/
- コード: https://github.com/5spgm/tcga-expression-explorer

---

## なぜこのデータを寄託するのか

TCGAのRNA-seq定量パイプラインは3回入れ替わっています。**中期(HTSeq)の値は
2022年のData Release 32でSTAR-Countsに置き換えられており、現行のGDCから
同じものを取得することはできません。**

| 世代 | パイプライン | 値 | 遺伝子行数 | 取得元 |
|---|---|---|---|---|
| 旧世代 | UNC IlluminaHiSeq RNASeqV2 | normalized count | 20,531 | 旧TCGA Data Coordinating Center (http://cancergenome.nih.gov/、閉鎖) |
| 中期 | GDC HTSeq (DR15–31) | FPKM / FPKM-UQ | 60,483 | 当時のGDC。**現行のGDCには無い** |
| 新世代 | GDC STAR-Counts (DR32+, 2022–) | TPM | 60,660 | GDC(現行) |

GDCはData Releaseのたびに過去世代の定量値を置き換えるため、DR15–31の期間に
取得した人の手元にしか中期のデータは残っていません。3世代の変遷を検証する
手段が失われつつある、というのが寄託の主な動機です。

中期について例外が1つあります。本パッケージの乳がん中期は、UCSC Xena の
GDC hub から **2024-04-13** に `TCGA-BRCA.htseq_fpkm.tsv.gz` として取得した
ものです(取得日は手元の .gz のタイムスタンプ。gzipヘッダには配布側の
ファイル時刻 2019-07-18 22:22:54 UTC が残っています。展開した tsv の SHA-256 はこの .gz を
展開したストリームと一致します)。
配布されているのは log2(FPKM+1) に変換された値で、生のFPKMではありません
(下記「収録内容」を参照)。それ以外のがん種の中期は、当時ダウンロードした
ものが唯一の入手経路です。

その後、この hub のカタログはSTAR期(gencode v36)に入れ替わり、
**コホート単位の HTSeq FPKM (`TCGA-BRCA.htseq_fpkm`) は削除されました。**
2026-09-08時点で hub に残っている HTSeq FPKM は pan-cancer の
`GDC-PANCAN.htseq_fpkm` (gencode v22) だけです。本パッケージの乳がん中期
1,217検体はすべてこれに含まれ、値も hub のAPIが返す桁の範囲で一致することを
確認しました。したがって同じ値は現在も次の経路で取得できます。

    hub      https://gdc.xenahubs.net
    dataset  GDC-PANCAN.htseq_fpkm.tsv   (gencode v22, log2(FPKM+1))
    page     https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.htseq_fpkm.tsv
             &host=https%3A%2F%2Fgdc.xenahubs.net
    手順     pan-cancer行列から、本パッケージの乳がん検体の列を抜き出す
             (検体一覧は sample_annotation/BRCA_samples.tsv)

旧世代は、**すべて旧TCGAポータル (http://cancergenome.nih.gov/) から
取得しました。** このポータルは閉鎖され、TCGAデータはNIH GDCに移管されています。
同等の値はUCSC Xenaなどが再配布していますが、遺伝子IDが記号のみになっており、
本パッケージが保持している "SYMBOL|EntrezID" とは異なります。
新世代は現行のGDCから取得できます。旧世代と新世代をここに含めているのは、
**公開図を作った正確な入力**を残すためです。

---

## 収録内容

```
matrix/<がん種>/<がん種>_<世代>_<値の種類>.tsv.gz
    発現行列。取得したままの形で置いています。
    世代ごとに遺伝子ID体系も検体集合も違うため、1つの表には統合していません
    (統合するとNAが大半を占め、ID変換という解釈が元データに混ざります)。

      旧世代  1列目が "SYMBOL|EntrezID"、20,531行 (hg19、GENCODE以前)
      中期    1列目がEnsembl ID、60,483行 (GENCODE v22)
      新世代  gene_id / gene_name / gene_type の3列 + 検体列、60,660行 (GENCODE v36)

    遺伝子行数は5がん種すべてで一致するので、ファイル単体でも世代を
    識別できます。

    乳がんの中期だけは値が log2(FPKM+1) です。UCSC Xena GDC hub が
    log2変換した形で配布しているものを、変換せずそのまま置いています。
    ファイル名も BRCA_mid_FPKM_log2p1.tsv.gz として中身に合わせています。
    生のFPKMに戻すには 2^x - 1 を適用してください。他のがん種の中期は
    変換されていないFPKM-UQです。

    この1ファイルだけは、UCSC Xenaから取得した .gz を再圧縮も書き換えも
    せずbyte-for-byteで置いています(SHA-256 b5856c84…)。そのため
      - 1列目の名前が他の14ファイルの gene_id ではなく Ensembl_ID
      - gzipヘッダに配布側のファイル名 TCGA-BRCA.htseq_fpkm.tsv と
        時刻 2019-07-18 22:22:54 UTC が残っている (gzip -l で確認可。
        gzip はヘッダにUTCの秒数を格納し、gzip -l は実行環境のローカル時刻で
        表示するため、JSTでは 2019-07-19 07:22:54 と出る)
    という点が他と異なります。公開サイトの図は、この .gz を展開して
    xlsxに変換したものを入力に作りました。両者の値の差は先頭100遺伝子
    ×全1,217検体で最大 3.2e-14 (倍精度の丸め)で、表示桁には現れません。

sample_annotation/<がん種>_samples.tsv
    検体ごとに1行。どの世代に存在するか、どの群に属するか、
    どの分類が割り当てられたかを1枚で追えます。
    世代によってバーコードの粒度が違う(旧世代はvial文字なし)ため、
    「患者ID + sample type code」で対応付けています。

exclude/
    解析から外した検体の一覧と、その根拠。

provenance/<がん種>/<世代>/
    ダウンロードした時点のGDC sample sheet と manifest。旧世代については
    旧TCGA DCCの一式 (file_manifest.txt / FILE_SAMPLE_MAP.txt / README_DCC.txt、
    胃がんと膵臓がんではさらに file_annotations.txt)。
    file_annotations.txt はDCCが付けた注記の記録で、現在のGDCからは同じ形では
    参照できない。中身は2種類ある。
      胃がん    aliquot単位の Item flagged DNU (1,212行 = 11 aliquot)。
                DCCが「解析に使うな」とフラグを立てたもの
      膵臓がん  患者単位の study protocol 逸脱の注記 (10行)
    GDCはData Releaseのたびに配布内容を差し替えるため、
    「その日にどのファイルが存在し、何を取得したか」を残せるのはこれだけになる。
    現在の公式からは取得できない中期(HTSeq)のデータについては特に重要。
    取得日は inventory.json に記録される。

    15層(5がん種 x 3世代)のうち14層について収録している。欠けているのは
    乳がんの中期だけで、これはGDCから取得したものではない(UCSC Xena GDC hub
    由来)ため、GDC sample sheetがそもそも存在しない。

inventory.json      収録した行列の一覧(行数・列数・元ファイル名)
checksums.sha256    全ファイルのSHA-256
```

### sample_annotation の列

| 列 | 内容 |
|---|---|
| `patient_id` | TCGA-XX-XXXX |
| `sample_type_code` | 01=原発巣, 02=再発巣, 06=転移巣, 11=非がん部 など。**先頭ゼロ付きの文字列** |
| `group` | Tumor(原発巣) / TumorExtra(再発・転移) / Normal(非がん部) |
| `generations` | その検体が存在する世代 |
| `barcode_old` / `_mid` / `_new` | 各世代での実際のバーコード |
| `subtype_*` | 分類の割り当て(下記の出典を参照) |

読み込む際は `sample_type_code` を文字列として扱ってください。

```python
pd.read_csv("PAAD_samples.tsv", sep="\t", dtype={"sample_type_code": str})
```

---

## 検体の取捨選択

数字を並べる前に落としているものがあります。根拠はすべて `exclude/` と
`sample_annotation/` から検証できます。

- **転移巣・再発巣を原発巣と分けています。** 混ぜた分布を既定にすると、
  気づかないまま解釈されるためです。`group` 列で区別できます。
- **DCCがフラグを立てた検体は除外していません。** 胃がんでDNUフラグの付いた
  11 aliquot は3世代の行列すべてに残っており(旧世代11 / 中期9 / 新世代11検体)、
  膵臓がんで study protocol 逸脱の注記が付いた10名も3世代すべてに含まれています。
  公開図はこれらを含めた状態で作ったものです。除外して解析したい場合は
  `provenance/STAD/old/file_annotations.txt` と
  `provenance/PAAD/old/file_annotations.txt` を使ってください。

- **FFPE由来検体を除外しています。** バーコードのvial文字(01A/01B)では
  FFPEを判別できないため、GDCメタデータの `is_ffpe` /
  `Preservation Method` を根拠にしました。乳がんでは13検体が該当します
  (`exclude/` に全件を記録)。**うち `TCGA-PL-A8LZ-01A` は `01A` です。**
  vial文字を根拠にできない理由がそのまま出ている例です。
  実際に除外される数は世代によって違い、中期(1,217検体)では13検体すべて、
  新世代(1,208検体)では2検体、旧世代(1,154検体)では0検体です。
  残る11検体は `01B` / `01C` のvialで、新世代と旧世代の行列には
  同一検体の `01A` しか配布されていないため、もともと入っていません。
- **同一患者の重複vialを1検体に統合しています。** 統合しないと同じ患者が
  2回カウントされます。vial文字が若い方を採用しています。
- **膵臓がんでは28検体を除外しています。** 膵管腺癌でないと判定された19検体と、
  腫瘍細胞含有率が極端に低い9検体です(Raphael et al. 2017, Table S1)。

---

## 出典

発現データは以下から取得しています。**利用にあたっては各提供元の
利用規約に従ってください。**

- **旧TCGA Data Coordinating Center** (http://cancergenome.nih.gov/) —
  旧世代(RNASeqV2 normalized count)は5がん種すべてこのポータルから取得しました。
  現在は閉鎖され、TCGAデータはNIH GDCに移管されています。
- **NCI Genomic Data Commons** (https://portal.gdc.cancer.gov/) —
  中期(HTSeq FPKM / FPKM-UQ、乳がんを除く)と新世代(STAR-Counts TPM)。
- **UCSC Xena** (https://xenabrowser.net/) — 乳がん中期のみ。GDC hub が
  再配布している HTSeq FPKM (log2(FPKM+1)) を用いています。

分類(`subtype_*` 列)の出典は2種類あります。**これらを利用する場合は
それぞれの原著・提供元を引用してください。**

トリプルネガティブ(`subtype_TNBC`)は、TCGA BCR の biospecimen 臨床ファイル
(`nationwidechildrens.org_clinical_patient_brca.txt`) の ER / PR / HER2 の
判定から求めています。論文の補足資料ではありません。

以下は論文の補足資料に基づきます。

- Hoadley KA, et al. *Cell* 2018 — pan-cancer iCluster
- Liu Y, et al. *Cancer Cell* 2018 — 消化管腺癌の分子サブタイプ / MSI / CIMP / CMS
- Thorsson V, et al. *Immunity* 2018 — 免疫サブタイプ C1–C6、各がん種の既発表サブタイプ(PAM50を含む)
- Raphael BJ, et al. *Cancer Cell* 2017 — 膵管腺癌の Moffitt / Bailey / Collisson / 腫瘍純度

---

## 再現する

前処理スクリプトはGitHubリポジトリの `scripts/` にあります。
このパッケージの `matrix/` を入力にすれば、公開サイトと同じJSONが生成されます。

```bash
python3 scripts/preprocess_tcga.py \
    --cancer-type PAAD \
    --new-tpm       matrix/PAAD/PAAD_new_TPM.tsv.gz \
    --mid-fpkmuq    matrix/PAAD/PAAD_mid_FPKM_UQ.tsv.gz \
    --old-normcount matrix/PAAD/PAAD_old_normalized_count.tsv.gz \
    --exclude-samples exclude/PAAD_exclude_all.txt \
    --out-dir ./data
```

乳がんは中期のスケールが違うので、`--mid-fpkm` と `--mid-log2p1` の
**両方**を渡してください。`--mid-log2p1` を忘れると、log2値をFPKMとして
読んだことになります。この場合スクリプトはエラーで停止するので、
気づかないまま誤った図ができることはありません。

```bash
python3 scripts/preprocess_tcga.py \
    --cancer-type BRCA \
    --new-tpm       matrix/BRCA/BRCA_new_TPM.tsv.gz \
    --mid-fpkm      matrix/BRCA/BRCA_mid_FPKM_log2p1.tsv.gz \
    --mid-log2p1 \
    --old-normcount matrix/BRCA/BRCA_old_normalized_count.tsv.gz \
    --exclude-samples exclude/BRCA_ffpe_exclude_brca.txt \
    --out-dir ./data
```

手順の詳細は `PIPELINE.md` を参照してください。

## 整合性の確認

```bash
sha256sum -c checksums.sha256
```

## ライセンス

このパッケージの構成・注釈は CC-BY 4.0 です。
元となるTCGAデータおよび各論文の補足資料は、それぞれの提供元の条件に従います。
