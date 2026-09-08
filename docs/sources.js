// ============================================================
// データの出典(がん種ごと)
// ============================================================
// sources.html から読み込む。がん種を追加したらここに1件足す。
//
// acquired は「取得した時期」。当時のポータルは既に閉鎖されている場合が
// あるため、いつ・どこから取ったかを残しておくことが再現性の担保になる。

// --- 取得元の定義 -----------------------------------------------------
const SOURCE_DEFS = {
  tcga_portal: {
    name: "TCGA Data Coordinating Center (DCC) — Data Matrix",
    url: null,                    // 閉鎖済みのためリンクは張らない
    retired: true,
    ja: "旧TCGAポータル (http://cancergenome.nih.gov) のData Matrix " +
        "(https://tcga-data.nci.nih.gov/)。本ツールの旧世代は5がん種すべてを" +
        "この経路から取得している。2016年にNCI Genomic Data Commonsへ移行し、現在は閉鎖。" +
        "当時は一度にダウンロードできる容量に上限があり、" +
        "複数の端末に分けて取得した。ダウンロード時には " +
        "file_manifest.txt / FILE_SAMPLE_MAP.txt / file_annotations.txt が" +
        "同梱され、後者にはDCCが解析対象から除外するよう指示した検体" +
        "(Item flagged DNU)の記録が含まれる。" +
        "この経路の行列は遺伝子行が SYMBOL|EntrezID 形式(例: TP53|7157)であり、" +
        "遺伝子記号のみを用いる他の再配布元と区別できる。",
    en: "The Data Matrix of the former TCGA portal (http://cancergenome.nih.gov) " +
        "at https://tcga-data.nci.nih.gov/. All five old-generation cohorts in this " +
        "tool were obtained by this route. It migrated to the NCI Genomic Data " +
        "Commons in 2016 and is now retired. The portal capped the size of a " +
        "single download at the time, so files were retrieved across several " +
        "machines. Each download included file_manifest.txt, FILE_SAMPLE_MAP.txt " +
        "and file_annotations.txt; the last of these records aliquots the DCC " +
        "flagged for exclusion from analysis (Item flagged DNU). The matrices obtained " +
        "by this route carry gene rows in SYMBOL|EntrezID form (e.g. TP53|7157), which " +
        "distinguishes them from redistributions that use gene symbols alone.",
  },
  gdc_htseq: {
    name: "NCI Genomic Data Commons (HTSeq era)",
    url: "https://portal.gdc.cancer.gov/",
    retired: false,
    ja: "GDCのData Release 15〜31の期間に取得したHTSeq由来のFPKM / FPKM-UQ。" +
        "2022年のDR32でSTAR-Countsへ置き換えられたため、" +
        "現在のGDCからは同じものを取得できない。",
    en: "FPKM / FPKM-UQ from the HTSeq pipeline, downloaded while GDC Data " +
        "Releases 15–31 were current. Replaced by STAR-Counts at DR32 in 2022, " +
        "so the same files can no longer be obtained from the GDC.",
  },
  xena_gdc_hub: {
    name: "UCSC Xena — GDC hub (GDC-PANCAN.htseq_fpkm)",
    url: "https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net",
    retired: false,
    ja: "GDCのHTSeq期のFPKMを、UCSC XenaのGDC hubが再配布しているもの。" +
        "GDCはコホート単位の発現行列を配布せず、HTSeq期の値は2022年のDR32以降" +
        "取得できないため、乳がんの中期のみこの経路で取得した。" +
        "取得日は2024-04-13で、当時は TCGA-BRCA.htseq_fpkm.tsv.gz として" +
        "配布されていた(手元の .gz のタイムスタンプがこの日時。gzipヘッダには" +
        "配布側のファイル時刻 2019-07-19 が残っている)。" +
        "その後hubのカタログはSTAR期(gencode v36)に入れ替わり、" +
        "コホート単位のHTSeq FPKMは削除された。2026-09-08時点で残っている" +
        "HTSeq FPKMは pan-cancerの GDC-PANCAN.htseq_fpkm (gencode v22) だけで、" +
        "本サイトの1,217検体はすべてこれに含まれ、値も一致する。" +
        "Xenaはこのhubの発現値を log2(x+1) に変換して配布しているので、" +
        "読み込み時に 2^x - 1 で線形スケールへ戻している。" +
        "値の種類はFPKMで、他がん種の中期(FPKM-UQ)とは異なる。",
    en: "HTSeq-era FPKM from the GDC, redistributed by the GDC hub of UCSC Xena. " +
        "The GDC does not distribute cohort-level expression matrices and the " +
        "HTSeq-era values became unavailable at DR32 in 2022, so this route was " +
        "used for the middle generation of the breast cohort only. It was " +
        "downloaded on 2024-04-13, when it was distributed as " +
        "TCGA-BRCA.htseq_fpkm.tsv.gz (the timestamp of the local .gz; the gzip " +
        "header still carries the distribution-side file time of 2019-07-19). " +
        "The hub catalogue has since moved " +
        "to the STAR era (gencode v36) and the cohort-level HTSeq FPKM datasets " +
        "have been removed. As of 2026-09-08 the only HTSeq FPKM left is the " +
        "pan-cancer GDC-PANCAN.htseq_fpkm (gencode v22); all 1,217 samples used " +
        "here are present in it and the values agree. Xena distributes " +
        "the values of this hub as log2(x+1); they are back-transformed with " +
        "2^x - 1 on load. The quantity is FPKM, unlike the FPKM-UQ used for the " +
        "other cohorts.",
  },
  gdc_star: {
    name: "NCI Genomic Data Commons (STAR-Counts)",
    url: "https://portal.gdc.cancer.gov/",
    retired: false,
    ja: "現行のGDCから取得したSTAR-Countsのgene_counts。" +
        "同梱のスクリプトで検体ごとのファイルを1つの行列にまとめている。",
    en: "STAR-Counts gene-level quantifications from the current GDC. " +
        "The per-sample files are assembled into a single matrix by the " +
        "scripts included in the repository.",
  },
};

// --- がん種ごとの取得状況 ----------------------------------------------
// generation: old / mid / new
//   pipeline  実際に走ったパイプライン
//   value     値の種類(サイト上の表示と一致させる)
//   source    SOURCE_DEFS のキー
//   acquired  取得時期
//   note      そのがん種固有の事情(任意)
const DATA_SOURCES = {
  BRCA: [
    { generation: "old", pipeline: "UNC IlluminaHiSeq RNASeqV2", value: "normalized_count",
      source: "tcga_portal", acquired: "2014–2015" },
    { generation: "mid", pipeline: "GDC HTSeq", value: "FPKM",
      source: "xena_gdc_hub", acquired: "2024-04",
      ja: "5がん種の中で唯一、GDCから直接取得できなかった層。値はFPKM(他は FPKM-UQ)で、" +
          "log2(x+1)で配布されているため読み込み時に 2^x - 1 で線形へ戻している。",
      en: "The only layer of the five cohorts that could not be obtained from the GDC " +
          "directly. The quantity is FPKM (FPKM-UQ elsewhere) and is distributed as " +
          "log2(x+1), so it is back-transformed with 2^x - 1 on load." },
    { generation: "new", pipeline: "GDC STAR-Counts", value: "TPM",
      source: "gdc_star", acquired: "2025-06",
      ja: "FFPEと判定された13検体のうち、同一患者×同一sample typeの別vialを持たない" +
          "2検体が解析対象から外れた。残り11検体は重複vialの解消の段階で既に除去されている。",
      en: "Of the 13 samples identified as FFPE, only 2 were actually dropped from the " +
          "analysis, having no alternative vial of the same patient and sample type; the " +
          "remaining 11 had already been removed at the vial-deduplication step." },
  ],
  COAD: [
    { generation: "old", pipeline: "UNC IlluminaHiSeq RNASeqV2", value: "normalized_count",
      source: "tcga_portal", acquired: "2014–2015" },
    { generation: "mid", pipeline: "GDC HTSeq", value: "FPKM_UQ",
      source: "gdc_htseq", acquired: "2021" },
    { generation: "new", pipeline: "GDC STAR-Counts", value: "TPM",
      source: "gdc_star", acquired: "2022-12",
      ja: "取得時のシートには頭頸部がん548検体が混在しており(大腸524検体)、" +
          "Project ID列で選別している。",
      en: "The sample sheet used at the time also listed 548 head-and-neck samples " +
          "alongside the 524 colon samples; they are separated by the Project ID column." },
  ],
  PAAD: [
    { generation: "old", pipeline: "UNC IlluminaHiSeq RNASeqV2", value: "normalized_count",
      source: "tcga_portal", acquired: "2015",
      ja: "当時の登録は44検体のみ(腫瘍43 / 非がん部1)。",
      en: "Only 44 samples were available at the time (43 tumour, 1 normal)." },
    { generation: "mid", pipeline: "GDC HTSeq", value: "FPKM_UQ",
      source: "gdc_htseq", acquired: "2021" },
    { generation: "new", pipeline: "GDC STAR-Counts", value: "TPM",
      source: "gdc_star", acquired: "2023-02",
      ja: "膵管腺癌でない19検体と、腫瘍細胞含有率が極端に低い9検体を除外。",
      en: "19 samples not classified as pancreatic ductal adenocarcinoma and 9 with very low neoplastic cellularity were excluded." },
  ],
  LUAD: [
    { generation: "old", pipeline: "UNC IlluminaHiSeq RNASeqV2", value: "normalized_count",
      source: "tcga_portal", acquired: "2014–2015" },
    { generation: "mid", pipeline: "GDC HTSeq", value: "FPKM_UQ",
      source: "gdc_htseq", acquired: "2021" },
    { generation: "new", pipeline: "GDC STAR-Counts", value: "TPM",
      source: "gdc_star", acquired: "2026-08" },
  ],
  STAD: [
    { generation: "old", pipeline: "UNC IlluminaHiSeq RNASeqV2", value: "normalized_count",
      source: "tcga_portal", acquired: "2015-07以降",
      ja: "同梱の file_annotations.txt の最新記録が2015-07-08であることから、" +
          "取得はそれ以降。450検体・Level 3。",
      en: "The most recent entry in the bundled file_annotations.txt is dated " +
          "2015-07-08, so the download was made after that date. 450 samples, Level 3." },
    { generation: "mid", pipeline: "GDC HTSeq", value: "FPKM_UQ",
      source: "gdc_htseq", acquired: "2021-02-05" },
    { generation: "new", pipeline: "GDC STAR-Counts", value: "TPM",
      source: "gdc_star", acquired: "2024-10",
      ja: "検体ごとのファイルは2024-10-14に取得し、行列は2024-10-18に作成した。" +
          "2026-09-05に再出力した取得シートも、検体集合とFile Name集合が一致する。",
      en: "The per-sample files were downloaded on 2024-10-14 and the matrix was built " +
          "on 2024-10-18. A sample sheet re-exported on 2026-09-05 lists an identical " +
          "set of samples and file names." },
  ],
};

// --- 分類(サブタイプ)の出典 --------------------------------------------
// どのがん種でどれが使えるかも併記する。
const SUBTYPE_SOURCES = [
  { label: "iCluster (Hoadley 2018)",
    citation: "Hoadley KA, et al. Cell 2018;173(2):291-304.e6",
    url: "https://doi.org/10.1016/j.cell.2018.03.022",
    used_in: ["BRCA", "COAD", "LUAD", "STAD"] },
  { label: "Molecular Subtype / MSI / CIMP / CMS (Liu 2018)",
    citation: "Liu Y, et al. Cancer Cell 2018;33(4):721-735.e8",
    url: "https://doi.org/10.1016/j.ccell.2018.03.010",
    used_in: ["COAD", "STAD"] },
  { label: "Immune Subtype C1–C6 / PAM50 (Thorsson 2018)",
    citation: "Thorsson V, et al. Immunity 2018;48(4):812-830.e14",
    url: "https://doi.org/10.1016/j.immuni.2018.03.023",
    used_in: ["BRCA", "COAD", "PAAD", "LUAD", "STAD"] },
  { label: "Moffitt / Bailey / Collisson / tumour purity (Raphael 2017)",
    citation: "Raphael BJ, et al. Cancer Cell 2017;32(2):185-203.e13",
    url: "https://doi.org/10.1016/j.ccell.2017.07.007",
    used_in: ["PAAD"] },
  { label: "TNBC (ER/PR/HER2 by IHC)",
    citation: "TCGA Biospecimen Core Resource, clinical patient file",
    url: null,
    used_in: ["BRCA"] },
];
