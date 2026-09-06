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
    ja: "旧TCGA DCCのData Matrix (https://tcga-data.nci.nih.gov/)。" +
        "2016年にNCI Genomic Data Commonsへ移行し、現在は閉鎖。" +
        "当時は一度にダウンロードできる容量に上限があり、" +
        "複数の端末に分けて取得した。ダウンロード時には " +
        "file_manifest.txt / FILE_SAMPLE_MAP.txt / file_annotations.txt が" +
        "同梱され、後者にはDCCが解析対象から除外するよう指示した検体" +
        "(Item flagged DNU)の記録が含まれる。",
    en: "The Data Matrix of the former TCGA Data Coordinating Center " +
        "(https://tcga-data.nci.nih.gov/), migrated to the NCI Genomic Data " +
        "Commons in 2016 and now retired. The portal capped the size of a " +
        "single download at the time, so files were retrieved across several " +
        "machines. Each download included file_manifest.txt, FILE_SAMPLE_MAP.txt " +
        "and file_annotations.txt; the last of these records aliquots the DCC " +
        "flagged for exclusion from analysis (Item flagged DNU).",
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
      source: "gdc_htseq", acquired: "2021" },
    { generation: "new", pipeline: "GDC STAR-Counts", value: "TPM",
      source: "gdc_star", acquired: "2025-06",
      ja: "FFPE由来の13検体を除外(うち1検体はvial文字が01A)。",
      en: "13 FFPE-derived samples excluded; one of them carried the 01A vial letter." },
  ],
  COAD: [
    { generation: "old", pipeline: "UNC IlluminaHiSeq RNASeqV2", value: "normalized_count",
      source: "tcga_portal", acquired: "2014–2015" },
    { generation: "mid", pipeline: "GDC HTSeq", value: "FPKM_UQ",
      source: "gdc_htseq", acquired: "2021" },
    { generation: "new", pipeline: "GDC STAR-Counts", value: "TPM",
      source: "gdc_star", acquired: "2022-12" },
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
      source: "gdc_star", acquired: "2026-09" },
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
