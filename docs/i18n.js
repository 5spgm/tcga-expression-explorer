// ============================================================
// 表示文字列と、がん種コードの対訳表
// ============================================================
// app.js から参照する。文言を直したいときは基本ここだけ触ればよい。

// --- がん種コード -> 表示名 -------------------------------------------
// 略称(PAAD, LUAD…)はドライ解析の慣用表記で、専門外には通じにくい。
// 疾患名を主・略称を従にして「膵臓がん (PAAD)」の形で表示する。
// ここに無いコードは、コードそのものをそのまま表示する(害はない)。
const CANCER_LABELS = {
  ACC:  { ja: "副腎皮質がん",             en: "Adrenocortical carcinoma" },
  BLCA: { ja: "膀胱がん",                 en: "Bladder cancer" },
  BRCA: { ja: "乳がん",                   en: "Breast cancer" },
  CESC: { ja: "子宮頸がん",               en: "Cervical cancer" },
  CHOL: { ja: "胆管がん",                 en: "Cholangiocarcinoma" },
  COAD: { ja: "大腸がん(結腸)",           en: "Colon cancer" },
  DLBC: { ja: "びまん性大細胞型B細胞リンパ腫", en: "Diffuse large B-cell lymphoma" },
  ESCA: { ja: "食道がん",                 en: "Esophageal cancer" },
  GBM:  { ja: "膠芽腫",                   en: "Glioblastoma" },
  HNSC: { ja: "頭頸部がん",               en: "Head and neck cancer" },
  KICH: { ja: "腎嫌色素細胞がん",         en: "Kidney chromophobe" },
  KIRC: { ja: "腎細胞がん(淡明細胞型)",   en: "Kidney clear cell carcinoma" },
  KIRP: { ja: "腎細胞がん(乳頭状)",       en: "Kidney papillary cell carcinoma" },
  LAML: { ja: "急性骨髄性白血病",         en: "Acute myeloid leukemia" },
  LGG:  { ja: "低悪性度グリオーマ",       en: "Lower grade glioma" },
  LIHC: { ja: "肝細胞がん",               en: "Liver hepatocellular carcinoma" },
  LUAD: { ja: "肺腺がん",                 en: "Lung adenocarcinoma" },
  LUSC: { ja: "肺扁平上皮がん",           en: "Lung squamous cell carcinoma" },
  MESO: { ja: "中皮腫",                   en: "Mesothelioma" },
  OV:   { ja: "卵巣がん",                 en: "Ovarian cancer" },
  PAAD: { ja: "膵臓がん",                 en: "Pancreatic cancer" },
  PCPG: { ja: "褐色細胞腫・傍神経節腫",   en: "Pheochromocytoma and paraganglioma" },
  PRAD: { ja: "前立腺がん",               en: "Prostate cancer" },
  READ: { ja: "大腸がん(直腸)",           en: "Rectal cancer" },
  SARC: { ja: "肉腫",                     en: "Sarcoma" },
  SKCM: { ja: "皮膚メラノーマ",           en: "Skin melanoma" },
  STAD: { ja: "胃がん",                   en: "Stomach cancer" },
  TGCT: { ja: "精巣腫瘍",                 en: "Testicular germ cell tumor" },
  THCA: { ja: "甲状腺がん",               en: "Thyroid cancer" },
  THYM: { ja: "胸腺腫",                   en: "Thymoma" },
  UCEC: { ja: "子宮体がん",               en: "Uterine corpus endometrial carcinoma" },
  UCS:  { ja: "子宮がん肉腫",             en: "Uterine carcinosarcoma" },
  UVM:  { ja: "ぶどう膜メラノーマ",       en: "Uveal melanoma" },
};

function cancerLabel(code, lang) {
  const entry = CANCER_LABELS[code];
  return entry ? `${entry[lang]} (${code})` : code;
}

// --- UI文言 ------------------------------------------------------------
const STRINGS = {
  ja: {
    htmlLang: "ja",
    pageTitle: "TCGA Expression Explorer — 新旧パイプライン比較",
    subtitle: "同一検体・同一遺伝子を、3世代の定量パイプラインで並べて見る",
    controlsLabel: "検索条件",
    cancerLabel: "がん種",
    cancerPlaceholder: "がん種を選択…",
    cancerLoading: "読み込み中…",
    geneLabel: "遺伝子",
    genePlaceholder: "遺伝子名を入力 (例: TP53)",
    subtypeLabel: "腫瘍の分類",
    subtypeNone: "なし(腫瘍 / 非がん部)",
    extraTumor: "転移巣・再発巣を含む",
    extraTumorWithCount: (n) => `転移巣・再発巣を含む (${n}検体)`,
    extraTumorNone: "このがん種には転移巣・再発巣の検体がありません",
    extraTumorHint: (detail) => `既定では原発巣のみを表示します。内訳: ${detail}`,
    logScale: "log2(x+1)表示",
    emptyState: "がん種と遺伝子を選ぶと、腫瘍部と非がん部の発現分布が3世代分並んで表示されます。",
    genOld: "旧世代",
    genMid: "中期",
    genNew: "新世代",
    footerSource: "データソース: GDC / TCGA public RNA-seq matrices",

    valueTypeTitle: "値の種類について",
    clear: "クリア",
    clearTitle: "がん種・遺伝子・分類の指定をすべて消して最初の状態に戻します",
    changelogTitle: "更新履歴",
    citeTitle: "引用",
    citeBody: "本ツールを研究に利用された場合は、下記リポジトリをご参照ください。" +
              "分類(サブタイプ)を利用した場合は、下に挙げた原著論文も併せて引用してください。",
    dataSourceTitle: "データの出典",
    dataSourceBody: "発現データは NCI Genomic Data Commons、UCSC Xena、cBioPortal から取得しています。" +
                    "分類の割り当ては以下の論文の補足資料に基づきます。",
    dataSourceNote: "本ツールは研究目的で提供しています。診断や治療の判断には使用できません。" +
                    "TCGAデータの利用にあたっては、各提供元の利用規約に従ってください。",
    licenseNote: "コードは MIT License。表示内容の正確性は保証されません。",

    tumor: "腫瘍",
    normalTissue: "非がん部",
    normalTissueAlt: "非がん部(正常組織)",

    statusReady: (n) => `${n} がん種が利用可能です`,
    statusManifestError: (msg) => `manifest.json の読み込みに失敗しました (${msg})`,
    statusGeneListLoading: "遺伝子リストを読み込み中…",
    statusGeneListError: (msg) => `遺伝子リストの読み込みに失敗しました (${msg})`,
    statusGeneList: (type, n) => `${type}: ${n.toLocaleString()} 遺伝子が利用可能です`,
    statusGeneLoading: (g) => `${g} を読み込み中…`,
    statusGeneError: (msg) => `遺伝子データの読み込みに失敗しました (${msg})`,
    statusGeneNotFound: (g) => `"${g}" は見つかりませんでした。候補一覧から選んでください。`,
    statusShowing: (type, sym, ens) => `${type} / ${sym} (${ens}) を表示中`,
    statusGeneCleared: "遺伝子名を入力してください。",
    statusGeneNotInCancer: (g, type) =>
      `"${g}" は ${type} のデータには含まれていません。別の遺伝子を入力してください。`,

    plotNoData: "このがん種では未取得のデータです",
    plotNoGene: "この遺伝子は対応するIDが見つかりませんでした",

    footNormal: (n) => `非がん部 n=${n}`,
    footTumor: (n) => `腫瘍 n=${n}`,
    footExtraIncluded: (n) => ` [転移巣・再発巣 ${n}件を含む]`,
    footExtraExcluded: (n) => ` [転移巣・再発巣 ${n}件を除外]`,
    smallNWarning: (groups) =>
      ` ⚠ 検体数が少ない群があります(${groups})。箱ひげ図の四分位数は不安定です。`,

    download: "保存",
    downloadPng: "PNG",
    downloadSvg: "SVG",
    downloadCsv: "CSV",
    downloadPngTitle: "画像として保存(発表資料向け・高解像度)",
    downloadSvgTitle: "ベクター形式で保存(論文図版向け)",
    downloadCsvTitle: "描画中の数値をCSVで保存",
    csvHeaderGroup: "群",
    csvHeaderValue: "値",
  },

  en: {
    htmlLang: "en",
    pageTitle: "TCGA Expression Explorer — pipeline generation comparison",
    subtitle: "The same samples and gene, side by side across three quantification pipelines",
    controlsLabel: "Query",
    cancerLabel: "Cancer type",
    cancerPlaceholder: "Select a cancer type…",
    cancerLoading: "Loading…",
    geneLabel: "Gene",
    genePlaceholder: "Enter a gene symbol (e.g. TP53)",
    subtypeLabel: "Tumor grouping",
    subtypeNone: "None (Tumor / Normal tissue)",
    extraTumor: "Include metastatic / recurrent",
    extraTumorWithCount: (n) => `Include metastatic / recurrent (${n} samples)`,
    extraTumorNone: "No metastatic or recurrent samples for this cancer type",
    extraTumorHint: (detail) => `Primary tumors only by default. Breakdown: ${detail}`,
    logScale: "log2(x+1) scale",
    emptyState: "Choose a cancer type and a gene to see tumor and normal-tissue expression across three pipeline generations.",
    genOld: "Old",
    genMid: "Mid",
    genNew: "New",
    footerSource: "Data source: GDC / TCGA public RNA-seq matrices",

    valueTypeTitle: "Units of the plotted values",
    clear: "Clear",
    clearTitle: "Reset the cancer type, gene and grouping to the initial state",
    changelogTitle: "Changelog",
    citeTitle: "Citation",
    citeBody: "If you use this tool in your research, please refer to the repository below. " +
              "If you use the molecular subtype assignments, please also cite the original papers listed here.",
    dataSourceTitle: "Data sources",
    dataSourceBody: "Expression data were obtained from the NCI Genomic Data Commons, UCSC Xena and cBioPortal. " +
                    "Subtype assignments are derived from the supplementary material of the following papers.",
    dataSourceNote: "This tool is provided for research purposes only and must not be used for diagnosis or treatment decisions. " +
                    "Use of TCGA data is subject to the terms of the respective data providers.",
    licenseNote: "Code is released under the MIT License. No warranty is given as to the accuracy of the displayed content.",

    tumor: "Tumor",
    normalTissue: "Normal tissue",
    normalTissueAlt: "Normal tissue (non-cancerous)",

    statusReady: (n) => `${n} cancer types available`,
    statusManifestError: (msg) => `Failed to load manifest.json (${msg})`,
    statusGeneListLoading: "Loading gene list…",
    statusGeneListError: (msg) => `Failed to load the gene list (${msg})`,
    statusGeneList: (type, n) => `${type}: ${n.toLocaleString()} genes available`,
    statusGeneLoading: (g) => `Loading ${g}…`,
    statusGeneError: (msg) => `Failed to load gene data (${msg})`,
    statusGeneNotFound: (g) => `"${g}" was not found. Please pick one from the suggestions.`,
    statusShowing: (type, sym, ens) => `Showing ${type} / ${sym} (${ens})`,
    statusGeneCleared: "Enter a gene symbol.",
    statusGeneNotInCancer: (g, type) =>
      `"${g}" is not present in the ${type} dataset. Please try another gene.`,

    plotNoData: "Not collected for this cancer type",
    plotNoGene: "No matching gene ID in this dataset",

    footNormal: (n) => `Normal tissue n=${n}`,
    footTumor: (n) => `Tumor n=${n}`,
    footExtraIncluded: (n) => ` [including ${n} metastatic/recurrent]`,
    footExtraExcluded: (n) => ` [excluding ${n} metastatic/recurrent]`,
    smallNWarning: (groups) =>
      ` ⚠ Small group size (${groups}). Box quartiles are unstable.`,

    download: "Save",
    downloadPng: "PNG",
    downloadSvg: "SVG",
    downloadCsv: "CSV",
    downloadPngTitle: "Save as an image (high resolution, for slides)",
    downloadSvgTitle: "Save as vector graphics (for publication figures)",
    downloadCsvTitle: "Save the plotted values as CSV",
    csvHeaderGroup: "group",
    csvHeaderValue: "value",
  },
};


// --- 引用すべき原著 ---------------------------------------------------
// 分類(サブタイプ)の割り当ては、いずれもこれらの論文の補足資料に由来する。
// 表示は言語によらず英語表記のまま(書誌情報のため)。
const REFERENCES = [
  { text: "Hoadley KA, et al. Cell 2018 — pan-cancer iCluster",
    url: "https://doi.org/10.1016/j.cell.2018.03.022" },
  { text: "Liu Y, et al. Cancer Cell 2018 — gastrointestinal adenocarcinoma subtypes, MSI, CIMP, CMS",
    url: "https://doi.org/10.1016/j.ccell.2018.03.010" },
  { text: "Thorsson V, et al. Immunity 2018 — immune subtypes C1–C6, published per-cancer subtypes incl. PAM50",
    url: "https://doi.org/10.1016/j.immuni.2018.03.023" },
  { text: "Raphael BJ, et al. Cancer Cell 2017 — pancreatic ductal adenocarcinoma (Moffitt / Bailey / Collisson / purity)",
    url: "https://doi.org/10.1016/j.ccell.2017.07.007" },
];

// --- 値の種類の説明 ---------------------------------------------------
// 「normalized count とは具体的に何か」が分からないという指摘への対応。
// パネル上部のバッジにツールチップとして出し、フッターにも一覧を置く。
//
// 要点は「遺伝子長で割っているかどうか」。
//   normalized_count は遺伝子長で割っていないため、**同じ遺伝子を検体間で
//   比べる**のには使えるが、1検体の中で別の遺伝子どうしを比べるのには適さない。
//   TPM / FPKM は遺伝子長で割っているので、その比較にも使える。
const VALUE_TYPE_INFO = {
  normalized_count: {
    ja: "RSEMの推定カウントを、その検体の上位四分位で割って1000倍した値" +
        "(TCGA RNASeqV2 の normalized_results)。遺伝子長で割っていないため、" +
        "同じ遺伝子を検体間で比べる用途に向く。1つの検体の中で別の遺伝子と" +
        "大小を比べるのには適さない。",
    en: "RSEM expected counts divided by the sample's upper quartile and scaled by 1000 " +
        "(TCGA RNASeqV2 normalized_results). Not length-normalised, so it is suited to " +
        "comparing the same gene across samples, but not to comparing different genes " +
        "within one sample.",
  },
  FPKM: {
    ja: "Fragments Per Kilobase of transcript per Million mapped reads。" +
        "遺伝子長と総リード数で割った値。",
    en: "Fragments Per Kilobase of transcript per Million mapped reads; " +
        "normalised by gene length and library size.",
  },
  FPKM_UQ: {
    ja: "FPKMの分母を総リード数ではなく上位四分位に置き換えたもの。" +
        "極端に発現の高い遺伝子の影響を受けにくい。",
    en: "FPKM with the library-size denominator replaced by the upper quartile, " +
        "reducing the influence of very highly expressed genes.",
  },
  TPM: {
    ja: "Transcripts Per Million。遺伝子長で割ってから総和が100万になるよう" +
        "揃えた値。検体間・遺伝子間のどちらの比較にも使える。",
    en: "Transcripts Per Million: length-normalised first, then scaled so that each " +
        "sample sums to one million. Comparable both across samples and across genes.",
  },
};

function valueTypeInfo(vt, lang) {
  const e = VALUE_TYPE_INFO[vt];
  return e ? e[lang] || e.en : "";
}


