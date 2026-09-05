// ============================================================
// 更新履歴
// ============================================================
// 新しいものを先頭に足していく。日付は YYYY-MM-DD。
// ja / en の両方を書くこと(片方だけだと言語切替で空欄になる)。
//
// version は任意。Zenodoに寄託したデータのバージョンと対応させておくと、
// 「どの図がどのデータで作られたか」を後から追える。

const CHANGELOG = [
  {
    date: "2026-09-05",
    ja: [
      "胃がん (STAD) を追加(分子サブタイプ / MSI / CIMP / 免疫サブタイプ / iCluster)",
      "指定をまとめて消す「クリア」ボタンを追加しました",
    ],
    en: [
      "Added gastric cancer (STAD) with molecular subtype, MSI, CIMP, immune subtype and iCluster",
      "Added a Clear button that resets all selections at once",
    ],
  },
  {
    date: "2026-08-29",
    ja: [
      "遺伝子名が空欄のときに、前に表示していた遺伝子の図が残ることがある問題を修正しました",
      "各世代の値の種類(normalized count / FPKM / FPKM-UQ / TPM)の説明を追加しました",
      "引用・データの出典・ライセンスをページ下部に明記しました",
      "更新履歴の表示を追加しました",
    ],
    en: [
      "Fixed a case where the previously displayed gene remained plotted after the gene field was cleared",
      "Added explanations of the value types (normalized count / FPKM / FPKM-UQ / TPM)",
      "Added citation, data source and licence information in the footer",
      "Added this changelog",
    ],
  },
  {
    date: "2026-08-27",
    version: "data v1.0",
    ja: [
      "肺腺がん (LUAD) を追加(3世代 / 免疫サブタイプ・iCluster)",
      "がん種を切り替えても入力した遺伝子名を引き継ぐようにしました",
      "言語を切り替えたときに状態表示が追従しない問題を修正しました",
    ],
    en: [
      "Added lung adenocarcinoma (LUAD) with all three generations",
      "The gene name is now kept when switching between cancer types",
      "Fixed the status line not following the language switch",
    ],
  },
  {
    date: "2026-08-25",
    ja: [
      "膵臓がん (PAAD) を追加(Moffitt / Bailey / 腫瘍純度 / 免疫サブタイプ)",
      "3世代でグラフの箱の並び順を揃えました",
      "英数字の書体を Arial にしました(図の書き出しにも反映されます)",
      "グラフを PNG / SVG / CSV で保存できるようにしました",
    ],
    en: [
      "Added pancreatic cancer (PAAD) with Moffitt, Bailey, tumour purity and immune subtypes",
      "Box order is now consistent across the three generations",
      "Latin text is now set in Arial, including exported figures",
      "Plots can be saved as PNG, SVG or CSV",
    ],
  },
  {
    date: "2026-08-23",
    ja: [
      "乳がん (BRCA) を追加(PAM50 / TNBC / 免疫サブタイプ / iCluster)",
      "転移巣・再発巣を既定で除外し、チェックボックスで合流できるようにしました",
      "FFPE由来検体を除外しました",
      "日本語・英語の表示切り替えを追加しました",
    ],
    en: [
      "Added breast cancer (BRCA) with PAM50, TNBC, immune subtypes and iCluster",
      "Metastatic and recurrent samples are now excluded by default, with a checkbox to include them",
      "FFPE-derived samples are excluded",
      "Added Japanese / English language switching",
    ],
  },
  {
    date: "2026-08-22",
    ja: ["公開(大腸がん / COAD)"],
    en: ["First release (colon cancer, COAD)"],
  },
];
