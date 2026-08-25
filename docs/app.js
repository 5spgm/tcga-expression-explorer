// ============================================================
// TCGA Expression Explorer
// ============================================================
// データの置き場所。自前サーバーにデータを置く場合はここを絶対URLに変更する。
// 末尾にスラッシュを付けないこと(下の fetch 側が `${DATA_BASE_URL}/...` の形で組み立てる)。
// (自前サーバー側で Access-Control-Allow-Origin ヘッダの設定と、HTTPS化が必要。
//  詳細はREADMEの「3. データを自前サーバーに置く場合」を参照)
const DATA_BASE_URL = "https://importantly-ministers-inquiry-through.trycloudflare.com/data";

// 非がん部の箱に使う名前。単に "Normal" にすると、PAM50の "Normal"
// (Normal-like)のようにsubtype側に同名のラベルがある場合、Plotlyが同じ
// カテゴリとして束ねてしまい、2つの箱が重なって描画される。
// 衝突しない名前を使い、必ず右端に置く。
// グラフ内のフォント。英数字はArial、日本語は和文フォントへ自動で落ちる。
// PNG/SVGで書き出した図にもこの指定がそのまま乗る。
const PLOT_FONT = 'Arial, Helvetica, "Hiragino Sans", "Yu Gothic", "Noto Sans JP", sans-serif';

function normalTraceName(boxes) {
  // subtype側に同名ラベル(PAM50の "Normal" 等)があると Plotly が同じ
  // カテゴリとして束ね、箱が重なってしまう。衝突する場合だけ名前をずらす。
  const base = t.normalTissue;
  return boxes && boxes.some((b) => b.name === base) ? t.normalTissueAlt : base;
}

const COLORS = {
  tumor: "#c1440e",
  tumorSoft: "rgba(193, 68, 14, 0.35)",
  normal: "#2e7d6b",
  normalSoft: "rgba(46, 125, 107, 0.35)",
  ink: "#1b231f",
  rule: "#d8dad3",
};

// Tumorをsubtypeで分けて表示するときの配色(カテゴリカル、6色を使い回す)
const SUBTYPE_PALETTE = [
  "#c1440e", "#8a5fc9", "#c98a2b", "#b2456e", "#5c7a99", "#8a8f5c",
];

// この人数を下回る群には脚注で注意を出す(膵臓がんの正常組織は4-5検体しかない)
const SMALL_N_THRESHOLD = 10;

const GENERATIONS = [
  { key: "old", plotId: "plot-old", rowId: "old-value-type-row", footId: "footnote-old", dlId: "download-old" },
  { key: "mid", plotId: "plot-mid", rowId: "mid-value-type-row", footId: "footnote-mid", dlId: "download-mid" },
  { key: "new", plotId: "plot-new", rowId: "new-value-type-row", footId: "footnote-new", dlId: "download-new" },
];

const state = {
  manifest: null,
  cancerType: null,
  geneIndex: [],          // [{ensembl, symbol}]
  symbolToEnsembl: new Map(),
  currentGeneJson: null,
  activeValueType: { old: null, mid: null, new: null },
  logScale: true,          // log2(x+1) 表示がデフォルト(裾の長い分布で箱が潰れるのを防ぐ)
  subtypeScheme: "__none__", // Tumorの分類方式。"__none__" ならTumor/Normalの2群のみ
  includeExtraTumor: false,  // 転移巣・再発巣(sample type code 02/06など)を含めるか。既定は除外。
  lang: "ja",                // "ja" | "en"
  boxOrder: null,            // 3世代で共通の箱の並び順(renderAllPanelsで決める)
};

// 現在の言語の文字列辞書。t.xxx で参照する。
let t = STRINGS.ja;

const LANG_KEY = "tcga-explorer-lang";

function detectLang() {
  try {
    const saved = localStorage.getItem(LANG_KEY);
    if (saved && STRINGS[saved]) return saved;
  } catch (e) { /* プライベートモード等でlocalStorageが使えない場合は無視 */ }
  return (navigator.language || "en").toLowerCase().startsWith("ja") ? "ja" : "en";
}

function setLang(lang) {
  if (!STRINGS[lang]) return;
  state.lang = lang;
  t = STRINGS[lang];
  try { localStorage.setItem(LANG_KEY, lang); } catch (e) { /* 保存できなくても動作に支障はない */ }
  applyStaticStrings();
  refreshDynamicStrings();
}

// data-i18n 属性を持つ要素に文字列を流し込む
function applyStaticStrings() {
  document.documentElement.lang = t.htmlLang;
  document.title = t.pageTitle;
  document.querySelectorAll("[data-i18n]").forEach((node) => {
    const value = t[node.dataset.i18n];
    if (typeof value === "string") node.textContent = value;
  });
  el("gene-input").placeholder = t.genePlaceholder;
  document.querySelectorAll(".lang-switch button").forEach((b) => {
    b.classList.toggle("active", b.dataset.lang === state.lang);
  });
}

// 言語切替時に、既に描画済みの動的部分を作り直す
function refreshDynamicStrings() {
  const cancerSelect = el("cancer-select");
  if (state.manifest) {
    const selected = cancerSelect.value;
    fillCancerSelect();
    cancerSelect.value = selected;
    if (!selected) setStatus(t.statusReady(Object.keys(state.manifest.cancer_types || {}).length));
  }
  const subtypeSelect = el("subtype-scheme-select");
  if (subtypeSelect && subtypeSelect.options.length) {
    const first = subtypeSelect.options[0];
    if (first && first.value === "__none__") first.textContent = t.subtypeNone;
  }
  if (state.cancerType) updateExtraTumorControl(state.cancerType);
  if (state.currentGeneJson) {
    renderAllPanels();
    setStatus(t.statusShowing(cancerLabel(state.cancerType, state.lang),
                              state.currentGeneJson.gene_symbol,
                              state.currentGeneJson.gene_id));
  } else {
    hidePanels();
  }
}

function fillCancerSelect() {
  const cancerSelect = el("cancer-select");
  const types = Object.keys(state.manifest.cancer_types || {});
  // 表示名(疾患名)の五十音/アルファベット順に並べる。コード順だと
  // 「BRCA, COAD, PAAD」のようになり、専門外には探しにくい。
  types.sort((a, b) => cancerLabel(a, state.lang).localeCompare(cancerLabel(b, state.lang), state.lang));
  cancerSelect.innerHTML =
    `<option value="">${t.cancerPlaceholder}</option>` +
    types.map((c) => `<option value="${c}">${cancerLabel(c, state.lang)}</option>`).join("");
}

const el = (id) => document.getElementById(id);
const statusLine = el("status-line");
const emptyState = el("empty-state");
const panelsSection = el("panels");

function setStatus(msg, isError = false) {
  statusLine.textContent = msg;
  statusLine.classList.toggle("error", isError);
}

async function fetchJson(path) {
  const res = await fetch(path);
  if (!res.ok) throw new Error(`${path} -> HTTP ${res.status}`);
  return res.json();
}

// ---------- 初期化 ----------
async function init() {
  state.lang = detectLang();
  t = STRINGS[state.lang];
  applyStaticStrings();
  document.querySelectorAll(".lang-switch button").forEach((b) => {
    b.addEventListener("click", () => setLang(b.dataset.lang));
  });

  panelsSection.style.display = "none";
  try {
    state.manifest = await fetchJson(`${DATA_BASE_URL}/manifest.json`);
  } catch (err) {
    setStatus(t.statusManifestError(err.message), true);
    return;
  }

  const cancerSelect = el("cancer-select");
  fillCancerSelect();
  const types = Object.keys(state.manifest.cancer_types || {});

  cancerSelect.addEventListener("change", onCancerTypeChange);
  el("gene-input").addEventListener("change", onGeneCommit);
  el("gene-input").addEventListener("keydown", (e) => {
    if (e.key === "Enter") onGeneCommit();
  });
  const logToggle = el("log-scale-toggle");
  if (logToggle) {
    logToggle.checked = state.logScale;
    logToggle.addEventListener("change", () => {
      state.logScale = logToggle.checked;
      if (state.currentGeneJson) renderAllPanels();
    });
  }
  const subtypeSelect = el("subtype-scheme-select");
  if (subtypeSelect) {
    subtypeSelect.addEventListener("change", () => {
      state.subtypeScheme = subtypeSelect.value;
      if (state.currentGeneJson) renderAllPanels();
    });
  }
  const extraToggle = el("extra-tumor-toggle");
  if (extraToggle) {
    extraToggle.checked = state.includeExtraTumor;
    extraToggle.addEventListener("change", () => {
      state.includeExtraTumor = extraToggle.checked;
      if (state.currentGeneJson) renderAllPanels();
    });
  }
  window.addEventListener("resize", () => {
    GENERATIONS.forEach((g) => {
      const node = el(g.plotId);
      if (node && node.dataset.rendered) Plotly.Plots.resize(node);
    });
  });

  setStatus(t.statusReady(types.length));
}

async function onCancerTypeChange(e) {
  const cancerType = e.target.value;
  state.cancerType = cancerType || null;
  const geneInput = el("gene-input");

  if (!cancerType) {
    geneInput.disabled = true;
    geneInput.value = "";
    return;
  }

  geneInput.disabled = true;
  geneInput.value = "";
  setStatus(t.statusGeneListLoading);
  hidePanels();

  try {
    state.geneIndex = await fetchJson(`${DATA_BASE_URL}/${cancerType}/_index.json`);
  } catch (err) {
    setStatus(t.statusGeneListError(err.message), true);
    return;
  }

  state.symbolToEnsembl = new Map();
  const datalist = el("gene-list");
  datalist.innerHTML = "";
  const frag = document.createDocumentFragment();
  for (const g of state.geneIndex) {
    if (!g.symbol) continue;
    const key = g.symbol.toUpperCase();
    if (!state.symbolToEnsembl.has(key)) state.symbolToEnsembl.set(key, g.ensembl);
    const opt = document.createElement("option");
    opt.value = g.symbol;
    frag.appendChild(opt);
  }
  datalist.appendChild(frag);

  geneInput.disabled = false;
  setStatus(t.statusGeneList(cancerLabel(cancerType, state.lang), state.geneIndex.length));

  // このがん種で使えるTumor分類方式をプルダウンに反映
  const subtypeSelect = el("subtype-scheme-select");
  if (subtypeSelect) {
    const schemes = (state.manifest.cancer_types[cancerType] || {}).subtype_schemes || [];
    subtypeSelect.innerHTML =
      `<option value="__none__">${t.subtypeNone}</option>` +
      schemes.map((s) => `<option value="${s}">${s}</option>`).join("");
    state.subtypeScheme = "__none__";
    subtypeSelect.value = "__none__";
    subtypeSelect.disabled = schemes.length === 0;
  }

  updateExtraTumorControl(cancerType);
}

// このがん種に転移巣・再発巣があるかをmanifestから見て、チェックボックスの
// 有効/無効と件数ラベルを更新する。がん種を変えるたびにOFFへ戻す。
function updateExtraTumorControl(cancerType) {
  const toggle = el("extra-tumor-toggle");
  const label = el("extra-tumor-label");
  if (!toggle || !label) return;

  const info = ((state.manifest.cancer_types || {})[cancerType] || {}).extra_tumor || {};
  const n = info.n_samples || 0;
  const byCode = info.by_code || {};

  state.includeExtraTumor = false;
  toggle.checked = false;
  toggle.disabled = n === 0;
  label.classList.toggle("is-disabled", n === 0);

  const textNode = label.querySelector("[data-i18n]");
  if (n === 0) {
    if (textNode) textNode.textContent = t.extraTumor;
    label.title = t.extraTumorNone;
  } else {
    const detail = Object.entries(byCode).map(([code, cnt]) => `${code}: ${cnt}`).join(", ");
    if (textNode) textNode.textContent = t.extraTumorWithCount(n);
    label.title = t.extraTumorHint(detail);
  }
}

async function onGeneCommit() {
  const raw = el("gene-input").value.trim();
  if (!raw || !state.cancerType) return;

  const ensembl = state.symbolToEnsembl.get(raw.toUpperCase());
  if (!ensembl) {
    setStatus(t.statusGeneNotFound(raw), true);
    hidePanels();
    return;
  }

  setStatus(t.statusGeneLoading(raw));
  try {
    state.currentGeneJson = await fetchJson(`${DATA_BASE_URL}/${state.cancerType}/${ensembl}.json`);
  } catch (err) {
    setStatus(t.statusGeneError(err.message), true);
    return;
  }

  renderAllPanels();
  setStatus(t.statusShowing(cancerLabel(state.cancerType, state.lang),
                            state.currentGeneJson.gene_symbol, ensembl));
}

function hidePanels() {
  panelsSection.style.display = "none";
  emptyState.style.display = "block";
}

function renderAllPanels() {
  emptyState.style.display = "none";
  panelsSection.style.display = "grid";

  const gj = state.currentGeneJson;
  // 3世代を通した並び順をここで1回だけ決め、以降のパネル描画で共有する
  state.boxOrder = computeBoxOrder(gj);
  for (const gen of GENERATIONS) {
    const genData = gj[gen.key]; // { TPM: {tumor:[],normal:[]}, FPKM_UQ: {...}, ... } または undefined
    const row = el(gen.rowId);
    const footnote = el(gen.footId);
    row.innerHTML = "";

    if (!genData || Object.keys(genData).length === 0) {
      Plotly.purge(el(gen.plotId));
      el(gen.plotId).dataset.rendered = "";
      // 「このがん種でそもそも取得していない」のか「この遺伝子が対応付かなかった」
      // のかは原因も対処も違うので、manifestを見て区別して伝える。
      const versions =
        (((state.manifest.cancer_types || {})[state.cancerType] || {}).versions || {})[gen.key] || [];
      renderEmptyPlot(
        gen.plotId,
        versions.length === 0 ? t.plotNoData : t.plotNoGene
      );
      footnote.textContent = "";
      el(gen.dlId).innerHTML = "";
      continue;
    }

    const valueTypes = Object.keys(genData);
    if (!state.activeValueType[gen.key] || !valueTypes.includes(state.activeValueType[gen.key])) {
      state.activeValueType[gen.key] = valueTypes[0];
    }

    if (valueTypes.length > 1) {
      for (const vt of valueTypes) {
        const btn = document.createElement("button");
        btn.textContent = vt;
        btn.className = state.activeValueType[gen.key] === vt ? "active" : "";
        btn.addEventListener("click", () => {
          state.activeValueType[gen.key] = vt;
          renderAllPanels();
        });
        row.appendChild(btn);
      }
    } else {
      // 1種類だけでも何の値かを必ず出す(がん種によって FPKM / FPKM-UQ が違うため)
      const chip = document.createElement("span");
      chip.className = "value-type-static";
      chip.textContent = valueTypes[0];
      row.appendChild(chip);
    }

    const vt = state.activeValueType[gen.key];
    const group = genData[vt];
    drawBoxPlot(gen.plotId, group, vt);
    renderDownloadRow(gen, vt);

    const nNormal = group.normal.length;
    const boxes = tumorBoxes(group);
    const nExtra = (group.tumor_extra || []).length;

    const head =
      boxes.length === 1 && boxes[0].name === "__TUMOR__"
        ? t.footTumor(boxes[0].values.length)
        : `${t.tumor}: ${boxes.map((b) => `${b.name} n=${b.values.length}`).join(" / ")}`;

    let extraNote = "";
    if (nExtra > 0) {
      extraNote = state.includeExtraTumor
        ? t.footExtraIncluded(nExtra)
        : t.footExtraExcluded(nExtra);
    }

    // 検体数が少ない群があれば注意を促す。膵臓がんのように正常組織が
    // 数検体しかないがん種で、箱ひげ図が過剰に信頼されるのを防ぐ。
    const small = [...boxes.map((b) => ({ name: b.name, n: b.values.length })),
                   { name: t.normalTissue, n: nNormal }]
      .filter((g) => g.n > 0 && g.n < SMALL_N_THRESHOLD);
    const smallNote = small.length
      ? t.smallNWarning(small.map((g) => `${g.name} n=${g.n}`).join(", "))
      : "";

    footnote.textContent =
      `${head} / ${t.footNormal(nNormal)}${extraNote}${smallNote}`;
    footnote.classList.toggle("has-warning", small.length > 0);
  }
}

// ---------- グラフの保存 ----------
// PNG(発表資料)、SVG(論文図版)、CSV(元の数値)の3種類を出せるようにする。
// Plotlyの標準ツールバーは非表示にしているため、自前でボタンを置く。
function renderDownloadRow(gen, valueLabel) {
  const row = el(gen.dlId);
  row.innerHTML = "";

  const label = document.createElement("span");
  label.className = "download-label";
  label.textContent = t.download;
  row.appendChild(label);

  const base = [
    state.cancerType,
    state.currentGeneJson ? state.currentGeneJson.gene_symbol : "gene",
    gen.key,
    valueLabel,
    state.subtypeScheme === "__none__" ? null : state.subtypeScheme,
  ].filter(Boolean).join("_").replace(/[^A-Za-z0-9_.-]+/g, "-");

  const buttons = [
    { text: t.downloadPng, title: t.downloadPngTitle, run: () => saveImage(gen.plotId, base, "png") },
    { text: t.downloadSvg, title: t.downloadSvgTitle, run: () => saveImage(gen.plotId, base, "svg") },
    { text: t.downloadCsv, title: t.downloadCsvTitle, run: () => saveCsv(gen.plotId, base) },
  ];
  for (const b of buttons) {
    const btn = document.createElement("button");
    btn.textContent = b.text;
    btn.title = b.title;
    btn.addEventListener("click", b.run);
    row.appendChild(btn);
  }
}

function saveImage(plotId, filename, format) {
  const node = el(plotId);
  if (!node || !node.dataset.rendered) return;
  // PNGは印刷や拡大に耐えるよう3倍解像度で書き出す。SVGは解像度非依存。
  Plotly.downloadImage(node, {
    format,
    filename,
    width: 900,
    height: 650,
    scale: format === "png" ? 3 : 1,
  });
}

function saveCsv(plotId, filename) {
  const node = el(plotId);
  if (!node || !node._exportRows) return;
  // 対数変換前の生の値を出す(log表示にしていても、保存されるのは元の値)
  const lines = [`${t.csvHeaderGroup},${t.csvHeaderValue}`];
  for (const [name, values] of node._exportRows) {
    const safe = /[",]/.test(name) ? `"${name.replace(/"/g, '""')}"` : name;
    for (const v of values) lines.push(`${safe},${v}`);
  }
  // Excelで開いたときに日本語が化けないよう BOM を付ける
  const blob = new Blob(["\uFEFF" + lines.join("\n")], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `${filename}.csv`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// 3世代に共通の「箱の並び順」を決める。
//
// 世代ごとに件数順で並べると、検体数の違い(膵臓がんは旧世代43 / 新世代150)で
// パネルごとに順序が変わってしまい、横に並べて比較できない。
// そこで全世代分を合算して順序を1回だけ決め、3パネルすべてに適用する。
//
// 並べ方は2通り:
//   (a) ラベルが共通の接頭辞+番号なら番号順  … CMS1-4, C1-C6, Cluster 1-3 など
//   (b) それ以外は合計件数の多い順          … PAM50, Moffitt, MSI など
// "Other" は常に右端に置く。
const OTHER_LABELS = new Set(["Other", "その他"]);

function computeBoxOrder(gj) {
  const totals = new Map();
  const prev = state.boxOrder;
  state.boxOrder = null;   // 集計中は件数順の既定動作に任せる
  for (const gen of GENERATIONS) {
    const genData = gj[gen.key];
    if (!genData) continue;
    for (const vt of Object.keys(genData)) {
      for (const box of tumorBoxes(genData[vt])) {
        totals.set(box.name, (totals.get(box.name) || 0) + box.values.length);
      }
    }
  }
  state.boxOrder = prev;

  const names = [...totals.keys()];
  const main = names.filter((n) => !OTHER_LABELS.has(n));

  // 共通の接頭辞+番号か判定する ("CMS1" -> 接頭辞 "CMS" / 番号 1)
  const parsed = main.map((n) => {
    const m = String(n).match(/^(\D*?)(\d+)/);
    return m ? { name: n, prefix: m[1].trim(), num: parseInt(m[2], 10) } : null;
  });
  const numbered = parsed.length > 1 && parsed.every(Boolean) &&
                   new Set(parsed.map((p) => p.prefix)).size === 1;

  const ordered = numbered
    ? parsed.slice().sort((a, b) => a.num - b.num).map((p) => p.name)
    : main.slice().sort((a, b) => (totals.get(b) - totals.get(a)) || a.localeCompare(b));

  return [...ordered, ...names.filter((n) => OTHER_LABELS.has(n))];
}

// 描画すべきTumor側の箱を [{name, values}] で返す。
// 分類方式の選択と「転移巣・再発巣を含む」チェックの状態をここで一元的に解決し、
// boxplot本体と脚注の両方が同じ内訳を見るようにしている。
function tumorBoxes(group) {
  const includeExtra = state.includeExtraTumor;
  const scheme = state.subtypeScheme;
  const bySubtype = group.tumor_by_scheme && group.tumor_by_scheme[scheme];
  const extraBySubtype = group.tumor_extra_by_scheme && group.tumor_extra_by_scheme[scheme];

  if (bySubtype || (includeExtra && extraBySubtype)) {
    const merged = {};
    for (const [name, vals] of Object.entries(bySubtype || {})) merged[name] = vals.slice();
    if (includeExtra) {
      for (const [name, vals] of Object.entries(extraBySubtype || {})) {
        merged[name] = (merged[name] || []).concat(vals);
      }
    }
    const order = state.boxOrder;
    const rank = (name) => {
      if (!order) return null;
      const i = order.indexOf(name);
      return i === -1 ? order.length : i;   // 想定外のラベルは末尾へ
    };
    return Object.entries(merged)
      .sort((a, b) => (order ? rank(a[0]) - rank(b[0]) : b[1].length - a[1].length))
      .map(([name, values]) => ({ name, values }));
  }

  let values = group.tumor;
  if (includeExtra && group.tumor_extra) values = values.concat(group.tumor_extra);
  // 分類なしの場合の腫瘍群。名前は描画直前に言語に応じて解決する。
  return [{ name: "__TUMOR__", values }];
}

function drawBoxPlot(plotId, group, valueLabel) {
  const useLog = state.logScale;
  // log2(x+1) 変換: 裾が長い分布(FPKM-UQ, raw countなど)で片方の群の箱が
  // もう片方の外れ値に押しつぶされて見えなくなるのを防ぐ。hoverには元の値を表示する。
  const transform = (arr) => (useLog ? arr.map((v) => Math.log2(v + 1)) : arr);
  const hoverText = (arr) => arr.map((v) => `${valueLabel}: ${v}`);

  const rawBoxes = tumorBoxes(group);
  // "__TUMOR__" は「分類なし」のときの腫瘍群を表すプレースホルダ。
  // 言語に応じてここで実際のラベルに置き換える。
  const boxes = rawBoxes.map((b) => ({
    name: b.name === "__TUMOR__" ? t.tumor : b.name,
    values: b.values,
  }));
  const traces = [];
  const useSubtypes = !(rawBoxes.length === 1 && rawBoxes[0].name === "__TUMOR__");

  boxes.forEach((box, i) => {
    const color = useSubtypes ? SUBTYPE_PALETTE[i % SUBTYPE_PALETTE.length] : COLORS.tumor;
    traces.push({
      y: transform(box.values),
      text: hoverText(box.values),
      hovertemplate: "%{text}<extra>" + box.name + "</extra>",
      type: "box",
      name: box.name,
      boxpoints: "outliers",
      marker: { color },
      fillcolor: useSubtypes ? color + "59" /* ~35% alpha */ : COLORS.tumorSoft,
      line: { color },
    });
  });

  const normalName = normalTraceName(boxes);

  traces.push({
    y: transform(group.normal),
    text: hoverText(group.normal),
    hovertemplate: "%{text}<extra>" + normalName + "</extra>",
    type: "box",
    name: normalName,
    boxpoints: "outliers",
    marker: { color: COLORS.normal },
    fillcolor: COLORS.normalSoft,
    line: { color: COLORS.normal },
  });

  const yTitle = useLog ? `log2(${valueLabel} + 1)` : valueLabel;

  // 腫瘍側の箱を左から順に並べ、非がん部を必ず右端に固定する。
  // categoryorder を指定しないと、Plotlyがトレースの出現順や値で並べ替える。
  const categoryArray = [...boxes.map((b) => b.name), normalName];

  const layout = {
    margin: { l: 46, r: 12, t: 10, b: 46 },
    yaxis: { title: yTitle, gridcolor: COLORS.rule, zeroline: false },
    xaxis: {
      gridcolor: COLORS.rule,
      tickangle: useSubtypes ? -20 : 0,
      categoryorder: "array",
      categoryarray: categoryArray,
    },
    // 腫瘍側と非がん部の境目に薄い縦線を引き、群が別物であることを示す
    shapes: [
      {
        type: "line",
        xref: "x",
        yref: "paper",
        x0: boxes.length - 0.5,
        x1: boxes.length - 0.5,
        y0: 0,
        y1: 1,
        line: { color: COLORS.rule, width: 1, dash: "dot" },
      },
    ],
    showlegend: false,
    // 図の英数字はArialで統一する(投稿規定で指定されることが多い)。
    // 日本語ラベルはArialに無いので、和文フォントへ自動的に落ちる。
    font: { family: PLOT_FONT, size: 11, color: COLORS.ink },
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
  };

  Plotly.newPlot(plotId, traces, layout, { displayModeBar: false, responsive: true });
  const node = el(plotId);
  node.dataset.rendered = "1";
  // CSV出力用に、変換前の生の値を群ごとに保持しておく
  node._exportRows = [
    ...boxes.map((b) => [b.name, b.values]),
    [normalName, group.normal],
  ];
  node._exportLabel = valueLabel;
}

function renderEmptyPlot(plotId, message) {
  const layout = {
    margin: { l: 10, r: 10, t: 10, b: 10 },
    xaxis: { visible: false },
    yaxis: { visible: false },
    annotations: [
      {
        text: message,
        showarrow: false,
        font: { family: PLOT_FONT, size: 12, color: COLORS.ink },
      },
    ],
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
  };
  Plotly.newPlot(plotId, [], layout, { displayModeBar: false, responsive: true });
}

hidePanels();
init();
