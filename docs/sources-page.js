// ============================================================
// データの出典ページ
// ============================================================
// 本体(app.js)とは独立して動く。言語切替の仕組みだけ共有する。

const state = { lang: "ja" };
let t = STRINGS.ja;
const LANG_KEY = "tcga-explorer-lang";
const el = (id) => document.getElementById(id);

function detectLang() {
  try {
    const saved = localStorage.getItem(LANG_KEY);
    if (saved && STRINGS[saved]) return saved;
  } catch (e) { /* プライベートモード等では無視 */ }
  return (navigator.language || "en").toLowerCase().startsWith("ja") ? "ja" : "en";
}

function setLang(lang) {
  if (!STRINGS[lang]) return;
  state.lang = lang;
  t = STRINGS[lang];
  try { localStorage.setItem(LANG_KEY, lang); } catch (e) { /* 無視 */ }
  render();
}

function localized(entry, key) {
  // { ja: "...", en: "..." } を持つ項目から現在の言語のものを返す
  return entry[state.lang] || entry.en || entry[key] || "";
}

// --- がん種ごとの表 ----------------------------------------------------
function renderCancerTables() {
  const host = el("cancer-tables");
  host.innerHTML = "";
  const genLabel = { old: t.genOld, mid: t.genMid, new: t.genNew };

  const codes = Object.keys(DATA_SOURCES).sort((a, b) =>
    cancerLabel(a, state.lang).localeCompare(cancerLabel(b, state.lang), state.lang));

  for (const code of codes) {
    const block = document.createElement("section");
    block.className = "cancer-block";

    const h = document.createElement("h3");
    h.textContent = cancerLabel(code, state.lang);
    block.appendChild(h);

    const table = document.createElement("table");
    table.className = "source-table";
    const thead = document.createElement("thead");
    thead.innerHTML =
      `<tr><th>${t.colGeneration}</th><th>${t.colPipeline}</th>` +
      `<th>${t.colValue}</th><th>${t.colSource}</th><th>${t.colAcquired}</th></tr>`;
    table.appendChild(thead);

    const tbody = document.createElement("tbody");
    for (const row of DATA_SOURCES[code]) {
      const tr = document.createElement("tr");
      const def = SOURCE_DEFS[row.source] || {};
      const srcCell = def.url
        ? `<a href="${def.url}" rel="noopener">${def.name}</a>`
        : `${def.name || row.source}<span class="retired"> (${t.retired})</span>`;
      tr.innerHTML =
        `<td>${genLabel[row.generation] || row.generation}</td>` +
        `<td>${row.pipeline}</td>` +
        `<td><code>${row.value}</code></td>` +
        `<td>${srcCell}</td>` +
        `<td>${row.acquired}</td>`;
      tbody.appendChild(tr);

      const note = row[state.lang] || row.en;
      if (note) {
        const nr = document.createElement("tr");
        nr.className = "note-row";
        nr.innerHTML = `<td colspan="5">${note}</td>`;
        tbody.appendChild(nr);
      }
    }
    table.appendChild(tbody);
    block.appendChild(table);
    host.appendChild(block);
  }
}

// --- 取得元ごとの補足 --------------------------------------------------
function renderSourceNotes() {
  const dl = el("source-notes");
  dl.innerHTML = "";
  for (const key of Object.keys(SOURCE_DEFS)) {
    const def = SOURCE_DEFS[key];
    const dt = document.createElement("dt");
    dt.textContent = def.name + (def.retired ? ` (${t.retired})` : "");
    const dd = document.createElement("dd");
    dd.textContent = localized(def);
    dl.appendChild(dt);
    dl.appendChild(dd);
  }
}

// --- 分類の出典 --------------------------------------------------------
function renderSubtypeSources() {
  const ul = el("subtype-sources");
  ul.innerHTML = "";
  for (const s of SUBTYPE_SOURCES) {
    const li = document.createElement("li");
    const label = document.createElement("strong");
    label.textContent = s.label;
    li.appendChild(label);
    li.appendChild(document.createElement("br"));
    if (s.url) {
      const a = document.createElement("a");
      a.href = s.url;
      a.rel = "noopener";
      a.textContent = s.citation;
      li.appendChild(a);
    } else {
      li.appendChild(document.createTextNode(s.citation));
    }
    const used = document.createElement("span");
    used.className = "used-in";
    used.textContent = " — " + t.usedIn + ": " +
      s.used_in.map((c) => cancerLabel(c, state.lang)).join(", ");
    li.appendChild(used);
    ul.appendChild(li);
  }
}

// --- 値の種類 ----------------------------------------------------------
function renderValueTypeGlossary() {
  const ul = el("value-type-glossary");
  if (!ul || typeof VALUE_TYPE_INFO === "undefined") return;
  ul.innerHTML = "";
  for (const vt of Object.keys(VALUE_TYPE_INFO)) {
    const li = document.createElement("li");
    const code = document.createElement("code");
    code.textContent = vt;
    li.appendChild(code);
    li.appendChild(document.createTextNode(" — " + valueTypeInfo(vt, state.lang)));
    ul.appendChild(li);
  }
}

function render() {
  document.documentElement.lang = t.htmlLang;
  document.title = t.sourcesPageTitle + " — TCGA Expression Explorer";
  document.querySelectorAll("[data-i18n]").forEach((node) => {
    const v = t[node.dataset.i18n];
    if (typeof v === "string") node.textContent = v;
  });
  document.querySelectorAll(".lang-switch button").forEach((b) => {
    b.classList.toggle("active", b.dataset.lang === state.lang);
  });
  renderCancerTables();
  renderSourceNotes();
  renderSubtypeSources();
  renderValueTypeGlossary();
}

state.lang = detectLang();
t = STRINGS[state.lang];
document.querySelectorAll(".lang-switch button").forEach((b) => {
  b.addEventListener("click", () => setLang(b.dataset.lang));
});
render();
