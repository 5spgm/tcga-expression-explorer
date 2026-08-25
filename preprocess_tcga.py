#!/usr/bin/env python3
"""
preprocess_tcga.py
===================
TCGAの新旧3世代RNA-seqデータ(旧: Firehose normalized_results, 中期: FPKM-UQ,
新: STAR-Counts TPM)を、がん種ごと・遺伝子ごとの小さなJSONファイルに変換する
前処理スクリプト。

出力されたJSONは、GitHub Pages上の静的サイトから fetch() して
boxplot(腫瘍 vs 正常、旧/中期/新を並べて表示)を描くために使う。

## 想定入力ファイル
- --new-tpm      : gene_id(Ensembl, バージョン付), gene_name, gene_type, サンプル列(TCGA barcode)
                    例: 20221216_COAD_TPM.xls (実体はタブ区切りのことが多い)
- --new-fpkm / --new-fpkmuq : 上と同フォーマットで value が違うだけ(将来追加用、任意)
- --mid-fpkmuq   : 1列目 = Ensembl ID(バージョンなし), サンプル列(TCGA barcode)
                    例: COAD_FPKM_UQ.xls
- --mid-fpkm     : 同上フォーマットで FPKM 版(将来追加用、任意)
- --old-normcount: 1列目 = "SYMBOL|EntrezID", サンプル列(TCGA barcode)
                    例: COAD-HiseqV2-20150129.txt (GDC Firehose Level_3)

## 使い方(例)
    python3 preprocess_tcga.py \\
        --cancer-type COAD \\
        --new-tpm /path/to/20221216_COAD_TPM.xls \\
        --mid-fpkmuq /path/to/COAD_FPKM_UQ.xls \\
        --old-normcount /path/to/COAD-HiseqV2-20150129.txt \\
        --out-dir ./data

## 遺伝子IDの統一方針
新データ(--new-*)の gene_id(Ensembl, バージョン付)からバージョンを除去した
ものを「マスターID」とする。
- 中期データはもとよりバージョンなしEnsembl IDなので、そのままjoin。
- 旧データは "SYMBOL|EntrezID" 形式なので、まず gene_name(symbol)でjoinを試み、
  一致しないもの("?"含む)は Entrez ID を mygene.info に問い合わせて
  Ensembl IDを取得し、再度join。それでも解決しないものは
  `<out-dir>/<cancer_type>/_unmapped_old_genes.json` に記録して除外する。
  (mygene.info への問い合わせはこのスクリプトを実行するマシンのネット環境が必要)

## サンプル分類
TCGAバーコードの sample type code(例 TCGA-AA-3870-01A の "01")で3群に判定:
  Tumor      : 01/03/05/09  原発巣。既定でTumor群として描画される。
  TumorExtra : 02/04/06/07/08/40  再発巣・転移巣など。JSONには "tumor_extra"
               として別枠で書き出し、**既定では描画に含めない**。サイト側の
               チェックボックス「転移巣・再発巣を含む」をONにしたときだけ
               Tumor群に合流する。
  Normal     : 10〜14  非腫瘍部。
想定外のコードが出た場合は警告を出したうえで、腫瘍系(<10)なら安全側に倒して
TumorExtra(既定で除外)に入れる。判定結果と出現コード一覧は実行時に画面出力
するので、毎回目視確認すること。
"""

import argparse
import json
import re
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

BARCODE_RE = re.compile(r"^TCGA-\w{2}-\w{4}-(\d{2})([A-Z]?)$")


# ----------------------------------------------------------------------
# ファイル読み込み(拡張子がxlsでも中身がタブ区切りテキストのことがある)
# ----------------------------------------------------------------------
def read_matrix_file(path: str) -> pd.DataFrame:
    path = Path(path)
    tried = []

    if path.suffix.lower() in (".xlsx", ".xls"):
        try:
            return _fix_implicit_index(pd.read_excel(path), path)
        except Exception as e:  # noqa: BLE001
            tried.append(f"read_excel: {e}")

    for sep in ("\t", ","):
        try:
            df = pd.read_csv(path, sep=sep)
            if df.shape[1] > 1:
                return _fix_implicit_index(df, path)
        except Exception as e:  # noqa: BLE001
            tried.append(f"read_csv(sep={sep!r}): {e}")

    raise ValueError(f"{path} を読み込めませんでした。試行: {tried}")


def _fix_implicit_index(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Rの write.table() のデフォルト書式(1列目=row.namesにヘッダー名を
    付けない)で書き出されたTSV/CSVでは、ヘッダー行の項目数がデータ行より
    1つ少なくなる。この場合pandasは「先頭列はインデックスだろう」と自動
    推測し、遺伝子ID列がDataFrameの列ではなくインデックスに吸い込まれて
    しまう(結果、先頭のサンプル列が誤って遺伝子ID列として扱われ、
    多数の遺伝子が消失/破損する)。ここでそれを検出し、インデックスを
    列に戻す。
    """
    if not isinstance(df.index, pd.RangeIndex):
        print(f"  [注記] {path.name}: ヘッダー行とデータ行の列数不一致を検出したため、"
              f"先頭列(インデックスに吸い込まれていた遺伝子ID)を列として復元します")
        df = df.reset_index()
        # reset_index() で作られる列名は 'index' になることが多いので、
        # 後続の「1列目=遺伝子ID列」という前提に合わせて明示的に確認できるようにしておく
        first_col = df.columns[0]
        print(f"    -> 復元した列名: {first_col!r} (先頭数件: {df[first_col].head(3).tolist()})")
    return df


# ----------------------------------------------------------------------
# バーコード処理
# ----------------------------------------------------------------------
def barcode_columns(df: pd.DataFrame) -> list[str]:
    """TCGAバーコード形式に一致する列名だけを返す。一致しない列は警告表示。"""
    cols = []
    dropped = []
    for c in df.columns:
        if BARCODE_RE.match(str(c)):
            cols.append(c)
        elif str(c).upper().startswith("TCGA"):
            dropped.append(c)
    if dropped:
        print(f"  [警告] TCGAらしき列だがバーコード形式に一致せず除外: {dropped}",
              file=sys.stderr)
    return cols


# TCGAのsample type codeを3群に分ける。
#   Tumor      : 原発巣。既定でTumor群として描画される。
#   TumorExtra : 再発巣・転移巣など。JSONには別枠で書き出し、サイト側の
#                チェックボックス「転移巣・再発巣を含む」をONにしたときだけ
#                Tumor群に合流する(既定は除外)。
#   Normal     : 非腫瘍部。
PRIMARY_TUMOR_CODES = {1, 3, 5, 9}
EXTRA_TUMOR_CODES = {2, 4, 6, 7, 8, 40}
NORMAL_CODES = {10, 11, 12, 13, 14}

SAMPLE_TYPE_LABELS = {
    1: "Primary Solid Tumor",
    2: "Recurrent Solid Tumor",
    3: "Primary Blood Derived Cancer (Peripheral Blood)",
    4: "Recurrent Blood Derived Cancer (Bone Marrow)",
    5: "Additional - New Primary",
    6: "Metastatic",
    7: "Additional Metastatic",
    8: "Human Tumor Original Cells",
    9: "Primary Blood Derived Cancer (Bone Marrow)",
    10: "Blood Derived Normal",
    11: "Solid Tissue Normal",
    12: "Buccal Cell Normal",
    13: "EBV Immortalized Normal",
    14: "Bone Marrow Normal",
    40: "Recurrent Blood Derived Cancer (Peripheral Blood)",
}

_UNKNOWN_CODES: set = set()


def group_for_code(code: str) -> str:
    """sample type code(2桁文字列)を Tumor / TumorExtra / Normal に振り分ける。"""
    if not str(code).isdigit():
        return "Unknown"
    c = int(code)
    if c in PRIMARY_TUMOR_CODES:
        return "Tumor"
    if c in EXTRA_TUMOR_CODES:
        return "TumorExtra"
    if c in NORMAL_CODES:
        return "Normal"
    # 未知のコード。腫瘍系(<10)なら安全側に倒して「既定では除外」の枠に入れる。
    _UNKNOWN_CODES.add(str(code))
    return "TumorExtra" if c < 10 else "Normal"


def sample_type_group(barcode: str) -> str:
    m = BARCODE_RE.match(barcode)
    if not m:
        return "Unknown"
    return group_for_code(m.group(1))


def summarize_sample_types(cols: list[str], label: str) -> None:
    counts = defaultdict(int)
    for c in cols:
        m = BARCODE_RE.match(c)
        code = m.group(1) if m else "??"
        counts[code] += 1
    print(f"  [{label}] sample type code 内訳:")
    for code in sorted(counts):
        grp = group_for_code(code)
        name = SAMPLE_TYPE_LABELS.get(int(code), "不明なコード") if code.isdigit() else "バーコード不一致"
        note = "  ← 既定では除外(サイト側のチェックボックスで合流可)" if grp == "TumorExtra" else ""
        print(f"      {code} {name:<48} -> {grp:<10} {counts[code]:>5} 件{note}")
    if _UNKNOWN_CODES:
        print(f"  [警告] 想定外のsample type code {sorted(_UNKNOWN_CODES)} が含まれています。"
              f"割り当てが妥当か確認してください。", file=sys.stderr)


def patient_id(barcode: str) -> str:
    return barcode[:12]


def vial_letter(barcode: str) -> str:
    """バーコード末尾のvial文字(TCGA-A7-A13D-01A の 'A')。無ければ空文字。"""
    m = BARCODE_RE.match(str(barcode))
    return m.group(2) if m else ""


# ----------------------------------------------------------------------
# 検体フィルタ(重複vialの除去・明示的な検体除外)
# ----------------------------------------------------------------------
# main() で引数から設定する。ファイルごとに独立して適用する
# (どのvialが含まれるかはファイル=世代によって違うため)。
FILTER_OPTS = {
    "dedup_vials": "first",    # "first" | "last" | "none"
    "exclude_vials": set(),    # 例 {"B"}: このvial文字を無条件に除外
    "exclude_samples": set(),  # フルバーコードまたは患者IDの集合
}


def filter_sample_columns(sample_cols: list[str], label: str) -> list[str]:
    """1) 明示的に除外指定された検体を落とし、
    2) 同一患者×同一sample typeで複数vialがある場合に1つへ絞る。

    注意: vial文字(A/B/C)はFFPE由来かどうかを表すものではなく、同一検体から
    採られた何本目のバイアルかを示すだけ。FFPEかどうかを根拠に落としたい場合は
    GDCのsample sheetの `is_ffpe` からバーコード一覧を作り、
    --exclude-samples で渡すこと。
    """
    excluded = []
    kept = []
    for c in sample_cols:
        name = str(c)
        if name in FILTER_OPTS["exclude_samples"] or patient_id(name) in FILTER_OPTS["exclude_samples"]:
            excluded.append(name)
            continue
        v = vial_letter(name)
        if v and v in FILTER_OPTS["exclude_vials"]:
            excluded.append(name)
            continue
        kept.append(c)

    if excluded:
        print(f"  [{label}] 除外指定により {len(excluded)} 検体を除外しました "
              f"(例: {', '.join(excluded[:5])}{' ...' if len(excluded) > 5 else ''})")

    if FILTER_OPTS["dedup_vials"] == "none":
        return kept

    by_key = defaultdict(list)
    for c in kept:
        m = BARCODE_RE.match(str(c))
        by_key[(patient_id(str(c)), m.group(1)) if m else (str(c), "")].append(c)

    dropped = []
    keep_set = set()
    for key, group in by_key.items():
        if len(group) == 1:
            keep_set.add(group[0])
            continue
        ordered = sorted(group, key=lambda c: vial_letter(c))
        chosen = ordered[0] if FILTER_OPTS["dedup_vials"] == "first" else ordered[-1]
        keep_set.add(chosen)
        dropped.append((key, [vial_letter(c) for c in ordered], vial_letter(chosen)))

    if dropped:
        print(f"  [{label}] 同一患者×同一sample typeの重複vialを {len(dropped)} 組"
              f"検出し、vial '{FILTER_OPTS['dedup_vials']}' を採用しました:")
        for (pid, code), vials, chosen in dropped[:10]:
            others = [v for v in vials if v != chosen]
            print(f"      {pid}-{code}: {'/'.join(vials)} -> {chosen} を採用"
                  f"({'/'.join(others)} を除外)")
        if len(dropped) > 10:
            print(f"      ... 他 {len(dropped) - 10} 組")
        return [c for c in kept if c in keep_set]

    return kept


# ----------------------------------------------------------------------
# Tumorサブタイプ表(任意、複数指定可)
# ----------------------------------------------------------------------
def apply_min_subtype_n(patient_to_subtype: dict, tumor_patients: set, min_n: int) -> dict:
    """このがん種のTumor検体だけに絞った上で、出現数が min_n 未満のsubtypeを
    'Other' にまとめて返す(patient_id -> subtype の辞書、このがん種のTumor
    患者のみを含む)。"""
    from collections import Counter

    restricted = {p: s for p, s in patient_to_subtype.items() if p in tumor_patients}
    counts = Counter(restricted.values())
    small = {s for s, n in counts.items() if n < min_n}
    if small:
        print(f"    出現数 {min_n} 未満のsubtypeを 'Other' にまとめます: {sorted(small)}")
        restricted = {p: ("Other" if s in small else s) for p, s in restricted.items()}
    print(f"    Tumor検体 {len(tumor_patients)} 件中 {len(restricted)} 件にsubtypeが割り当てられました"
          f" / 内訳: {Counter(restricted.values())}")
    return restricted


def parse_subtype_scheme_arg(raw: str) -> tuple:
    """'表示名=path/to/csv.csv' または '表示名=path/to/csv.csv:min_n' 形式をパースする。
    末尾の ':数字' はこのスキーム専用の min-subtype-n 閾値(省略時は None、
    呼び出し側でグローバルな --min-subtype-n にフォールバックする)。
    """
    if "=" not in raw:
        sys.exit(f"--tumor-subtype-table は 'NAME=path.csv' 形式で指定してください: {raw}")
    name, rest = raw.split("=", 1)
    name, rest = name.strip(), rest.strip()

    m = re.match(r"^(.*):(\d+)$", rest)
    if m:
        path, min_n = m.group(1), int(m.group(2))
    else:
        path, min_n = rest, None
    return name, path, min_n


class ParseSchemeAction(argparse.Action):
    """--tumor-subtype-table 'iCluster (Hoadley 2018)=icluster_subtypes.csv' や
    '...=path.csv:3' (スキーム専用のmin-subtype-n付き)のように、複数回指定
    できるようにする argparse アクション。
    """

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest, None) or []
        try:
            items.append(parse_subtype_scheme_arg(values))
        except SystemExit as e:
            parser.error(str(e))
        setattr(namespace, self.dest, items)


# ----------------------------------------------------------------------
# 各世代のローダー
# ----------------------------------------------------------------------
def load_new(path: str, value_label: str):
    """新データ(gene_id, gene_name, gene_type, サンプル列...)を読み込む。
    戻り値: (matrix辞書{ensembl_ids, sample_cols, values(2次元numpy配列)}, gene_master_df[ensembl, symbol, gene_type])

    以前はここで wide -> long (melt) 変換していたが、数万遺伝子×数百検体規模だと
    long-format は数千万行に膨れ上がり、文字列列(遺伝子ID・バーコード)が
    行ごとに複製されるためメモリを大量に消費し、segmentation faultの原因に
    なっていた。wide行列のままnumpy 2次元配列として保持し、遺伝子ごとの
    アクセスは「行を1回スライスする」だけにする(メモリはほぼ元のファイル分だけ)。
    """
    print(f"[new:{value_label}] 読み込み中: {path}")
    df = read_matrix_file(path)

    id_col = next((c for c in df.columns if str(c).lower() in ("gene_id", "ensembl_id")), df.columns[0])
    name_col = next((c for c in df.columns if str(c).lower() in ("gene_name", "symbol", "gene_symbol")), None)
    type_col = next((c for c in df.columns if str(c).lower() in ("gene_type",)), None)

    df["_ensembl"] = df[id_col].astype(str).str.replace(r"\.\d+(_PAR_Y)?$", lambda m: m.group(1) or "", regex=True)

    sample_cols = barcode_columns(df)
    sample_cols = filter_sample_columns(sample_cols, f"new:{value_label}")
    summarize_sample_types(sample_cols, f"new:{value_label}")

    master = pd.DataFrame({"ensembl": df["_ensembl"]})
    master["symbol"] = df[name_col] if name_col else None
    master["gene_type"] = df[type_col] if type_col else None
    master = master.drop_duplicates(subset="ensembl")

    n_before = len(df)
    df = df.drop_duplicates(subset="_ensembl", keep="first")
    if len(df) < n_before:
        print(f"  [注記] 重複する遺伝子ID {n_before - len(df)} 件は最初の行のみ採用")

    values = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    matrix = {"ensembl_ids": df["_ensembl"].to_numpy(), "sample_cols": sample_cols, "values": values}

    print(f"  遺伝子数: {master.shape[0]}, サンプル数: {len(sample_cols)}")
    return matrix, master


def load_mid(path: str, value_label: str):
    print(f"[mid:{value_label}] 読み込み中: {path}")
    df = read_matrix_file(path)
    id_col = df.columns[0]
    df["_ensembl"] = df[id_col].astype(str).str.replace(r"\.\d+(_PAR_Y)?$", lambda m: m.group(1) or "", regex=True)

    sample_cols = barcode_columns(df)
    sample_cols = filter_sample_columns(sample_cols, f"mid:{value_label}")
    summarize_sample_types(sample_cols, f"mid:{value_label}")

    n_before = len(df)
    df = df.drop_duplicates(subset="_ensembl", keep="first")
    if len(df) < n_before:
        print(f"  [注記] 重複する遺伝子ID {n_before - len(df)} 件は最初の行のみ採用")

    values = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    matrix = {"ensembl_ids": df["_ensembl"].to_numpy(), "sample_cols": sample_cols, "values": values}

    print(f"  遺伝子数: {len(df)}, サンプル数: {len(sample_cols)}")
    return matrix


def load_old(path: str):
    print(f"[old] 読み込み中: {path}")
    df = read_matrix_file(path)

    # 旧世代のファイルには2つの流儀がある。
    #   (a) Xena / Firehose 形式: 1列目が "SYMBOL|EntrezID"
    #   (b) cBioPortal 形式    : "Hugo_Symbol" と "Entrez_Gene_Id" の2列に分かれている
    lower = {str(c).strip().lower(): c for c in df.columns}
    hugo_col = lower.get("hugo_symbol")
    entrez_col = lower.get("entrez_gene_id")

    if hugo_col is not None and entrez_col is not None:
        print("  形式: cBioPortal (Hugo_Symbol + Entrez_Gene_Id の2列)")
        # pandasのバージョンによって astype(str) がNaNを "nan" にする場合と
        # NaNのまま残す場合がある。notna() で明示的に判定してから変換する。
        sym = df[hugo_col]
        df["_symbol"] = sym.where(sym.notna(), "?").astype(str).str.strip()
        # 空欄は Xena 形式に合わせて "?" に寄せる(後段のEntrez解決に乗せるため)
        df.loc[df["_symbol"].isin(["", "nan", "NaN", "None", "<NA>"]), "_symbol"] = "?"

        ent = df[entrez_col]
        df["_entrez"] = (ent.where(ent.notna(), "").astype(str).str.strip()
                         .str.replace(r"\.0$", "", regex=True))
        df.loc[df["_entrez"].isin(["", "nan", "NaN", "None", "<NA>", "0"]), "_entrez"] = None
    else:
        id_col = df.columns[0]
        print(f"  形式: Xena/Firehose (1列目 '{id_col}' が SYMBOL|EntrezID)")
        split = df[id_col].astype(str).str.split("|", n=1, expand=True)
        df["_symbol"] = split[0]
        df["_entrez"] = split[1] if split.shape[1] > 1 else None

    sample_cols = barcode_columns(df)
    sample_cols = filter_sample_columns(sample_cols, "old")
    summarize_sample_types(sample_cols, "old")

    n_before = len(df)
    df = df.drop_duplicates(subset=["_symbol", "_entrez"], keep="first")
    if len(df) < n_before:
        print(f"  [注記] 重複する遺伝子キー {n_before - len(df)} 件は最初の行のみ採用")

    values = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    matrix = {
        "symbol": df["_symbol"].to_numpy(),
        "entrez": df["_entrez"].to_numpy(),
        "sample_cols": sample_cols,
        "values": values,
    }

    n_unknown = (df["_symbol"] == "?").sum()
    print(f"  遺伝子数: {len(df)} (symbol不明 '?': {n_unknown}件), サンプル数: {len(sample_cols)}")
    return matrix


# ----------------------------------------------------------------------
def resolve_symbols_from_ensembl(master: pd.DataFrame, cache_path: Path | None = None) -> pd.DataFrame:
    """中期データをマスターにする場合、gene_symbol の列が無いので
    mygene.info で Ensembl -> symbol を引く。サイトの遺伝子名検索が
    symbol依存なので、これが無いと何も検索できなくなる。
    問い合わせに失敗した場合は Ensembl ID を symbol の代わりに使い、
    最低限検索できる状態にしておく。"""
    ids = [e for e in master["ensembl"].astype(str).tolist() if e]
    cache = {}
    if cache_path and cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except Exception:  # noqa: BLE001
            cache = {}

    todo = [e for e in ids if e not in cache]
    if todo:
        try:
            import mygene
        except ImportError:
            print("  [警告] mygene がインストールされていません。"
                  "遺伝子シンボルが引けないため、Ensembl IDで代用します。", file=sys.stderr)
            master["symbol"] = master["ensembl"]
            return master

        print(f"  mygene.info に {len(todo):,} 件のEnsembl IDを問い合わせ中…"
              f"(件数が多いので数分かかります)")
        mg = mygene.MyGeneInfo()
        try:
            results = mg.querymany(todo, scopes="ensembl.gene", fields="symbol",
                                   species="human", returnall=False)
            for r in results:
                cache[r.get("query")] = r.get("symbol")
            if cache_path:
                cache_path.write_text(json.dumps(cache, ensure_ascii=False, indent=2))
        except Exception as e:  # noqa: BLE001
            print(f"  [警告] mygene.info への問い合わせに失敗しました: {e}", file=sys.stderr)

    master["symbol"] = master["ensembl"].map(lambda e: cache.get(e) or e)
    n_resolved = sum(1 for e in ids if cache.get(e))
    print(f"  シンボル解決: {n_resolved:,} / {len(ids):,} 件"
          f"(未解決分はEnsembl IDのまま検索できます)")
    return master


# ----------------------------------------------------------------------
# Entrez -> Ensembl/Symbol マッピング (mygene.info, 要ネット接続)
# ----------------------------------------------------------------------
def resolve_entrez_to_ensembl(entrez_ids: list[str], cache_path: Path) -> dict:
    cache = {}
    if cache_path.exists():
        cache = json.loads(cache_path.read_text())

    todo = [e for e in entrez_ids if e and e not in cache]
    if todo:
        try:
            import mygene
        except ImportError:
            print("  [警告] mygene がインストールされていません。"
                  "`pip install mygene` を実行してください。Entrez解決をスキップします。",
                  file=sys.stderr)
            return cache

        print(f"  mygene.info に {len(todo)} 件のEntrez IDを問い合わせ中...(ネット接続が必要)")
        mg = mygene.MyGeneInfo()
        try:
            results = mg.querymany(todo, scopes="entrezgene", fields="ensembl.gene,symbol",
                                    species="human", returnall=False)
        except Exception as e:  # noqa: BLE001
            print(f"  [警告] mygene.info への問い合わせに失敗しました: {e}", file=sys.stderr)
            return cache

        for r in results:
            eid = r.get("query")
            ens = r.get("ensembl")
            if isinstance(ens, list):
                ens = ens[0] if ens else None
            ens_id = ens.get("gene") if isinstance(ens, dict) else None
            cache[eid] = {"ensembl": ens_id, "symbol": r.get("symbol")}

        cache_path.write_text(json.dumps(cache, ensure_ascii=False, indent=2))

    return cache


# ----------------------------------------------------------------------
# メイン処理
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cancer-type", required=True, help="例: COAD")
    ap.add_argument("--new-tpm", help="新データ TPM ファイル")
    ap.add_argument("--new-fpkm", help="新データ FPKM ファイル(任意)")
    ap.add_argument("--new-fpkmuq", help="新データ FPKM-UQ ファイル(任意)")
    ap.add_argument("--mid-fpkmuq", help="中期データ FPKM-UQ ファイル")
    ap.add_argument("--mid-fpkm", help="中期データ FPKM ファイル(任意)")
    ap.add_argument("--old-normcount", help="旧データ normalized_count ファイル")
    ap.add_argument("--dedup-vials", choices=["first", "last", "none"], default="first",
                     help="同一患者×同一sample typeに複数vial(01A/01Bなど)がある場合の扱い。"
                          "first(既定)=vial文字が若い方を採用 / last=遅い方を採用 / "
                          "none=重複を除去しない。同じ患者が2回カウントされるのを防ぐための処理。")
    ap.add_argument("--exclude-vials", default="",
                     help="無条件に除外するvial文字をカンマ区切りで指定(例: B)。"
                          "重複の有無に関わらず落とすので、そのvialしか無い患者は"
                          "まるごと消える。実行前に必ず件数を確認すること。")
    ap.add_argument("--exclude-samples",
                     help="除外する検体のバーコード(TCGA-XX-XXXX-01A)または患者ID"
                          "(TCGA-XX-XXXX)を1行1件で並べたテキスト/CSVファイル。"
                          "FFPE由来検体を落としたい場合は、GDCのsample sheetの is_ffpe から"
                          "このファイルを作って渡す(vial文字はFFPEかどうかを表さない)。")
    ap.add_argument("--out-dir", default="./data", help="出力先ディレクトリ")
    ap.add_argument("--entrez-cache", default="./entrez_cache.json",
                     help="Entrez->Ensembl マッピングのキャッシュファイル(全がん種で共有可)")
    ap.add_argument("--tumor-subtype-table", action=ParseSchemeAction, dest="tumor_subtype_tables", default=[],
                     metavar="NAME=path.csv[:min_n]",
                     help="Tumor検体を分類するための patient_id,subtype 2列CSVを "
                          "'表示名=path.csv' の形式で指定。複数回指定すると、サイト上で"
                          "プルダウンから分類方式を切り替えられる。末尾に ':数字' を付けると"
                          "そのスキーム専用の min-subtype-n を指定できる(省略時は"
                          "--min-subtype-n の値を使う)。例: "
                          "--tumor-subtype-table \"iCluster (Hoadley 2018)=icluster_subtypes.csv\" "
                          "--tumor-subtype-table \"MSI Status (Liu 2018)=subtypes/MSI_Status.csv:1\" "
                          "省略時はTumor/Normalの2群のみ(分類方式のプルダウンは「なし」だけになる)。")
    ap.add_argument("--min-subtype-n", type=int, default=5,
                     help="このがん種内での出現数がこの値未満のsubtypeは'Other'にまとめる既定の閾値"
                          "(既定: 5)。スキームごとに変えたい場合は --tumor-subtype-table 側の "
                          "':数字' で上書きできる。")
    args = ap.parse_args()

    # --- 検体フィルタ設定 ---
    FILTER_OPTS["dedup_vials"] = args.dedup_vials
    FILTER_OPTS["exclude_vials"] = {
        v.strip().upper() for v in args.exclude_vials.split(",") if v.strip()
    }
    if args.exclude_samples:
        raw = Path(args.exclude_samples).read_text().splitlines()
        ids = set()
        for line in raw:
            token = line.split(",")[0].strip()
            if not token or token.lower() in ("barcode", "sample_id", "patient_id", "sample"):
                continue
            ids.add(token)
        FILTER_OPTS["exclude_samples"] = ids
        print(f"[filter] --exclude-samples: {len(ids)} 件の除外対象を読み込みました "
              f"({args.exclude_samples})")
    if FILTER_OPTS["exclude_vials"]:
        print(f"[filter] --exclude-vials: vial {sorted(FILTER_OPTS['exclude_vials'])} を"
              f"無条件に除外します(そのvialしか無い患者も消えます)")
    print(f"[filter] --dedup-vials: {args.dedup_vials}")

    # 分類方式(scheme)ごとに patient_id -> subtype の辞書と、専用min_nを読み込む
    scheme_maps: dict[str, dict] = {}
    scheme_min_n: dict[str, int] = {}
    for scheme_name, path, min_n_override in args.tumor_subtype_tables:
        raw_df = pd.read_csv(path)
        cols_lower = {c.lower(): c for c in raw_df.columns}
        if "patient_id" not in cols_lower or "subtype" not in cols_lower:
            sys.exit(f"{path} には patient_id, subtype の2列が必要です。見つかった列: {list(raw_df.columns)}")
        raw_df = raw_df.rename(columns={cols_lower["patient_id"]: "patient_id", cols_lower["subtype"]: "subtype"})
        scheme_maps[scheme_name] = dict(zip(raw_df["patient_id"], raw_df["subtype"]))
        scheme_min_n[scheme_name] = min_n_override if min_n_override is not None else args.min_subtype_n
        min_n_note = f"(このスキーム専用: {min_n_override})" if min_n_override is not None else f"(共通既定値: {args.min_subtype_n})"
        print(f"[subtype] '{scheme_name}' <- {path} から {len(scheme_maps[scheme_name])} 患者分を読み込みました "
              f"/ min-subtype-n {min_n_note}")

    out_dir = Path(args.out_dir) / args.cancer_type
    out_dir.mkdir(parents=True, exist_ok=True)

    # マスター遺伝子リストは Ensembl ID を持つ世代からしか作れない
    # (旧データは SYMBOL|EntrezID 形式のため単独ではマスターになれない)。
    if not (args.new_tpm or args.new_fpkm or args.new_fpkmuq
            or args.mid_fpkmuq or args.mid_fpkm):
        sys.exit("新世代(--new-*)または中期(--mid-*)のデータが最低1つ必要です"
                 "(マスター遺伝子リストの元になります)。旧データ単独では実行できません。")

    # --- 新データ(複数value_type対応) ---
    new_frames = {}
    master = None
    new_specs = [("TPM", args.new_tpm), ("FPKM", args.new_fpkm), ("FPKM_UQ", args.new_fpkmuq)]
    for label, path in new_specs:
        if not path:
            continue
        long_df, m = load_new(path, label)
        new_frames[label] = long_df
        master = m if master is None else master  # 最初に読んだものをマスターにする

    # --- 中期データ ---
    mid_frames = {}
    mid_specs = [("FPKM_UQ", args.mid_fpkmuq), ("FPKM", args.mid_fpkm)]
    for label, path in mid_specs:
        if not path:
            continue
        mid_frames[label] = load_mid(path, label)

    # 新世代が無いがん種では、中期をマスター遺伝子リストとして使う。
    # 中期もEnsembl IDなので同じ役割を果たせる。ただしgene_symbol /
    # gene_typeの列は持たないため、symbolはmygeneで後から補う。
    if master is None and mid_frames:
        first_label, first_matrix = next(iter(mid_frames.items()))
        print(f"[master] 新世代データが無いため、中期({first_label})を"
              f"マスター遺伝子リストとして使います")
        master = pd.DataFrame({"ensembl": first_matrix["ensembl_ids"]})
        master["symbol"] = None
        master["gene_type"] = None
        master = master.drop_duplicates(subset="ensembl")
        master = resolve_symbols_from_ensembl(
            master, Path(args.entrez_cache).with_name("ensembl_symbol_cache.json"))

    # --- 旧データ ---
    old_df = None
    if args.old_normcount:
        old_df = load_old(args.old_normcount)

    # --- 旧データのsymbol -> master ensembl 対応 ---
    if old_df is not None:
        sym_to_ens = master.dropna(subset=["symbol"]).copy()
        sym_to_ens["_symlow"] = sym_to_ens["symbol"].str.lower()
        sym_to_ens = dict(zip(sym_to_ens["_symlow"], sym_to_ens["ensembl"]))

        old_symbol = old_df["symbol"]
        old_entrez = old_df["entrez"]
        old_gene_keys = pd.DataFrame({"symbol": old_symbol, "entrez": old_entrez}).drop_duplicates()
        unmatched_rows = old_gene_keys[
            ~old_gene_keys["symbol"].str.lower().isin(sym_to_ens.keys())
        ]

        print(f"[old->master] symbol直接一致: "
              f"{len(old_gene_keys) - len(unmatched_rows)} / {len(old_gene_keys)} 遺伝子")

        entrez_cache = {}
        if not unmatched_rows.empty:
            cache_path = Path(args.entrez_cache)
            entrez_cache = resolve_entrez_to_ensembl(
                unmatched_rows["entrez"].dropna().unique().tolist(), cache_path
            )

        def resolve_symbol(sym, entrez):
            s = str(sym).lower()
            if s in sym_to_ens:
                return sym_to_ens[s]
            info = entrez_cache.get(entrez)
            if info and info.get("ensembl"):
                return info["ensembl"]
            return None

        resolved_ensembl = np.array(
            [resolve_symbol(s, e) for s, e in zip(old_symbol, old_entrez)], dtype=object
        )
        is_mapped = resolved_ensembl != None  # noqa: E711

        n_unmapped_genes = pd.Series(old_symbol[~is_mapped]).nunique()
        if n_unmapped_genes:
            unmapped_syms = sorted(
                pd.DataFrame({"symbol": old_symbol[~is_mapped], "entrez": old_entrez[~is_mapped]})
                .drop_duplicates().to_dict("records"),
                key=lambda x: str(x["symbol"]),
            )
            (out_dir / "_unmapped_old_genes.json").write_text(
                json.dumps(unmapped_syms, ensure_ascii=False, indent=2)
            )
            print(f"  [警告] 旧データで未マッピングの遺伝子 {n_unmapped_genes} 件 -> "
                  f"{out_dir / '_unmapped_old_genes.json'} に記録")

        # マッピングできた行だけ残し、行列も同じ行だけに絞る
        old_df["ensembl_ids"] = resolved_ensembl[is_mapped]
        old_df["values"] = old_df["values"][is_mapped]
        del old_df["symbol"], old_df["entrez"]

    # --- サンプル分類(全データ横断) ---
    all_barcodes = set()
    for mat in list(new_frames.values()) + list(mid_frames.values()) + ([old_df] if old_df is not None else []):
        all_barcodes |= set(mat["sample_cols"])
    sample_group = {b: sample_type_group(b) for b in all_barcodes}

    # --- Tumorサブタイプ(任意、複数スキームに対応) ---
    # scheme_name -> {barcode: subtype}
    barcode_to_subtype_by_scheme: dict[str, dict] = {}
    if scheme_maps:
        tumor_patients = {
            patient_id(b) for b, g in sample_group.items() if g in ("Tumor", "TumorExtra")
        }
        for scheme_name, patient_map in scheme_maps.items():
            print(f"  [subtype] スキーム '{scheme_name}' を適用中...")
            patient_subtype = apply_min_subtype_n(patient_map, tumor_patients, scheme_min_n[scheme_name])
            barcode_to_subtype_by_scheme[scheme_name] = {
                b: patient_subtype[patient_id(b)]
                for b in all_barcodes
                if sample_group.get(b) in ("Tumor", "TumorExtra") and patient_id(b) in patient_subtype
            }

    # --- 遺伝子ごとにJSON出力 ---
    print(f"\n[出力] {out_dir} に遺伝子ごとのJSONを書き出し中...")

    def build_lookup_from_matrix(mat: dict, label: str):
        """wide行列(ensembl_ids, sample_cols, values)から、遺伝子を引くための
        インデックス(辞書・boolマスク)だけを事前に作る。実データ本体
        (values)はコピーしない(参照を持つだけ)ので、メモリはほぼ増えない。
        """
        sample_cols = mat["sample_cols"]
        col_group = np.array([sample_group.get(b) for b in sample_cols], dtype=object)
        tumor_mask = col_group == "Tumor"
        extra_mask = col_group == "TumorExtra"
        normal_mask = col_group == "Normal"
        ensembl_to_rowidx = {ens: i for i, ens in enumerate(mat["ensembl_ids"])}

        def scheme_masks_for(base_mask):
            """base_mask(原発巣 or 転移・再発)の列だけについて、subtypeラベル別の
            列マスクを作る。"""
            out: dict[str, dict] = {}
            for scheme_name, barcode_to_subtype in barcode_to_subtype_by_scheme.items():
                subtypes_for_cols = [
                    barcode_to_subtype.get(b) if in_base else None
                    for b, in_base in zip(sample_cols, base_mask)
                ]
                subtypes_arr = np.array(subtypes_for_cols, dtype=object)
                unique_labels = sorted({s for s in subtypes_for_cols if s is not None})
                if not unique_labels:
                    continue
                out[scheme_name] = {lbl: (subtypes_arr == lbl) for lbl in unique_labels}
            return out

        scheme_col_masks = scheme_masks_for(tumor_mask)
        extra_scheme_col_masks = scheme_masks_for(extra_mask)

        n_extra = int(extra_mask.sum())
        extra_note = f" / うち転移・再発 {n_extra} 件は既定で除外" if n_extra else ""
        print(f"    [{label}] 集計完了(遺伝子 {len(ensembl_to_rowidx):,} × 検体 {len(sample_cols)}){extra_note}",
              flush=True)
        return (ensembl_to_rowidx, mat["values"], tumor_mask, extra_mask, normal_mask,
                scheme_col_masks, extra_scheme_col_masks)

    def make_group_dict(ensembl, lookup) -> dict | None:
        (ensembl_to_rowidx, values, tumor_mask, extra_mask, normal_mask,
         scheme_col_masks, extra_scheme_col_masks) = lookup
        idx = ensembl_to_rowidx.get(ensembl)
        if idx is None:
            return None
        row = values[idx]

        def take(mask):
            vals = row[mask]
            return vals[~np.isnan(vals)].tolist()

        def by_scheme_dict(masks_by_scheme):
            out = {}
            for scheme_name, label_masks in masks_by_scheme.items():
                d = {}
                for lbl, mask in label_masks.items():
                    vals = take(mask)
                    if vals:
                        d[lbl] = vals
                if d:
                    out[scheme_name] = d
            return out

        tumor = take(tumor_mask)
        extra = take(extra_mask)
        normal = take(normal_mask)
        if not tumor and not extra and not normal:
            return None

        result = {"tumor": tumor, "normal": normal}
        # 転移巣・再発巣は別枠。サイト側でチェックONのときだけ tumor に合流させる。
        if extra:
            result["tumor_extra"] = extra

        by_scheme = by_scheme_dict(scheme_col_masks)
        if by_scheme:
            result["tumor_by_scheme"] = by_scheme
        extra_by_scheme = by_scheme_dict(extra_scheme_col_masks)
        if extra_by_scheme:
            result["tumor_extra_by_scheme"] = extra_by_scheme
        return result

    print("  各バージョンを事前集計中(この後の遺伝子ループが高速化されます)...", flush=True)
    new_lookup = {label: build_lookup_from_matrix(mat, f"new:{label}") for label, mat in new_frames.items()}
    mid_lookup = {label: build_lookup_from_matrix(mat, f"mid:{label}") for label, mat in mid_frames.items()}
    old_lookup = build_lookup_from_matrix(old_df, "old:normalized_count") if old_df is not None else None

    index_records = []
    n_written = 0
    for _, row in master.iterrows():
        ensembl = row["ensembl"]
        gene_json = {
            "gene_id": ensembl,
            "gene_symbol": row.get("symbol"),
            "gene_type": row.get("gene_type"),
        }

        new_out = {}
        for label, lookup in new_lookup.items():
            g = make_group_dict(ensembl, lookup)
            if g:
                new_out[label] = g
        if new_out:
            gene_json["new"] = new_out

        mid_out = {}
        for label, lookup in mid_lookup.items():
            g = make_group_dict(ensembl, lookup)
            if g:
                mid_out[label] = g
        if mid_out:
            gene_json["mid"] = mid_out

        if old_lookup is not None:
            g = make_group_dict(ensembl, old_lookup)
            if g:
                gene_json["old"] = {"normalized_count": g}

        # 新データにしか無い(=old/midどちらとも紐付かない)遺伝子でも出力はする。
        # ただしどのバージョンにも値が無ければスキップ。
        if not any(k in gene_json for k in ("new", "mid", "old")):
            continue

        (out_dir / f"{ensembl}.json").write_text(
            json.dumps(gene_json, ensure_ascii=False)
        )
        index_records.append({"ensembl": ensembl, "symbol": row.get("symbol")})
        n_written += 1

    index_records.sort(key=lambda r: (r["symbol"] or "").upper())
    (out_dir / "_index.json").write_text(json.dumps(index_records, ensure_ascii=False))

    # --- マニフェスト更新(がん種一覧) ---
    manifest_path = Path(args.out_dir) / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"cancer_types": {}}
    # 転移巣・再発巣がこのがん種に存在するか(サイト側チェックボックスの
    # 有効/無効と、ラベルに出す件数に使う)
    extra_counts: dict[str, int] = defaultdict(int)
    for b in all_barcodes:
        if sample_group.get(b) == "TumorExtra":
            m = BARCODE_RE.match(b)
            extra_counts[m.group(1) if m else "??"] += 1

    manifest["cancer_types"][args.cancer_type] = {
        "n_genes": n_written,
        "filters": {
            "dedup_vials": args.dedup_vials,
            "exclude_vials": sorted(FILTER_OPTS["exclude_vials"]),
            "n_exclude_samples": len(FILTER_OPTS["exclude_samples"]),
        },
        "extra_tumor": {
            "n_samples": sum(extra_counts.values()),
            "by_code": {c: extra_counts[c] for c in sorted(extra_counts)},
        },
        "versions": {
            "new": list(new_frames.keys()),
            "mid": list(mid_frames.keys()),
            "old": ["normalized_count"] if old_df is not None else [],
        },
        "subtype_schemes": list(barcode_to_subtype_by_scheme.keys()),
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2))

    print(f"\n完了: {n_written} 遺伝子分のJSONを書き出しました -> {out_dir}")
    print(f"インデックス: {out_dir / '_index.json'}")
    print(f"マニフェスト: {manifest_path}")


if __name__ == "__main__":
    main()
