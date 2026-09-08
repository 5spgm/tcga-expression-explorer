"""スプレッドシートとテキストの2つの発現行列を、遺伝子×検体で突き合わせる。

  python3 compare_matrix.py <ファイルA> <ファイルB>

行名・列名を揃えたうえで共通部分の最大相対差・最大絶対差を報告する。
Excel 往復で失われた桁があれば相対差として現れる。
"""
import sys, numpy as np, pandas as pd

def load(p):
    if p.endswith((".xlsx", ".xls")):
        df = pd.read_excel(p, index_col=0)
    else:
        df = pd.read_csv(p, sep="\t", index_col=0, low_memory=False)
    df.index = df.index.astype(str).str.strip()
    df.columns = [str(c).strip() for c in df.columns]
    return df

A, B = load(sys.argv[1]), load(sys.argv[2])
print(f"A: {A.shape[0]:,} 行 × {A.shape[1]:,} 列   {sys.argv[1].split('/')[-1]}")
print(f"B: {B.shape[0]:,} 行 × {B.shape[1]:,} 列   {sys.argv[2].split('/')[-1]}")
gi = A.index.intersection(B.index)
co = A.columns.intersection(B.columns)
print(f"共通: 行 {len(gi):,}  列 {len(co):,}")
print(f"A だけの列 {len(A.columns.difference(B.columns))} / B だけの列 {len(B.columns.difference(A.columns))}")
print(f"A だけの行 {len(A.index.difference(B.index))} / B だけの行 {len(B.index.difference(A.index))}")
if len(gi) == 0 or len(co) == 0:
    sys.exit("共通部分なし — 行名・列名の形式が違う")
a = A.loc[gi, co].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
b = B.loc[gi, co].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
both_nan = np.isnan(a) & np.isnan(b)
one_nan = (np.isnan(a) ^ np.isnan(b)).sum()
d = np.abs(a - b)
d[both_nan] = 0.0
denom = np.maximum(np.abs(a), np.abs(b))
rel = np.where(denom > 0, d / np.where(denom > 0, denom, 1), 0.0)
rel[both_nan] = 0.0
print(f"比較した値 {a.size:,}  片側だけ欠損 {one_nan:,}")
print(f"最大絶対差 {np.nanmax(d):.6e}   最大相対差 {np.nanmax(rel):.6e}")
print(f"完全一致した値 {int((d == 0).sum()):,} / {a.size:,}")
