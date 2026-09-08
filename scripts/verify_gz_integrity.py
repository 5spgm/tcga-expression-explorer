import os, sys, gzip, json
from multiprocessing import Pool
D = sys.argv[1]
def chk(f):
    p = os.path.join(D, f)
    try:
        with gzip.open(p + ".gz", "rb") as fh: g = fh.read()
    except Exception as e:
        return (f, "gz読み込み失敗", f"{type(e).__name__} {e}"[:80])
    with open(p, "rb") as fh: raw = fh.read()
    if g != raw:
        return (f, "gzと.jsonが不一致", f"{len(raw)} vs {len(g)} B")
    try:
        json.loads(g)
    except Exception as e:
        return (f, "JSONとして不正", f"{type(e).__name__} {e}"[:80])
    return None
files = sorted(x for x in os.listdir(D) if x.endswith(".json") and os.path.exists(os.path.join(D, x + ".gz")))
with Pool(16) as pool:
    bad = [r for r in pool.imap_unordered(chk, files, chunksize=64) if r]
print(f"{D}\n  .gz と対で存在 {len(files):,} / 異常 {len(bad)}")
for b in bad[:10]: print("   ", b)
