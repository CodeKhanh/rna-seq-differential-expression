#!/usr/bin/env python3
# scripts/merge_gsm_to_matrix.py
# Merge per-sample GSM files (3-column: gene_id, gene_len, count) into one counts matrix, 
# filter low-count genes, and write outputs to data/processed/.

from pathlib import Path
import glob, gzip, os, re
import pandas as pd
SCRIPT_DIR = Path(r"C:\Users\jilly\OneDrive\CodeKhanh\rna-seq-differential-expression\scripts\merge_gsm_to_matrix.py").resolve().parent

ROOT = SCRIPT_DIR.parent / "data"
DEV = ROOT / "dev" # local small files for development
PROCESSED = ROOT / "processed" # output directory for merged and filtered matrices 

print(f"Terminal is at: {Path.cwd()}")
print(f"Looking for data in: {ROOT.absolute()}")
print(f"Does data folder exist? {ROOT.exists()}")
PROCESSED.mkdir(parents=True, exist_ok=True)  

# 1) Find GSM files (adjust pattern if your filenames differ)
files = sorted(glob.glob(str(DEV / "GSM*.gz")) + glob.glob(str(DEV / "GSM*.txt")))
if not files:
    raise SystemExit("No GSM files found in data/dev. Put some GSM*.gz or GSM*.txt there.")

dfs = []
for f in files:
    fname = os.path.basename(f)
    # sample id: take leading GSMxxxxx part
    sample = re.split(r"[._\-]", fname)[0]
    print("Reading", fname, "as sample", sample)
    # open .gz or plain text
    opener = gzip.open if fname.endswith(".gz") else open
    with opener(f, "rt") as fh:
        # read three-column file, with header
        df = pd.read_csv(fh, sep="\t", header=0, names=["gene_id","gene_len","count"], index_col=0)
        # keep only counts and name the column after the sample id
        df = df[["count"]].rename(columns={"count": sample})
        dfs.append(df)

# 2) Outer join on gene_id and fill missing with 0
merged = pd.concat(dfs, axis=1, join="outer").fillna(0).astype(int)
print("Merged shape (genes, samples):", merged.shape)

# 3) Simple low-count filter
n_samples = merged.shape[1]
min_samples = min(3, max(1, n_samples//2))   # at least 3 or half of samples
min_count = 10
keep = (merged >= min_count).sum(axis=1) >= min_samples
filtered = merged.loc[keep]
print(f"Filtering: kept {filtered.shape[0]} genes out of {merged.shape[0]} (min_count={min_count}, min_samples={min_samples})")

# 4) Write outputs
merged.to_csv(PROCESSED / "counts_merged_unfiltered.csv")
filtered.to_csv(PROCESSED / "counts_cleaned.csv")
print("Wrote:", PROCESSED / "counts_merged_unfiltered.csv")
print("Wrote:", PROCESSED / "counts_cleaned.csv")