#!/usr/bin/env bash
# Usage: bash scripts/run_pipeline.sh
# Requires: conda env "bio-p2" (see README), plus R with DESeq2, apeglm, ggplot2, ggrepel, pheatmap
set -euo pipefail
ENV_NAME="bio-p2"

# Always run from the project root so relative paths (data/, results/) work
cd "$(dirname "${BASH_SOURCE[0]}")/.."

# Make `conda activate` available inside this script
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  . "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
else
  echo "Warning: conda not found; using python from PATH." >&2
fi

echo "Step 1: merge GSM count files into a matrix"
python scripts/merge_gsm_to_matrix.py

echo "Step 2: DESeq2 analysis"
Rscript scripts/deseq2.R

echo "Copying final results to results/"
mkdir -p results

# Copy a file if it exists
cp_if_exists() {
  if [ -f "$1" ]; then
    cp -f "$1" "$2"
    echo "Copied: $1 -> $2"
  else
    echo "Missing (skipped): $1"
  fi
}

D=data/processed/deseq2

# Shrunken DEG tables
cp_if_exists "$D/Basal_vs_LP_shrunken.csv" results/Basal_vs_LP_DEGs_full.csv
cp_if_exists "$D/Basal_vs_ML_shrunken.csv" results/Basal_vs_ML_DEGs_full.csv

# Concise tables, normalized counts, session info, and plots
for f in Basal_vs_LP_lfc_padj_table.csv Basal_vs_ML_lfc_padj_table.csv \
         Combined_lfc_padj_Basal_vs_LP_vs_ML.csv normalized_counts_vst.csv \
         sessionInfo.txt PCA.png MA_Basal_vs_LP.png MA_Basal_vs_ML.png \
         Heatmap_top50_Basal_vs_LP.png Heatmap_top50_Basal_vs_ML.png; do
  cp_if_exists "$D/$f" "results/$f"
done

echo "Pipeline finished successfully."
echo "Results available in: results/ and data/processed/deseq2/"
