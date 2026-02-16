#!/usr/bin/env bash 
# Usage: bash scripts/run_pipeline.sh
# Requires: conda, R, Python. Create env with: conda env create -f environment.yml
set -euo pipefail
ENV_NAME="bio-p2"

CONDABASE=""
if command -v conda >/dev/null 2>&1; then
  CONDABASE="$(conda info --base 2>/dev/null || echo '')"
fi

# common Windows paths for conda (adjust as needed)
if [ -z "$CONDABASE" ]; then
  if [ -d "$HOME/miniconda3" ]; then
    CONDABASE="$HOME/miniconda3"
  elif [ -d "$HOME/miniconda3" ]; then
    CONDABASE="$HOME/miniconda3"
  elif [ -d "$HOME/miniconda3" ]; then
    CONDABASE="$HOME/miniconda3"
  fi
fi

if [ -n "$CONDABASE" ] && [ -f "$CONDABASE/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  . "$CONDABASE/etc/profile.d/conda.sh"
else
  echo "Warning: conda not found automatically. If the script fails, run:"
  echo ". /c/Users/jilly/miniconda3/etc/profile.d/conda.sh && conda activate bio-p2"
fi

# activate the conda environment
conda activate bio-p2

echo "Running Phase 2: merge GSM files"
python scripts/merge_gsm_to_matrix.py

echo "Running Phase 3: DESeq2 analysis"
Rscript scripts/deseq2.R

echo "Copying final results to results/"
mkdir -p results

# Helper to copy if present 
cp_if_exists() {
  src="$1"
  dst="$2"
  if [ -f "$src" ]; then
    cp -f "$src" "$dst"
    echo "Copied: $src -> $dst"
  else
    echo "Missing (skipped): $src"
  fi
}

# DEG tables (shrunken)
cp_if_exists data/processed/deseq2/Basal_vs_LP_shrunken.csv \
             results/Basal_vs_LP_DEGs_full.csv

cp_if_exists data/processed/deseq2/Basal_vs_ML_shrunken.csv \
             results/Basal_vs_ML_DEGs_full.csv

# Extracted lfc/padj tables
cp_if_exists data/processed/deseq2/Basal_vs_LP_lfc_padj_table.csv \
             results/Basal_vs_LP_lfc_padj_table.csv

cp_if_exists data/processed/deseq2/Basal_vs_ML_lfc_padj_table.csv \
             results/Basal_vs_ML_lfc_padj_table.csv

cp_if_exists data/processed/deseq2/Combined_lfc_padj_Basal_vs_LP_vs_ML.csv \
             results/Combined_lfc_padj_Basal_vs_LP_vs_ML.csv

# Normalized counts and session info
cp_if_exists data/processed/deseq2/normalized_counts_vst.csv \
             results/normalized_counts_vst.csv

cp_if_exists data/processed/deseq2/sessionInfo.txt \
             results/sessionInfo.txt

# Plots: PCA, MA (LP/ML), Heatmaps (LP/ML)
cp_if_exists data/processed/deseq2/PCA.png \
             results/PCA.png

cp_if_exists data/processed/deseq2/MA_Basal_vs_LP.png \
             results/MA_Basal_vs_LP.png

cp_if_exists data/processed/deseq2/MA_Basal_vs_ML.png \
             results/MA_Basal_vs_ML.png

cp_if_exists data/processed/deseq2/Heatmap_top50_Basal_vs_LP.png \
             results/Heatmap_top50_Basal_vs_LP.png

cp_if_exists data/processed/deseq2/Heatmap_top50_Basal_vs_ML.png \
             results/Heatmap_top50_Basal_vs_ML.png

echo "Pipeline finished successfully."
echo "Results available in: results/ and data/processed/deseq2/"
