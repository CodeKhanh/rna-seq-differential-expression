# Scripts Explained

### inspect_one_gsm.py

Loads a single GSM file to confirm:
* 3-column structure
* gene ID format
* count column integrity

### merge_gsm_to_matrix.py

* Finds GSM count files in data/dev/
* Extracts raw read counts
* Merges samples into a single matrix (genes × samples)
* Applies simple low-count filtering
* Writes:
  - counts_merged_unfiltered.csv
  - counts_cleaned.csv

### deseq2.R

Performs differential expression analysis:
* Reads filtered count matrix + metadata
* Sets Basal as reference group
* Runs DESeq2 (negative binomial model)
* Computes:
  - Basal vs LP
  - Basal vs ML
* Applies log2 fold-change shrinkage using apeglm
* Adjusts p-values using Benjamini–Hochberg
* Generates:
  - PCA plot
  - MA plots
  - Heatmaps (top 50 genes per contrast)
* Saves final tables and plots

Outputs are written to:
* data/processed/deseq2/
* results/

### run_pipeline.sh
Reproducible script that:
1. Activates conda environment
2. Runs merge script
3. Runs DESeq2
4. Copies final results into results/

## Key Outputs
Inside results/:
* Basal_vs_LP_DEGs_full.csv
* Basal_vs_ML_DEGs_full.csv
* PCA.png
* MA_Basal_vs_LP.png
* MA_Basal_vs_ML.png
* Heatmap_top50_Basal_vs_LP.png
* Heatmap_top50_Basal_vs_ML.png
