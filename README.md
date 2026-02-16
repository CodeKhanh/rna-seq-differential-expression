# RNA-seq differential expression pipeline

Overview
A compact, reproducible pipeline that downloads a small GEO count matrix, preprocesses sample metadata, runs DESeq2 for differential expression, and produces PCA and volcano plots.

Data
Example dataset: GEO accession GSEXXXXX (small human / mouse dataset with simple two-group design). Place raw counts in `data/raw/`.

Quick start
1. Create environment: `conda env create -f environment.yml` (or `pip install -r requirements.txt`)
2. Place counts in `data/raw/counts.csv` and metadata in `data/raw/metadata.csv`
3. Preprocess: `python scripts/preprocess.py`
4. Run DESeq2: `Rscript scripts/run_deseq2.R`
5. View results: `results/plots/` and `results/tables/`
environment.yml 
```bash 

conda env create -f environment.yml
conda activate rna-seq-de

```

## Reproducibility

Create environment:
conda env create -f environment.yml
conda activate rna-seq-de

Run full pipeline:
bash scripts/run_pipeline.sh

Outputs:
data/processed/deseq2/
