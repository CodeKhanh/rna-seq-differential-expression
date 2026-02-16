# RNA-seq differential expression pipeline 
## Dataset: GEO GSE63310
## Author: Jillian Higbee   

**Project Overview**
This project aims to perform differential gene expression analysis on RNA-seq data from mouse mammary epithelial cell populations: Basal (control), Luminal Progenitor (LP), and Mature Luminal (ML). It was built as a beginner bioinformatics piepeline and demonstrates reproducable analysis on VSC using: 
* Python (data cleaning)
* R (DESeq2)
* Bash (pipeline automation)
* Git (version control) 

**Data**
Analyzed dataset: (GEO accession GSE63310) RNA-seq dataset featuring transcriptomic data from sorted mouse mammary epithelial cell populations. Raw count files were downloaded from the Arizona State University high performance computing cluster and transferred locally using FileZilla.

***Question***
How does gene expression differ between basal mammary epithelial cells and more differentiated luminal cell populations?

**Statistical Methods**
* Negative binomial model (DESeq2)
* Benjamini–Hochberg multiple testing correction
* Log2 fold-change shrinkage using apeglm
* Low-count gene filtering before modeling 

**Reproducibility**
1. Open Visual Studio Code
2.  Select Interpreter: `~/miniconda3/envs/bio-p2/python.exe`
3. Recreate environment:
```bash
conda create -n bio-p2 python=3.11
conda activate bio-p2
pip install -r requirements.txt 
```
4. One time insatll:
```bash
Rscript -e "if (!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('DESeq2','apeglm')); install.packages(c('ggplot2','ggrepel','pheatmap'), repos='https://cloud.r-project.org')"
```
8. Place selected GSM count files from GSE63310 into data/dev
9. From project root run full pipeline:
```bash
Rscript scripts/run_pipeline.sh
```
15. Final ouputs in results/ and data/processed/deseq2

## Citation
Log2 fold-change shrinkage performed using apeglm:
Zhu, A., Ibrahim, J.G., Love, M.I. (2018).
Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.
Bioinformatics.
https://doi.org/10.1093/bioinformatics/bty895

