# Data Description

## Dataset: GEO GSE63310
This project uses the publicly available RNA-seq dataset: GEO Accession: **GSE63310**

The dataset contains transcriptomic data from sorted mouse mammary epithelial cell populations.

## Experimental Groups
Samples are grouped into the following biological conditions:

### Basal (Reference / Control)
- GSM1545538  
- GSM1545539  
- GSM1545542  

### Luminal Progenitor (LP)
- GSM1545535  
- GSM1545541  
- GSM1545545  

### Mature Luminal (ML)
- GSM1545536  
- GSM1545540  
- GSM1545544  

Each main group contains **three replicates**.

## Additional Samples (Excluded)
The dataset also contains:
- GSM1545537 — Mammosphere  
- GSM1545543 — CommaDB  

These samples were excluded from the primary differential expression analysis to maintain clear group comparisons between Basal, LP, and ML populations.

## Data Availability

Raw GSM count files are not included in this repository due to size constraints and `.gitignore` settings.
To reproduce the analysis:
1. Download raw count files from GEO (GSE63310).
2. Place GSM count files into: data/dev

## Metadata 
A cleaned metadata is included in the repository under the processed folder. The file helps organizes the raw data allowing for reproducability in DESeq2 analysis.
