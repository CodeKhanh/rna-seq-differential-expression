# Interpretation of Results

# Results Summary
This folder contains the final outputs from the RNA-seq differential expression analysis of GEO dataset GSE63310. 

# Research Question
How does gene expression differ between basal mammary epithelial cells and more differentiated luminal cell populations?

## Comparisons Performed
The following contrasts were analyzed using DESeq2:
- Luminal Progenitor (LP) vs Basal
- Mature Luminal (ML) vs Basal
Basal cells were used as the reference group.

## Differential Expression Tables
### 1. Basal_vs_LP_DEGs_full.csv
Contains shrunken log2 fold-change estimates and adjusted p-values for LP vs Basal.

### 2. Basal_vs_ML_DEGs_full.csv
Contains shrunken log2 fold-change estimates and adjusted p-values for ML vs Basal.  

Each table includes:
- `baseMean` (average normalized expression)
- `log2FoldChange`
- `lfcSE`
- `pvalue`
- `padj` (Benjamini–Hochberg adjusted p-value)

*Note: Gene identifiers in these tables are numeric (they appear to be Entrez Gene IDs), not gene symbols.*


## PCA.png 
Principal Component Analysis of variance-stabilized expression data.

**Purpose:** 
- Check whether samples cluster by biological condition  
- Detect possible outliers  
- Visualize overall differences between groups  

**Observations:**
- **PC1 (75.5% of variance)** clearly separates Basal cells from LP and ML samples.  
- Basal replicates cluster tightly. LP and ML replicates are more spread out along PC2, but each condition still forms its own group.  
- LP and ML samples separate further along **PC2 (15.6% of variance)**. 

**Relation to Biological Question** 
Most of the variation in the dataset is explained by differences between basal and luminal cell types. This shows that basal cells and differentiated luminal cells have clearly different overall gene expression profiles.

## MA Plots
### MA_Basal_vs_LP.png
### MA_Basal_vs_ML.png

**Purpose**  
- Show how gene expression changes between groups  
- Identify differentially expressed genes  
- Visualize fold changes across expression levels  
  
**Features**
- X-axis: mean expression
- Y-axis: shrunken log2 fold change between conditions
- Each point represents one gene.

**Observations**
- Many genes show strong positive or negative fold changes.  
- Differences occur across a wide range of expression levels.  
- The fold changes are not centered tightly at zero.  

**Relation to Biological Question** 
Gene expression differences between basal and luminal populations are widespread. The transition from basal to luminal states involves coordinated changes in many genes, not just a small number of markers.

## Heatmaps
### Heatmap_top50_Basal_vs_LP.png
### Heatmap_top50_Basal_vs_ML.png

These heatmaps display the top 50 genes ranked by absolute log2 fold change for each comparison.

**Features**
- Rows: genes
- Columns: samples
- Values: variance-stabilized expression (z-score scaled)
- Clustering performed on both genes and samples

**Purpose**
- Visualize expression patterns
- Confirm separation between cell populations

**Observations** 
- Basal samples cluster together.  
- LP and ML samples cluster separately.  
- Genes highly expressed in Basal tend to be lower in LP and ML.  
- Genes enriched in luminal cells tend to be lower in Basal.  

**Relation to Biological Question** 
The top differentially expressed genes clearly distinguish basal and luminal identities. This supports the conclusion that basal and differentiated luminal cells maintain distinct transcriptional programs.
 

## Limitations
- Sample size is limited (n = 3 per group).  
- Bulk RNA-seq cannot detect variation within individual cells.  
- Results are computational and not experimentally validated.  
- Only comparisons relative to Basal were performed.  
- Genes are identified by numeric IDs and have not been annotated with gene symbols or tested for pathway enrichment.
- No batch or sequencing-lane variable is included in the model (design is `~ condition`).


## Conclusion
This reproducible RNA-seq workflow identifies meaningful gene expression differences between Basal, Luminal Progenitor, and Mature Luminal mammary epithelial cells.

Specifically:
- Global expression patterns separate clearly by lineage.  
- Many genes show significant changes between basal and luminal states.  
- Expression differences reflect lineage-specific regulation during differentiation.


## Direction of fold changes
A positive `log2FoldChange` means higher expression in LP or ML than in Basal. The files are named `Basal_vs_X`, but the fold change is X relative to Basal.

## Summary (adjusted p < 0.05, 16,545 genes tested)
- LP vs Basal: 9,539 significant genes (4,669 higher in LP, 4,870 lower); 5,831 also have |shrunken log2FC| > 1
- ML vs Basal: 9,757 significant genes (4,789 higher in ML, 4,968 lower); 6,026 also have |shrunken log2FC| > 1
- 7,390 genes are significant in both comparisons
