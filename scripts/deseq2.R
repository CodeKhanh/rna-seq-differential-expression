#!/usr/bin/env Rscript
# scripts/deseq2.R
# Minimal DESeq2 pipeline: read counts and metadata, run DESeq2,
# then save results and a few diagnostic plots.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(apeglm)   #LFC shrinkage (citation in README.md) 
  library(pheatmap) #heatmaps of top genes
})

outdir <- "data/processed/deseq2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Read inputs 
counts <- as.matrix(
  read.csv(
    "data/processed/counts_cleaned.csv",
    row.names = 1,
    check.names = FALSE
  )
)
meta <- read.csv("data/processed/metadata_cleaned.csv",
                 stringsAsFactors = TRUE)

# If you flagged some rows to exclude, keep only included samples.
if ("include" %in% colnames(meta)) {
  meta <- meta[meta$include == TRUE, , drop = FALSE]
}

# Subset counts to the samples listed in metadata (same order).
counts <- counts[, as.character(meta$sample), drop = FALSE]

# Build DESeq2 dataset 
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData  = meta,
  design   = ~ condition
)

# Prefilter low-count genes and run DESeq2.
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

# Write the default results table (unshrunken LFC).
res <- results(dds)
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered),
          file = file.path(outdir, "deseq2_results.csv"),
          row.names = TRUE)

# Shrink log2 fold changes with apeglm (better effect estimates).
# Prefer explicit contrast names if you have them:
# print(resultsNames(dds))
# resLFC <- lfcShrink(dds, coef = "condition_LP_vs_Basal", type = "apeglm")
resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")

write.csv(as.data.frame(resLFC[order(resLFC$padj), ]),
          file = file.path(outdir, "deseq2_results_shrunken.csv"),
          row.names = TRUE)

# PCA (variance-stabilizing transform) 
vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     meta)

png(filename = file.path(outdir, "pca.png"),
    width = 1200, height = 900, res = 150)

ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100, 1), "%")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100, 1), "%")) +
  theme_minimal()

dev.off()

# MA plot (shrunken LFC)
png(filename = file.path(outdir, "ma.png"),
    width = 1200, height = 900, res = 150)

plotMA(resLFC, main = "MA plot (shrunken LFC)", ylim = c(-5, 5))

dev.off() 

# Volcano plot (shrunken LFC)
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)
# protect against padj NA or 0
res_df$padj_plot <- ifelse(is.na(res_df$padj), 1, res_df$padj)
res_df$minusLog10Padj <- -log10(pmax(res_df$padj_plot, 1e-300))
# choose top labels by absolute LFC among significant genes (padj < 0.05) or overall top
top_label_n <- 20
top_genes <- with(res_df, head(gene[order(-abs(log2FoldChange))], top_label_n))

png(filename = file.path(outdir, "volcano.png"), width = 1200, height = 900, res = 150)
ggplot(res_df, aes(x = log2FoldChange, y = minusLog10Padj)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_minimal() +
  xlab("log2 fold change") + ylab("-log10(padj)") +
  geom_text_repel(data = subset(res_df, gene %in% top_genes),
                  aes(label = gene), size = 3, max.overlaps = 30)
dev.off()

# Heatmap: top N genes by absolute shrunken LFC
topN <- 50
res_for_rank <- as.data.frame(resLFC)
ranked_genes <- rownames(res_for_rank)[order(-abs(res_for_rank$log2FoldChange))]
top_genes_heat <- head(ranked_genes, topN)

# extract VST expression for the top genes, scale rows (z-score)
mat <- assay(vsd)[top_genes_heat, , drop = FALSE]
mat_z <- t(scale(t(mat)))   # scale rows

# annotation for columns from meta (must match column order of counts/vsd)
annotation_col <- data.frame(condition = meta$condition)
rownames(annotation_col) <- as.character(meta$sample)

# draw and save heatmap (heatmap uses VST-normalized values)
png(filename = file.path(outdir, "heatmap_top50.png"), width = 1400, height = 1000, res = 150)
pheatmap(mat_z,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         main = paste0("Top ", topN, " genes by |LFC|"))
dev.off()


# Write normalized counts (VST) 
norm_counts <- assay(vsd)
write.csv(as.data.frame(norm_counts),
          file = file.path(outdir, "normalized_counts_vst.csv"),
          row.names = TRUE)

# Session info and finish 
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "sessionInfo.txt"))

message("DESeq2 run complete. Outputs in: ", outdir)