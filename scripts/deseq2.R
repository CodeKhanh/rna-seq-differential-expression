#!/usr/bin/env Rscript
# scripts/deseq2.R

# This script runs a full DESeq2 workflow:
# - Reads cleaned count data and metadata
# - Runs differential expression
# - Extracts result tables for Basal vs LP and Basal vs ML
# - Adjusts p-values (BH) on extracted tables
# - Produces diagnostic plots: PCA, MA, Heatmap

# Outputs are written to:
# - data/processed/deseq2/
# - results/

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(apeglm)    # recommended shrinkage method
  library(pheatmap)  # for heatmaps
})

outdir <- file.path("data", "processed", "deseq2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Load input data
# Counts matrix: genes as rows, samples as columns
counts <- as.matrix(
  read.csv(
    file.path("data", "processed", "counts_cleaned.csv"),
    row.names = 1,
    check.names = FALSE
  )
)

# Metadata: contains columns "sample" and "condition"
meta <- read.csv(
  file.path("data", "processed", "metadata_cleaned.csv"),
  stringsAsFactors = TRUE
)

# Drop samples marked include == FALSE
if ("include" %in% colnames(meta)) {
  meta <- meta[meta$include == TRUE, , drop = FALSE]
}

# Ensure sample order matches between counts and metadata
meta$sample <- as.character(meta$sample)
counts <- counts[, meta$sample, drop = FALSE]

# Ensure Basal is reference where present
meta$condition <- factor(meta$condition)
if ("Basal" %in% levels(meta$condition)) {
  meta$condition <- relevel(meta$condition, ref = "Basal")
}

# Build and run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

# Remove extremely low-count genes to reduce noise
dds <- dds[rowSums(counts(dds)) > 1, ]

# Fit model
dds <- DESeq(dds)

# Save object
saveRDS(dds, file = file.path(outdir, "dds_object.rds"))

# Contrasts and results
message("Available contrasts: ",
        paste(resultsNames(dds), collapse = ", "))

res_LP <- results(dds, contrast = c("condition", "LP", "Basal"))
res_ML <- results(dds, contrast = c("condition", "ML", "Basal"))

# Save unshrunken tables ordered by padj
write.csv(as.data.frame(res_LP[order(res_LP$padj), ]),
          file = file.path(outdir, "Basal_vs_LP_unshrunken.csv"))

write.csv(as.data.frame(res_ML[order(res_ML$padj), ]),
          file = file.path(outdir, "Basal_vs_ML_unshrunken.csv"))

# Shrink log2 fold changes using apeglm
res_LP_shr <- lfcShrink(dds,
                        coef = "condition_LP_vs_Basal",
                        type = "apeglm")

res_ML_shr <- lfcShrink(dds,
                        coef = "condition_ML_vs_Basal",
                        type = "apeglm")

write.csv(as.data.frame(res_LP_shr[order(res_LP_shr$padj), ]),
          file = file.path(outdir, "Basal_vs_LP_shrunken.csv"))

write.csv(as.data.frame(res_ML_shr[order(res_ML_shr$padj), ]),
          file = file.path(outdir, "Basal_vs_ML_shrunken.csv"))

# Copy clean versions to results for GitHub viewing
file.copy(file.path(outdir, "Basal_vs_LP_shrunken.csv"),
          file.path("results", "Basal_vs_LP_DEGs_full.csv"),
          overwrite = TRUE)

file.copy(file.path(outdir, "Basal_vs_ML_shrunken.csv"),
          file.path("results", "Basal_vs_ML_DEGs_full.csv"),
          overwrite = TRUE)

# Extract concise LFC + adjusted p-value
extract_cols <- function(res_obj, name_prefix) {
  df <- as.data.frame(res_obj)
  df$gene <- rownames(df)
  keep <- c("gene", "log2FoldChange", "lfcSE",
            "stat", "pvalue", "padj")
  df_out <- df[, intersect(keep, colnames(df)), drop = FALSE]

  # guarantee column order
  col_order <- intersect(keep, colnames(df_out))
  df_out <- df_out[, col_order]

  # Adjust p-values using BH on raw pvalue
  if ("pvalue" %in% colnames(df_out)) {
    df_out$padj_manual <- p.adjust(df_out$pvalue, method = "BH")
  } else {
    df_out$padj_manual <- NA
  }

  write.csv(df_out[order(ifelse(is.na(df_out$padj), 1, df_out$padj)), ],
            file = file.path(outdir,
                             paste0(name_prefix, "_lfc_padj_table.csv")),
            row.names = FALSE)
  return(df_out)
}

lp_table <- extract_cols(res_LP_shr, "Basal_vs_LP")
ml_table <- extract_cols(res_ML_shr, "Basal_vs_ML")

# Combined table side-by-side LFC and padj
common_genes <- union(lp_table$gene, ml_table$gene)
combined <- data.frame(gene = common_genes, stringsAsFactors = FALSE)

combined <- merge(combined,
                  lp_table[, c("gene", "log2FoldChange",
                               "padj", "padj_manual")],
                  by = "gene", all.x = TRUE)
colnames(combined)[colnames(combined) == "log2FoldChange"] <- "log2FC_LP"
colnames(combined)[colnames(combined) == "padj"] <- "padj_LP"
colnames(combined)[colnames(combined) == "padj_manual"] <- "padj_manual_LP"

combined <- merge(combined,
                  ml_table[, c("gene", "log2FoldChange",
                               "padj", "padj_manual")],
                  by = "gene", all.x = TRUE)
colnames(combined)[colnames(combined) == "log2FoldChange"] <- "log2FC_ML"
colnames(combined)[colnames(combined) == "padj"] <- "padj_ML"
colnames(combined)[colnames(combined) == "padj_manual"] <- "padj_manual_ML"

write.csv(combined[order(combined$padj_LP, combined$padj_ML), ],
          file = file.path(outdir,
                           "Combined_lfc_padj_Basal_vs_LP_vs_ML.csv"),
          row.names = FALSE)

# PCA
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file = file.path(outdir, "vsd_object.rds"))

pca <- prcomp(t(assay(vsd)))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2],
                     meta)

png(file.path(outdir, "PCA.png"),
    width = 1200, height = 900, res = 150)

ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100, 1), "%")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100, 1), "%")) +
  theme_minimal()

dev.off()

file.copy(file.path(outdir, "PCA.png"),
          file.path("results", "PCA.png"),
          overwrite = TRUE)

# MA plots: LP vs Basal and ML vs Basal
png(file.path(outdir, "MA_Basal_vs_LP.png"),
    width = 1200, height = 900, res = 150)
plotMA(res_LP_shr, main = "MA: Basal vs LP (shrunken)", ylim = c(-5, 5))
dev.off()
file.copy(file.path(outdir, "MA_Basal_vs_LP.png"),
          file.path("results", "MA_Basal_vs_LP.png"),
          overwrite = TRUE)

png(file.path(outdir, "MA_Basal_vs_ML.png"),
    width = 1200, height = 900, res = 150)
plotMA(res_ML_shr, main = "MA: Basal vs ML (shrunken)", ylim = c(-5, 5))
dev.off()
file.copy(file.path(outdir, "MA_Basal_vs_ML.png"),
          file.path("results", "MA_Basal_vs_ML.png"),
          overwrite = TRUE)

# Heatmaps: per-contrast top N genes by absolute LFC
make_heatmap <- function(res_shr, prefix, topN = 50) {
  res_for_rank <- as.data.frame(res_shr)
  ranked_genes <- rownames(res_for_rank)[
    order(-abs(res_for_rank$log2FoldChange))
  ]
  top_genes_heat <- head(ranked_genes, topN)

  # extract VST expression for the top genes, scale rows (z-score)
  mat <- assay(vsd)[top_genes_heat, , drop = FALSE]
  mat_z <- t(scale(t(mat)))   # scale rows

  # annotation for columns from meta
  annotation_col <- data.frame(condition = meta$condition)
  rownames(annotation_col) <- as.character(meta$sample)

  out_png <- file.path(outdir, paste0("Heatmap_top", topN, "_", prefix, ".png"))
  png(filename = out_png, width = 1400, height = 1000, res = 150)
  pheatmap(mat_z,
           annotation_col = annotation_col,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           fontsize_col = 10,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "none",
           main = paste0("Top ", topN, " genes by |LFC| (", prefix, ")"))
  dev.off()

  file.copy(out_png,
            file.path("results", paste0("Heatmap_top", topN, "_", prefix, ".png")),
            overwrite = TRUE)
}

# generate heatmaps for LP and ML
make_heatmap(res_LP_shr, "Basal_vs_LP", topN = 50)
make_heatmap(res_ML_shr, "Basal_vs_ML", topN = 50)

# Save normalized counts + session info
write.csv(as.data.frame(assay(vsd)),
          file = file.path(outdir, "normalized_counts_vst.csv"))

writeLines(capture.output(sessionInfo()),
           file.path(outdir, "sessionInfo.txt"))

message("DESeq2 analysis complete.")
message("Technical outputs: ", outdir)
message("Clean results copied to: results/")