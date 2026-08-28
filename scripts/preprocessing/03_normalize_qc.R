#!/usr/bin/env Rscript
# =============================================================================
# 03_normalize_qc.R
# Filter low-count genes, apply DESeq2 variance-stabilizing transformation,
# and run sample-level QC (PCA, correlation) on the gene count matrix.
#   Filter: >=10 counts in >=25% of samples  (82,323 -> ~20,139 genes)
#   Transform: DESeq2 VST (blind = TRUE)
# =============================================================================
suppressMessages({ library(DESeq2) })

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
SCRIPT_DIR <- if (length(file_arg) > 0) dirname(normalizePath(file_arg[1])) else getwd()
PROJ <- normalizePath(file.path(SCRIPT_DIR, "../.."), mustWork = FALSE)
OUT  <- file.path(PROJ, "data/processed")
QC   <- file.path(PROJ, "data/interim/rna_qc"); dir.create(QC, showWarnings = FALSE, recursive = TRUE)

counts <- as.matrix(read.csv(file.path(OUT, "transcriptomics_counts.csv"), row.names = 1))
cat("Loaded counts:", nrow(counts), "genes x", ncol(counts), "samples\n")

# 1. Filter low-count genes: >=10 counts in >=25% of samples
keep <- rowSums(counts >= 10) >= (0.25 * ncol(counts))
counts_f <- counts[keep, ]
cat("Retained", nrow(counts_f), "genes after filtering\n")

# 2. Variance-stabilizing transformation (blind)
dds <- DESeqDataSetFromMatrix(countData = counts_f,
        colData = data.frame(row.names = colnames(counts_f),
                             sample = colnames(counts_f)),
        design = ~ 1)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_mat <- assay(vsd)
write.csv(vst_mat, file.path(OUT, "transcriptomics_vst.csv"))
cat("Saved VST matrix\n")

# 3. QC: PCA on the transformed matrix
pca <- prcomp(t(vst_mat), scale. = FALSE)
pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
cat(sprintf("PC1 %.1f%%  PC2 %.1f%%\n", pct[1], pct[2]))
write.csv(data.frame(sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2]),
          file.path(QC, "rna_pca.csv"), row.names = FALSE)

# 4. QC: sample-sample Spearman correlation
corr <- cor(vst_mat, method = "spearman")
write.csv(corr, file.path(QC, "rna_sample_correlation.csv"))
cat("QC written to", QC, "\n")
