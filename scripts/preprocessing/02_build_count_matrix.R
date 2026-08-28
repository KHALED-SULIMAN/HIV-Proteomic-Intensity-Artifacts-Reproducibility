#!/usr/bin/env Rscript
# =============================================================================
# 02_build_count_matrix.R
# Build a gene-level count matrix from all Salmon quant.sf files via tximport.
# The transcript->gene map is parsed from GENCODE-style names (split on "|"),
# so no GTF is required. Output: 82,323 genes x N samples.
# =============================================================================
suppressMessages({ library(tximport); library(readr) })

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
SCRIPT_DIR <- if (length(file_arg) > 0) dirname(normalizePath(file_arg[1])) else getwd()
PROJ <- normalizePath(file.path(SCRIPT_DIR, "../.."), mustWork = FALSE)
QUANT <- file.path(PROJ, "data/interim/rna_quantification")
OUT   <- file.path(PROJ, "data/processed")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 1. Find every sample that has a quant.sf
files  <- list.files(QUANT, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
sample <- basename(dirname(files)); names(files) <- sample
cat("Found", length(files), "samples\n")

# 2. Build transcript -> gene map from GENCODE names
#    (format: ENST...|ENSG...|...|SYMBOL|... ; ENSG is field 2)
first   <- read_tsv(files[1], show_col_types = FALSE)
ids     <- do.call(rbind, strsplit(first$Name, "\\|", perl = TRUE))
tx2gene <- data.frame(TXNAME = first$Name, GENEID = ids[, 2])

# 3. Import and sum transcripts to genes
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE)

# 4. Save matrices
counts <- round(txi$counts)
write.csv(counts,        file.path(OUT, "transcriptomics_counts.csv"))
write.csv(txi$abundance, file.path(OUT, "transcriptomics_tpm.csv"))
saveRDS(txi,             file.path(OUT, "transcriptomics_txi.rds"))

cat("Genes:", nrow(counts), " Samples:", ncol(counts), "\n")
cat("Saved to", OUT, "\n")
