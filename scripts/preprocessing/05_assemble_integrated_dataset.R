#!/usr/bin/env Rscript
# =============================================================================
# 05_assemble_integrated_dataset.R
# Assemble the aligned four-layer dataset for the 89 linked participants.
# All layers are restricted to the linked COCOMO_IDs and ordered identically;
# identical ordering is asserted before writing.
#   RNA  -> top 3,000 most-variable VST genes
#   Proteomics (Olink, 2,923), Metabolomics (877), Clinical
# =============================================================================
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
SCRIPT_DIR <- if (length(file_arg) > 0) dirname(normalizePath(file_arg[1])) else getwd()
PROJ <- normalizePath(file.path(SCRIPT_DIR, "../.."), mustWork = FALSE)
OUT  <- file.path(PROJ, "data/processed")
INT  <- file.path(OUT, "integrated"); dir.create(INT, showWarnings = FALSE, recursive = TRUE)
COC  <- file.path(PROJ, "external/Immunometabolism_HIV/processing")

link <- read.csv(file.path(OUT, "master_patient_linkage.csv"), stringsAsFactors = FALSE)
ids  <- as.character(link$COCOMO_ID)

# --- RNA: VST matrix (genes x SRR) -> restrict, rename to COCOMO_ID ----------
vst  <- as.matrix(read.csv(file.path(OUT, "transcriptomics_vst.csv"), row.names = 1))
srr2coc <- setNames(as.character(link$COCOMO_ID), link$Run)
vst  <- vst[, colnames(vst) %in% link$Run, drop = FALSE]
colnames(vst) <- srr2coc[colnames(vst)]
# top 3,000 most variable genes
v <- apply(vst, 1, var); top <- names(sort(v, decreasing = TRUE))[1:3000]
rna <- t(vst[top, ids])                            # samples x genes

# --- Proteomics / metabolomics / clinical (indexed by COCOMO_ID) ------------
prot <- read.csv(file.path(COC, "COCOMO_proteomics_olink_filt.csv"), row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(COC, "metabolomics.csv"),      row.names = 1, check.names = FALSE)
clin <- read.csv(file.path(COC, "clinical_data_clean_with_clusters_and_categories.csv"),
                 stringsAsFactors = FALSE); rownames(clin) <- as.character(clin$COCOMO_ID)
prot <- prot[ids, , drop = FALSE]
meta <- meta[ids, , drop = FALSE]
clin <- clin[ids, , drop = FALSE]

# --- assert identical ordering across all layers ----------------------------
stopifnot(identical(rownames(rna), ids),
          identical(rownames(prot), ids),
          identical(rownames(meta), ids),
          identical(rownames(clin), ids))

write.csv(rna,  file.path(INT, "RNA_top_variable.csv"))
write.csv(prot, file.path(INT, "COCOMO_proteomics_olink_filt.csv"))
write.csv(meta, file.path(INT, "metabolomics.csv"))
write.csv(clin, file.path(INT, "clinical.csv"))
cat("Integrated cohort assembled: n =", length(ids), "participants, 4 layers\n")
