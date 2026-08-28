#!/usr/bin/env Rscript
# =============================================================================
# 06_build_outcomes_table.R
# Build the clinical outcomes table (indexed by COCOMO_ID, aligned to the
# integrated cohort). Primary held-out endpoint: metabolic syndrome (METS).
# =============================================================================
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
SCRIPT_DIR <- if (length(file_arg) > 0) dirname(normalizePath(file_arg[1])) else getwd()
PROJ <- normalizePath(file.path(SCRIPT_DIR, "../.."), mustWork = FALSE)
INT  <- file.path(PROJ, "data/processed/integrated")

clin <- read.csv(file.path(INT, "clinical.csv"), row.names = 1, check.names = FALSE)
ids  <- rownames(clin)

# Select outcome + covariate columns actually present in the clinical table.
wanted <- c("METS","BMI","AGE","Tgl","Hdl","Ldl","CD4","CD8",
            "central_obesity","hypertension","VAT","SAT","cluster")
present <- intersect(wanted, colnames(clin))
out <- data.frame(COCOMO_ID = ids, clin[, present, drop = FALSE],
                  check.names = FALSE, row.names = NULL)

write.csv(out, file.path(INT, "outcomes.csv"), row.names = FALSE)
cat("Outcomes table:", nrow(out), "participants,",
    length(present), "variables ->", file.path(INT, "outcomes.csv"), "\n")
