#!/usr/bin/env Rscript

# =============================================================================
# 04_build_patient_linkage.R
#
# Reconstruct:
#   SRA Run -> LibraryName (P20109_nnn) -> COCOMO_ID
#
# Only successfully quantified RNA-seq runs are retained.
# Expected:
#   96 deposited SRA runs
#   95 successfully quantified
#   89 successfully quantified + clinically linked participants
# =============================================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])

SCRIPT_DIR <- if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg[1]))
} else {
    getwd()
}

PROJ <- normalizePath(file.path(SCRIPT_DIR, "../.."), mustWork = FALSE)

META <- file.path(PROJ, "data/metadata")
COC  <- file.path(PROJ, "external/Immunometabolism_HIV/processing")
OUT  <- file.path(PROJ, "data/processed")
QUANT <- file.path(PROJ, "data/interim/rna_quantification")
MAP_FILE <- file.path(PROJ, "results/figure_data/S1a_mapping_rates.csv")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# 1. SRA metadata
# -------------------------------------------------------------------------

sra <- read.csv(
    file.path(META, "PRJNA983231_SraRunInfo.csv"),
    stringsAsFactors = FALSE
)

stopifnot(all(c("Run", "LibraryName") %in% colnames(sra)))

sra <- sra[, c("Run", "LibraryName")]
sra$Run <- trimws(as.character(sra$Run))
sra$LibraryName <- trimws(as.character(sra$LibraryName))
sra$P20109 <- sra$LibraryName

cat("Deposited SRA runs:", length(unique(sra$Run)), "\n")

# -------------------------------------------------------------------------
# 2. Determine successfully quantified RNA runs
#
# Preferred source:
#   actual Salmon quant.sf files
#
# Fallback:
#   saved mapping-rate table from the completed analysis
# -------------------------------------------------------------------------

quant_files <- character(0)

if (dir.exists(QUANT)) {
    quant_files <- list.files(
        QUANT,
        pattern = "quant\\.sf$",
        recursive = TRUE,
        full.names = TRUE
    )
}

if (length(quant_files) > 0) {

    quantified_runs <- unique(basename(dirname(quant_files)))

    cat(
        "Quantified runs identified from Salmon outputs:",
        length(quantified_runs),
        "\n"
    )

} else {

    if (!file.exists(MAP_FILE)) {
        stop(
            "No Salmon quant.sf files and no S1a_mapping_rates.csv were found."
        )
    }

    mapping <- read.csv(MAP_FILE, stringsAsFactors = FALSE)

    if (!"sample" %in% colnames(mapping)) {
        stop("S1a_mapping_rates.csv lacks the required 'sample' column.")
    }

    quantified_runs <- unique(trimws(as.character(mapping$sample)))

    cat(
        "Quantified runs identified from saved mapping-rate table:",
        length(quantified_runs),
        "\n"
    )
}

if (length(quantified_runs) != 95) {
    stop(
        paste(
            "Expected 95 successfully quantified RNA runs but found",
            length(quantified_runs)
        )
    )
}

# Restrict SRA metadata to successfully quantified runs
sra_q <- sra[sra$Run %in% quantified_runs, , drop = FALSE]

# -------------------------------------------------------------------------
# 3. Clinical participant identifiers
# -------------------------------------------------------------------------

clin <- read.csv(
    file.path(
        COC,
        "clinical_data_clean_with_clusters_and_categories.csv"
    ),
    stringsAsFactors = FALSE
)

stopifnot(all(c("COCOMO_ID", "cluster", "User") %in% colnames(clin)))

clin$P20109 <- sub(
    "_R[0-9]+$",
    "",
    as.character(clin$User)
)

# Verify unique join keys
if (anyDuplicated(clin$P20109)) {
    stop("Clinical P20109 identifiers are not unique.")
}

# -------------------------------------------------------------------------
# 4. Join quantified RNA runs to COCOMO participants
# -------------------------------------------------------------------------

link <- merge(
    sra_q,
    clin[, c("COCOMO_ID", "cluster", "P20109")],
    by = "P20109",
    all = FALSE,
    sort = FALSE
)

link <- link[, c(
    "Run",
    "P20109",
    "COCOMO_ID",
    "cluster"
)]

cat("Successfully quantified runs:", nrow(sra_q), "\n")
cat("Linked participants:", nrow(link), "\n")
cat("Unique COCOMO_ID:", length(unique(link$COCOMO_ID)), "\n")

if (nrow(link) != 89) {
    stop(
        paste(
            "Expected 89 successfully quantified and linked participants;",
            "obtained",
            nrow(link)
        )
    )
}

if (length(unique(link$COCOMO_ID)) != 89) {
    stop("COCOMO_ID values are not unique.")
}

# -------------------------------------------------------------------------
# 5. Save
# -------------------------------------------------------------------------

write.csv(
    link,
    file.path(OUT, "master_patient_linkage.csv"),
    row.names = FALSE
)

cat(
    "Saved:",
    file.path(OUT, "master_patient_linkage.csv"),
    "\n"
)
