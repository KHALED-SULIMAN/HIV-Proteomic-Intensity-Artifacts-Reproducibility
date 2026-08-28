# Obtaining the source data

This repository does not duplicate raw sequencing files. Public source data are
retrieved from the repositories below, and analysis-ready intermediate matrices
can be released separately under `data/intermediate_public/` for peer-review
reproducibility.

## 1. Transcriptomics: NCBI SRA

BioProject **PRJNA983231**. Save the SRA run table as:

```text
data/metadata/PRJNA983231_SraRunInfo.csv
```

Its `LibraryName` field (format `P20109_nnn`) is required for participant
linkage. Place FASTQ files under:

```text
data/raw/cocomo158_rna/fastq/
```

## 2. COCOMO public source repository

Clone the original public repository:

```bash
bash scripts/get_public_sources.sh
```

Source repository:
`https://github.com/neogilab/Immunometabolism_HIV`

The original Olink NPX source file used to document the proteomic provenance is:

```text
external/Immunometabolism_HIV/data/Olink/Q-01546_Neogi_NPX_2022-06-09.csv
```

The public metabolomics source file is:

```text
external/Immunometabolism_HIV/data/metabolomics/cocomo_norm_data_filt.csv
```

The preprocessing scripts use the source repository's processed/clinical tables
required to construct `proteomics_olink.csv`, `metabolomics.csv`, and
`clinical_data_clean_with_clusters_and_categories.csv`. The exact files consumed
by the current pipeline are documented in the scripts themselves.

## 3. Participant linkage

`scripts/preprocessing/04_build_patient_linkage.R` reconstructs:

```text
SRA Run -> LibraryName (P20109_nnn)
        -> clinical User field after removal of trailing _R[0-9]+
        -> COCOMO_ID
```

Retaining the `_R` suffix yields zero matches.

## 4. Intermediate-data release for reviewers

See `data/intermediate_public/README.md`. Do not upload raw FASTQ files to GitHub.
