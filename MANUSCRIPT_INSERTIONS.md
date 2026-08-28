# Exact manuscript insertions for Proteomes editorial pre-check

Replace `[GITHUB_URL]` only after the repository is public.
If you later archive a release in Zenodo, add the optional DOI sentence.

## Replace the entire Section 2.1 with

### 2.1. Data Sources and Cohort Assembly

All data analyzed in this study were derived from previously generated COCOMO cohort resources, and no new human data were generated. Raw RNA-sequencing reads were retrieved from the National Center for Biotechnology Information Sequence Read Archive (SRA) under BioProject PRJNA983231. The proteomic layer consisted of plasma proteins measured using the Olink Explore 3072 proximity extension assay platform. The original Olink NPX source file, `Q-01546_Neogi_NPX_2022-06-09.csv`, together with the corresponding assay documentation, was obtained from the public `neogilab/Immunometabolism_HIV` repository (`data/Olink/`). The untargeted metabolomic source matrix, `cocomo_norm_data_filt.csv`, was obtained from `data/metabolomics/` in the same repository. The processed clinical and molecular tables used for participant linkage and multi-omics assembly, including `clinical_data_clean_with_clusters_and_categories.csv`, were obtained from the same COCOMO analysis resource and are documented in the study reproducibility repository at [GITHUB_URL].

Proteins may exist as multiple proteoforms arising from sequence variation, alternative splicing, and post-translational modification. In the present study, the Olink proximity extension assays were interpreted as measurements of relative abundance at the protein-target level; individual proteoforms were not resolved, and proteomic differences are therefore described as differences in protein abundance rather than protein expression.

As the SRA records lacked participant identifiers, transcriptomic samples were linked to the remaining layers via the sequencing facility library code (SRA `LibraryName` field, format `P20109_nnn`), which was matched to the `User` field of the clinical table after removal of the trailing run suffix (`R[0-9]+`) to recover the participant identifier (`COCOMO_ID`) indexing the proteomic and metabolomic matrices. Of 96 sequencing runs, 95 were quantified, and 89 were linked to complete proteomic, metabolomic, and clinical data, defining the analysis cohort (n = 89).

## Replace the entire Section 2.11 with

### 2.11. Software, Code, and Reproducibility

Analyses were performed using Python 3.11 (PyTorch, scikit-learn, pandas, NumPy, SciPy, and GSEApy) and R 4.3 (tximport and DESeq2). The complete reproducibility repository for this study is publicly available at [GITHUB_URL]. The repository contains a source-data manifest, scripts for retrieving the public source resources, RNA-sequencing preprocessing and participant-linkage code, analysis-ready intermediate matrices for the 89 linked participants, model configurations and fixed random seeds, scripts for multi-omics integration and benchmarking, proteomic artifact diagnostics, normalization analyses, endotype derivation, pathway and pharmacological analyses, missing-modality experiments, and figure-level source data. The repository is organized so that the study can be reproduced either from the original public source data, from the released analysis-ready intermediate matrices, or from the figure-level source data for figure regeneration. Exact software dependencies and the execution order are documented in the repository `README.md`, `environment.yml`, and `REPRODUCIBILITY.md` files.

Optional after Zenodo archiving: `A frozen version of the reproducibility package is archived at Zenodo under DOI: [ZENODO_DOI].`

## Replace the entire Data Availability Statement with

### Data Availability Statement

The RNA-sequencing data analyzed in this study are publicly available through the National Center for Biotechnology Information Sequence Read Archive under BioProject PRJNA983231. The Olink Explore 3072 proteomic source file (`data/Olink/Q-01546_Neogi_NPX_2022-06-09.csv`), the untargeted metabolomic source matrix (`data/metabolomics/cocomo_norm_data_filt.csv`), and the COCOMO clinical/processing resources used for participant linkage are publicly available through the original `neogilab/Immunometabolism_HIV` repository. The study-specific reproducibility repository, including the source-data manifest, participant-linkage workflow, analysis-ready transformed matrices, analysis scripts, model settings and random seeds, statistical outputs, and figure-level source data, is publicly available at [GITHUB_URL]. Raw RNA-sequencing files are not duplicated in the study repository because they remain accessible through the Sequence Read Archive.

Optional after Zenodo archiving: `A frozen archival release is available at Zenodo under DOI: [ZENODO_DOI].`
