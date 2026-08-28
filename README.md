# Proteomic Intensity Artifacts in HIV Multi-Omics Integration

Code and per-figure source data for the manuscript:

> **Proteomic Intensity Artifacts in HIV Multi-Omics Integration: Deep Learning Benchmarking and Corrected Metabolic-Risk Endotyping**
> K. M. Elamin, S. M. I. Elbashir, I. Adam. Submitted to *Proteomes* (MDPI), 2026.

This repository contains the analysis pipeline, artifact-diagnostic framework, and figure-generation code needed to reproduce the study. Public source-data provenance is documented in `data/source_manifest.tsv`. For peer-review reproducibility, the repository is designed to support three entry points: full reproduction from public source data, analysis from released intermediate matrices, and figure-only reproduction. See `REPRODUCIBILITY.md`.

---

## Overview

The study assembles a linked four-layer multi-omics resource (transcriptomics, Olink proteomics, metabolomics, clinical) for 89 antiretroviral-treated people living with HIV, benchmarks a product-of-experts variational autoencoder against linear integration, and shows that the initially derived endotypes were driven by a global per-sample intensity shift in the proteomic layer. The repository provides:

- the RNA-seq preprocessing pipeline (raw reads to gene counts),
- the participant-linkage reconstruction that joins the public deposits,
- the integration models and the outcome-anchored benchmark,
- the four-test artifact-diagnostic framework and the normalization-strategy benchmark,
- endotype derivation, pathway enrichment, and drug mapping,
- the missing-modality experiment,
- MATLAB figure-generation scripts and per-panel source data.

---

## Repository structure

```
.
├── README.md
├── LICENSE                     # MIT (code)
├── LICENSE-DATA                # CC-BY-4.0 (derived figure data)
├── CITATION.cff                # how to cite this repository
├── REPRODUCIBILITY.md          # reviewer entry points
├── UPLOAD_TO_GITHUB.md         # upload instructions
├── requirements.txt            # Python environment
├── environment.yml             # conda environment (equivalent)
├── data/
│   ├── README.md               # how to obtain the public source data
│   ├── source_manifest.tsv     # exact data provenance
│   └── intermediate_public/    # reviewer-facing transformed matrices
├── scripts/
│   ├── preprocessing/          # 00–06: RNA-seq + linkage + assembly
│   │   ├── 00_build_salmon_index.sh
│   │   ├── 01_quantify_rna_salmon.sh
│   │   ├── 02_build_count_matrix.R
│   │   ├── 03_normalize_qc.R
│   │   ├── 04_build_patient_linkage.R
│   │   ├── 05_assemble_integrated_dataset.R
│   │   └── 06_build_outcomes_table.R
│   ├── mosaic_stage1_core.py            # PoE-VAE
│   ├── mosaic_stage2_anchor.py          # anchor sweep + baselines
│   ├── mosaic_stage2b_rigor.py          # bootstrap CIs
│   ├── mosaic_stage3_endotypes.py       # PCA vs MOFA-like, endotypes
│   ├── mosaic_stage4_characterize.py    # differential features
│   ├── mosaic_stage4b_artifact_check.py # 4-test artifact diagnostic
│   ├── mosaic_stage5_fix_and_verify.py  # correction + verification
│   ├── mosaic_stage6_pathways.py        # enrichment
│   ├── mosaic_stage7_drugs.py           # LINCS reversal (negative result)
│   ├── mosaic_stage7b_drugs_refined.py  # refined reversal
│   ├── mosaic_stage8_pathway_drugs.py   # DGIdb pathway-targeted (used)
│   ├── exp_A_missing_modality.py        # Figure 6
│   ├── exp_B_batch_correction.py        # Figure 7
│   └── export_figure_data.py            # writes per-panel CSVs
├── matlab/
│   ├── make_all_figures.m      # Figures 1b–7, S1
│   ├── make_figure1a.m         # study-design schematic
│   ├── fig_style.m             # centralized styling
│   └── save_panel.m            # robust .fig/.tiff/.png saver
└── results/
    └── figure_data/            # per-panel source CSVs (one per figure panel)
```

---

## Data availability

Raw sequencing files and the complete original source repository are not duplicated here. Public primary sources are cited below. For peer-review reproducibility, selected analysis-ready transformed matrices may be released in `data/intermediate_public/` after verification that redistribution is permitted and that no direct identifiers are present:

| Layer | Source | Accession |
|---|---|---|
| Transcriptomics (raw reads) | NCBI Sequence Read Archive | BioProject **PRJNA983231** |
| Proteomics (Olink), metabolomics, clinical | Public COCOMO analysis repository | see `data/README.md` |

The participant-linkage step (`scripts/preprocessing/04_build_patient_linkage.R`) reconstructs the mapping between sequencing runs and the cohort identifier. Exact source-data provenance is listed in `data/source_manifest.tsv`. Reviewers can reproduce the full workflow or start from the intermediate matrices described in `data/intermediate_public/README.md`.

`results/figure_data/` contains the derived, non-identifying per-panel values used to render each figure; these are released under CC-BY-4.0.

---

## Reproducing the analysis

### 1. Environment

```bash
# option A: conda
conda env create -f environment.yml
conda activate hiv_integration

# option B: pip
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

The supplied Conda environment includes Python 3.11, R 4.3, Salmon, fastp, FastQC, MultiQC, SRA Toolkit, `tximport`, and `DESeq2`, together with the Python analysis dependencies.

### 2. Obtain the public data

Follow `data/README.md` to download the SRA reads and the COCOMO proteomic, metabolomic and clinical matrices into `data/`.

### 3. Preprocess and assemble

```bash
cd scripts/preprocessing
bash   00_build_salmon_index.sh
bash   01_quantify_rna_salmon.sh
Rscript 02_build_count_matrix.R
Rscript 03_normalize_qc.R
Rscript 04_build_patient_linkage.R
Rscript 05_assemble_integrated_dataset.R
Rscript 06_build_outcomes_table.R
```

### 4. Run the analysis stages

```bash
cd ../..
python scripts/mosaic_stage2_anchor.py         # benchmark (Fig 2)
python scripts/mosaic_stage4b_artifact_check.py # artifact diagnostics (Fig 3)
python scripts/mosaic_stage5_fix_and_verify.py  # correction
python scripts/mosaic_stage6_pathways.py        # pathways (Fig 5a,b)
python scripts/mosaic_stage8_pathway_drugs.py   # drugs (Fig 5c,d)
python scripts/exp_A_missing_modality.py        # Fig 6
python scripts/exp_B_batch_correction.py        # Fig 7
python scripts/export_figure_data.py            # writes results/figure_data/
```

### 5. Generate figures

In MATLAB (R2025b or compatible):

```matlab
cd matlab
make_all_figures   % Figures 1b–7 and S1
make_figure1a      % study-design schematic
```

Each panel is saved as `.fig` (editable), `.tiff` and `.png` at 1200 dpi.

---

## Notes and caveats

- **Paths.** Project-root paths are resolved automatically from the repository location; no user-specific absolute path is required.
- **Drug candidates** are computational hypotheses only and are not clinical recommendations; antiretroviral-interaction annotations are a screening aid requiring pharmacological verification.
- **Batch correction** is empirical (per-sample), as no plate or batch metadata accompany the deposited proteomic matrix.
- **Determinism.** Random seeds are fixed, but exact numerical values may vary slightly across BLAS/library versions.

---

## Citing this work

If you use this code or the diagnostic framework, please cite the manuscript (details in `CITATION.cff`). A DOI for this repository can be minted via Zenodo on release.

## License

Code is released under the MIT License (`LICENSE`). Derived per-figure data in `results/figure_data/` are released under CC-BY-4.0 (`LICENSE-DATA`).

## Contact

Corresponding author: Prof. Ishag Adam — ia.ahmed@qu.edu.sa
Code author: K. M. Elamin — Kumamoto University.
