# Reproducibility entry points

Reviewers can reproduce the study at three levels.

## A. Full workflow from public source data

1. Clone the source COCOMO repository with `bash scripts/get_public_sources.sh`.
2. Download the NCBI SRA run table and FASTQ files for PRJNA983231.
3. Run `scripts/preprocessing/00_...` through `06_...` in order.
4. Run the analysis-stage scripts.
5. Regenerate figure source data and figures.

## B. Start from analysis-ready intermediate matrices

Place the released intermediate matrices in `data/processed/integrated/` or copy
the contents of `data/intermediate_public/` there, then run the stage scripts.
This entry point avoids re-running RNA-seq quantification while preserving the
exact inputs used for integration, artifact diagnostics, normalization, endotype
analysis, pathway analysis, and drug mapping.

## C. Reproduce figures only

Use the CSV files under `results/figure_data/` with the MATLAB scripts in
`matlab/` to regenerate the manuscript figures.

## Determinism

Random seeds are fixed in the analysis scripts. Small numerical differences can
still occur across BLAS, PyTorch, scikit-learn, and operating-system versions.
