# Preprocessing scripts (00–06)

RNA-seq pipeline, participant-linkage reconstruction, and integrated-dataset
assembly. Run in order. Edit the `PROJ`/project-root path at the top of each
script to match your system before running.

| Script | Language | Purpose |
|---|---|---|
| `00_build_salmon_index.sh` | bash | Decoy-aware Salmon index (GENCODE v50 + GRCh38); 655,022 refs, 194 decoys, k=31 |
| `01_quantify_rna_salmon.sh` | bash | fastp trimming + Salmon selective-alignment quantification (mapping 85.9–92.9%) |
| `02_build_count_matrix.R` | R | tximport: transcripts → genes (82,323 genes); tx2gene parsed from GENCODE names |
| `03_normalize_qc.R` | R | Filter (≥10 counts in ≥25% samples → ~20,139 genes) + DESeq2 VST + PCA/correlation QC |
| `04_build_patient_linkage.R` | R | Reconstruct SRR → P20109 → COCOMO_ID linkage (remove `_R[0-9]+` suffix) |
| `05_assemble_integrated_dataset.R` | R | Align 4 layers to the 89-participant cohort; assert identical ordering |
| `06_build_outcomes_table.R` | R | Build clinical outcomes table (primary endpoint: metabolic syndrome) |

## Requirements
- Tools (conda): `salmon`, `fastp`, `fastqc`, `multiqc`, `pigz`
- R ≥ 4.3 with Bioconductor: `tximport`, `DESeq2`, plus `readr`

## Notes
- `00` and `01` use `GENCODE_RELEASE=50`; the two must match.
- `04` performs the linkage described in the manuscript Methods; retaining the
  `_R` suffix returns zero matches.
- `05` writes aligned matrices to `data/processed/integrated/`, which the
  analysis scripts in `../` consume.
