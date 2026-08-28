# Analysis-ready intermediate data for peer-review reproducibility

The Academic Editor requested transformed/intermediate data so that reviewers
can begin reproduction from an intermediate stage rather than repeating the
entire raw-data preprocessing workflow.

Before public release, inspect every table and confirm that it contains no direct
identifiers and that redistribution is permitted by the source-data license and
ethics conditions.

Recommended files to place here from the completed local analysis:

- `master_patient_linkage.csv` — pseudonymous linkage used by the analysis.
- `RNA_top_variable.csv` — 89 x 3,000 variance-stabilized transcriptomic matrix.
- `proteomics_olink.csv` — aligned Olink protein-target abundance matrix.
- `metabolomics.csv` — aligned metabolomic matrix.
- `outcomes.csv` — minimal clinical endpoints/covariates required by the scripts.

Do not upload raw FASTQ files here. The raw RNA-seq data remain available from
NCBI SRA under BioProject PRJNA983231.
