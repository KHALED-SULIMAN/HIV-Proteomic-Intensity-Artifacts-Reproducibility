#!/usr/bin/env bash
# =============================================================================
# 00_build_salmon_index.sh
# Build a decoy-aware Salmon index for human RNA-seq quantification. RUN ONCE.
#
# Salmon  = transcript quantification by selective alignment (light on memory)
# GENCODE = human reference transcriptome + annotation used for alignment
# decoy   = the genome, added so reads from unannotated regions are not
#           forced onto a transcript by mistake
#
# Requirements: ~16-25 GB RAM for the index build, ~20 GB free disk.
# One-time tool environment:
#   conda create -n hiv_rnaseq -c conda-forge -c bioconda -y \
#       salmon fastp fastqc multiqc pigz
#   conda activate hiv_rnaseq
# =============================================================================
set -euo pipefail

# ---- Paths (edit PROJ to your project root) ---------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$SCRIPT_DIR/../.." && pwd)"
REF="$PROJ/data/raw/reference"
mkdir -p "$REF"; cd "$REF"

# ---- GENCODE release (this study used release 50) ---------------------------
# Confirm the current human release at https://www.gencodegenes.org/human/
GENCODE_RELEASE=50

BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
TRANSCRIPTS="gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"
GENOME="GRCh38.primary_assembly.genome.fa.gz"

# ---- Download reference files (only if missing) -----------------------------
[ -f "$TRANSCRIPTS" ] || wget -c "${BASE}/${TRANSCRIPTS}"
[ -f "$GENOME" ]      || wget -c "${BASE}/${GENOME}"

# ---- Build decoy list, concatenate (transcripts first), build index ---------
zcat "$GENOME" | grep "^>" | cut -d ' ' -f1 | sed 's/>//g' > decoys.txt
cat "$TRANSCRIPTS" "$GENOME" > gentrome.fa.gz

salmon index --threads "$(nproc)" \
  --transcripts gentrome.fa.gz --decoys decoys.txt \
  --index "$REF/salmon_index_gencode_v${GENCODE_RELEASE}" --kmerLen 31

echo "Salmon index built at: $REF/salmon_index_gencode_v${GENCODE_RELEASE}"
# Verified build: 655,022 references, 194 decoys, k=31.
