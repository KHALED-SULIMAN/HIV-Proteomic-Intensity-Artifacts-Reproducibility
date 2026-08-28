#!/usr/bin/env bash
# =============================================================================
# 01_quantify_rna_salmon.sh
# Trim (fastp) and quantify (Salmon) every paired-end RNA sample.
# Safe to run mid-download; already-quantified samples are skipped.
#
# fastp = adapter/quality trimming ; salmon quant = per-transcript abundance
# Observed mapping rates in this study: 85.9-92.9% (median ~89.5%).
# =============================================================================
set -euo pipefail
#   conda activate hiv_rnaseq

# ---- Paths ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$SCRIPT_DIR/../.." && pwd)"
FASTQ="$PROJ/data/raw/cocomo158_rna/fastq"
TRIM="$PROJ/data/interim/rna_qc/trimmed"
QUANT="$PROJ/data/interim/rna_quantification"
GENCODE_RELEASE=50   # must match step 00
INDEX="$PROJ/data/raw/reference/salmon_index_gencode_v${GENCODE_RELEASE}"

mkdir -p "$TRIM" "$QUANT"; THREADS="$(nproc)"

# ---- Loop over every *_1.fastq.gz and find its mate *_2.fastq.gz ------------
for R1 in "$FASTQ"/*_1.fastq.gz; do
    SAMPLE="$(basename "$R1" _1.fastq.gz)"      # e.g. SRR24906078
    R2="$FASTQ/${SAMPLE}_2.fastq.gz"
    [ -f "$R2" ] || { echo "WARNING: mate missing for $SAMPLE, skipping"; continue; }
    [ -f "$QUANT/$SAMPLE/quant.sf" ] && { echo "Already quantified: $SAMPLE, skipping"; continue; }

    echo "=== $SAMPLE : trimming ==="
    fastp --in1 "$R1" --in2 "$R2" \
      --out1 "$TRIM/${SAMPLE}_1.trim.fastq.gz" --out2 "$TRIM/${SAMPLE}_2.trim.fastq.gz" \
      --json "$TRIM/${SAMPLE}.fastp.json" --html "$TRIM/${SAMPLE}.fastp.html" \
      --thread "$THREADS"

    echo "=== $SAMPLE : quantifying ==="
    salmon quant --index "$INDEX" --libType A \
      --mates1 "$TRIM/${SAMPLE}_1.trim.fastq.gz" --mates2 "$TRIM/${SAMPLE}_2.trim.fastq.gz" \
      --validateMappings --gcBias --seqBias --threads "$THREADS" \
      --output "$QUANT/$SAMPLE"
done
echo "Done. Per-sample results in: $QUANT/<SAMPLE>/quant.sf"
