#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$SCRIPT_DIR/.." && pwd)"
mkdir -p "$PROJ/external" "$PROJ/data/metadata" "$PROJ/data/raw/cocomo158_rna/fastq"

if [ ! -d "$PROJ/external/Immunometabolism_HIV/.git" ]; then
  git clone https://github.com/neogilab/Immunometabolism_HIV.git \
    "$PROJ/external/Immunometabolism_HIV"
else
  echo "COCOMO source repository already present."
fi

cat <<'MSG'
Next, obtain the NCBI SRA run table for BioProject PRJNA983231 and save it as:
  data/metadata/PRJNA983231_SraRunInfo.csv
Then download the RNA-seq FASTQ files into:
  data/raw/cocomo158_rna/fastq/
See data/README.md for details.
MSG
