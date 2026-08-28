#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC="$PROJ/data/intermediate_public"
DST="$PROJ/data/processed/integrated"
mkdir -p "$DST"

required=(RNA_top_variable.csv proteomics_olink.csv metabolomics.csv outcomes.csv)
for f in "${required[@]}"; do
  if [ ! -f "$SRC/$f" ]; then
    echo "Missing required intermediate file: $SRC/$f" >&2
    exit 1
  fi
  cp "$SRC/$f" "$DST/$f"
done

if [ -f "$SRC/clinical.csv" ]; then cp "$SRC/clinical.csv" "$DST/clinical.csv"; fi
mkdir -p "$PROJ/results/mosaic"
echo "Intermediate matrices prepared under data/processed/integrated/."
