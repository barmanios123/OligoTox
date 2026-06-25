#!/usr/bin/env bash
set -euo pipefail
cd ~/oligotox

if [ -f scripts/papargyri_extractor.py ] && [ -f papargyri_2020_LNA_512gapmers_TableS1.xlsx ]; then
  python scripts/papargyri_extractor.py \
    --input papargyri_2020_LNA_512gapmers_TableS1.xlsx \
    --output data/curated/papargyri_512gapmers_labeled.csv
else
  echo "Papargyri extractor or Table S1 workbook not present; skipping Papargyri ingest."
fi

if [ -f scripts/fda_nda_platelet_extractor.py ]; then
  echo "FDA extractor present. Review and update seeded values before running full extraction."
else
  echo "FDA platelet extractor not present."
fi
