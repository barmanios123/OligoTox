#!/usr/bin/env bash
set -euo pipefail

cd ~/oligotox

echo '== OligoTox v0.3 package integration =='

SRC_DIR="downloads"
[ -d "$SRC_DIR" ] || { echo "Missing $SRC_DIR/"; exit 1; }

mkdir -p data/curated
mkdir -p data/raw/fda_reviews
mkdir -p scripts
mkdir -p docs
mkdir -p notebooks
mkdir -p reports

copy_if_exists () {
  local src="$1"
  local dst="$2"
  if [ -f "$src" ]; then
    cp "$src" "$dst"
    echo "copied: $src -> $dst"
  else
    echo "skip: $src not found"
  fi
}

copy_if_exists "$SRC_DIR/curated_oligo_tox_master_v3.csv" data/curated/curated_oligo_tox_master_v3.csv
copy_if_exists "$SRC_DIR/endpoint_summary_v3.csv" reports/endpoint_summary_v3.csv
copy_if_exists "$SRC_DIR/oligotox_labeling_scheme.csv" docs/oligotox_labeling_scheme_v3.csv
copy_if_exists "$SRC_DIR/MODEL_STATUS_v3.md" docs/MODEL_STATUS_v3.md
copy_if_exists "$SRC_DIR/nih_oligotox_challenge_registration_email.txt" docs/nih_oligotox_challenge_registration_email.txt
copy_if_exists "$SRC_DIR/papargyri_extractor.py" scripts/papargyri_extractor.py
copy_if_exists "$SRC_DIR/fda_nda_platelet_extractor.py" scripts/fda_nda_platelet_extractor.py
copy_if_exists "$SRC_DIR/train_v3.py" scripts/train_v3.py
copy_if_exists "$SRC_DIR/endpoint_gap_analysis.ipynb" notebooks/endpoint_gap_analysis.ipynb
copy_if_exists "$SRC_DIR/curated_oligo_tox_template.csv" data/curated/curated_oligo_tox_template_v3.csv
copy_if_exists "$SRC_DIR/download_all_datasets.sh" scripts/download_all_datasets.sh
copy_if_exists "$SRC_DIR/download_datasets_v3.sh" scripts/download_datasets_v3.sh

copy_if_exists "$SRC_DIR/inotersen_NDA211172_PharmR.pdf" data/raw/fda_reviews/inotersen_NDA211172_PharmR.pdf
copy_if_exists "$SRC_DIR/mipomersen_NDA203568_PharmR.pdf" data/raw/fda_reviews/mipomersen_NDA203568_PharmR.pdf
copy_if_exists "$SRC_DIR/nusinersen_NDA209531_PharmR.pdf" data/raw/fda_reviews/nusinersen_NDA209531_PharmR.pdf

cat > docs/V03_PACKAGE_INTEGRATION.md <<'DOC'
# v0.3 package integration

This repo integrated the downloaded v0.3 package from `downloads/` into the canonical project structure.

## Canonical destinations

- Curated dataset: `data/curated/curated_oligo_tox_master_v3.csv`
- FDA review PDFs: `data/raw/fda_reviews/`
- Extractors and training helpers: `scripts/`
- Labeling docs and model notes: `docs/`
- Gap notebook: `notebooks/`
- Endpoint summary report: `reports/endpoint_summary_v3.csv`

## Recommended next actions

1. Compare v3 schema with the current master dataset before overwriting any existing canonical master file.
2. Run `scripts/papargyri_extractor.py` once Table S1 Excel is downloaded locally.
3. Update `scripts/fda_nda_platelet_extractor.py` with real platelet nadirs from the NDA PDFs.
4. Decide whether `scripts/train_v3.py` becomes the canonical trainer or is merged into `train.py` / `pipelines/`.
5. Regenerate model artifacts only after the unified curated dataset is finalized.
DOC

cat > scripts/run_v03_ingest.sh <<'RUN'
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
RUN
chmod +x scripts/run_v03_ingest.sh

cat > scripts/train_from_v03.sh <<'RUN'
#!/usr/bin/env bash
set -euo pipefail
cd ~/oligotox

if [ -f scripts/train_v3.py ]; then
  python scripts/train_v3.py
else
  echo "scripts/train_v3.py not found"
  exit 1
fi
RUN
chmod +x scripts/train_from_v03.sh

echo '== Done =='
echo 'Created docs/V03_PACKAGE_INTEGRATION.md and helper scripts.'
