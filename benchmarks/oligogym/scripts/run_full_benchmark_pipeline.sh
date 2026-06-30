#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 6 ]; then
  echo "Usage:"
  echo "  bash benchmarks/oligogym/scripts/run_full_benchmark_pipeline.sh \\"
  echo "    <raw_csv> <dataset_name> <endpoint> <label_type> <modality> <group_col>"
  exit 1
fi

RAW_CSV="$1"
DATASET_NAME="$2"
ENDPOINT="$3"
LABEL_TYPE="$4"
MODALITY="$5"
GROUP_COL="$6"

NORMALIZED="benchmarks/oligogym/processed/${DATASET_NAME}_normalized.csv"
SPLIT_OUT="benchmarks/oligogym/splits/${DATASET_NAME}_split.csv"
RESULT_PREFIX="benchmarks/oligogym/results/${DATASET_NAME}_baseline"
ANALYSIS_PREFIX="benchmarks/oligogym/analysis/${DATASET_NAME}_baseline"

echo "== Step 1: Normalize dataset =="
python benchmarks/oligogym/scripts/normalize_benchmark_dataset.py \
  "$RAW_CSV" \
  "$NORMALIZED" \
  "$DATASET_NAME" \
  "$ENDPOINT" \
  "$LABEL_TYPE" \
  "$MODALITY"

echo
echo "== Step 2: Generate leakage-aware splits =="
python benchmarks/oligogym/scripts/make_leakage_aware_splits.py \
  --input "$NORMALIZED" \
  --output "$SPLIT_OUT" \
  --group-col "$GROUP_COL"

echo
echo "== Step 3: Check split leakage =="
python benchmarks/oligogym/scripts/check_split_leakage.py \
  --input "$SPLIT_OUT" \
  --group-col "$GROUP_COL"

echo
echo "== Step 4: Run baseline benchmark =="
python benchmarks/oligogym/scripts/run_baseline_benchmark.py \
  --input "$SPLIT_OUT" \
  --output-prefix "$RESULT_PREFIX"

echo
echo "== Step 5: Analyze prediction errors =="
python benchmarks/oligogym/scripts/analyze_baseline_errors.py \
  --predictions "${RESULT_PREFIX}.predictions.csv" \
  --output-prefix "$ANALYSIS_PREFIX"

echo
echo "== Done =="
echo "Artifacts created under:"
echo "  benchmarks/oligogym/processed/"
echo "  benchmarks/oligogym/splits/"
echo "  benchmarks/oligogym/results/"
echo "  benchmarks/oligogym/analysis/"
