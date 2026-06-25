#!/usr/bin/env bash
set -euo pipefail

echo "=== OligoTox PMO cycle ==="
echo "CWD: $(pwd)"
echo

echo "[1/3] Rebuilding PMO training-ready + endpoint summary..."
PYTHONPATH=$HOME python scripts/curation/build_pmo_label_dataset.py
echo

echo "[2/3] Merging PMO batch into master curated set..."
PYTHONPATH=$HOME python scripts/curation/merge_into_master.py
echo

echo "[3/3] Running full training pipeline..."
PYTHONPATH=$HOME python train.py
echo
echo "=== Done ==="
