#!/usr/bin/env bash
set -Eeuo pipefail

pass() { echo "[PASS] $*"; }
warn() { echo "[WARN] $*"; }
section() { echo; echo "## $*"; }

echo "== OligoTox project health =="
echo "[INFO] root: $(pwd)"

section "1) core repo layout"
for p in pyproject.toml README.md src data artifacts configs validate_platform.py; do
  [[ -e "$p" ]] && pass "$p exists" || warn "$p missing"
done

section "2) python environment"
python -V
which python
python - <<'PY'
mods = ["pandas","numpy","sklearn","joblib"]
for m in mods:
    try:
        __import__(m)
        print(f"[PASS] import {m}")
    except Exception as e:
        print(f"[FAIL] import {m}: {e}")
PY

section "3) package import smoke test"
python - <<'PY'
import sys, os
sys.path.insert(0, os.path.abspath("src"))
for t in ["oligotox", "oligotox.training.train"]:
    try:
        __import__(t)
        print(f"[PASS] imported {t}")
    except Exception as e:
        print(f"[WARN] import {t} failed: {e}")
PY

section "4) training input smoke test"
python - <<'PY'
from pathlib import Path
for p in [
    "data/processed/unified",
    "data/processed/trainingdataset.csv",
    "data/processed/datasetmetadata.json",
    "artifacts/config/composite_config.json",
    "artifacts/config/endpoint_thresholds.json",
]:
    pp = Path(p)
    print(f"[CHECK] {p} exists: {pp.exists()}")
PY

section "5) model artifacts"
find artifacts -type f \( -iname "*.joblib" -o -iname "*.json" \) 2>/dev/null | sort | head -n 200

section "6) validate_platform.py"
if [[ -f validate_platform.py ]]; then
  python validate_platform.py || warn "validate_platform.py returned non-zero"
else
  warn "validate_platform.py not found"
fi

section "7) inference smoke test"
python - <<'PY'
import os, sys
from pathlib import Path
sys.path.insert(0, os.path.abspath("src"))

from oligotox.pipelines.score_candidate import score_candidate

processed_dir = Path("data/processed")
model_dir = Path("artifacts/models")
config_dir = Path("artifacts/config")

test_features = {
    "seq_length": 20,
    "gc_content": 0.5,
    "cpg_density": 0.1,
    "ps_fraction": 0.8,
    "dose_mg_kg": 5.0,
}

result = score_candidate(
    candidate_features=test_features,
    processed_dir=processed_dir,
    model_dir=model_dir,
    config_dir=config_dir
)

if isinstance(result, dict):
    print("[OK] top-level keys:", sorted(result.keys()))
    print("[OK] composite_score:", result.get("composite_score"))
    print("[OK] composite_band:", result.get("composite_band"))
    print("[OK] dominant_endpoint:", result.get("dominant_endpoint"))
    pe = result.get("per_endpoint", {})
    print("[OK] endpoint_count:", len(pe) if isinstance(pe, dict) else "n/a")
else:
    print("[WARN] unexpected inference output type:", type(result))
PY

echo
echo "== project health done =="
