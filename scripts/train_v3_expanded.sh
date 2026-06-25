#!/usr/bin/env bash
set -euo pipefail
cd ~/oligotox
source .venv/bin/activate
cp scripts/train_v3.py scripts/train_v3.py.expandedbak
python - <<'PY'
from pathlib import Path
p = Path("scripts/train_v3.py")
t = p.read_text()
old = '''def load_data():
    candidate_paths = [
        Path("data/curated/curated_oligo_tox_master_v3.csv"),
        Path("data/curated/curated_oligo_tox_template_v3.csv"),
        Path("downloads/curated_oligo_tox_master_v3.csv"),
        Path("curated_oligo_tox_master_v3.csv"),
    ]
    for csv_path in candidate_paths:
        if csv_path.exists():
            return pd.read_csv(csv_path)
    raise FileNotFoundError("No curated CSV found.")'''
new = '''def load_data():
    candidate_paths = [
        Path("data/curated/curated_oligo_tox_master_v3_expanded.csv"),
        Path("data/curated/curated_oligo_tox_master_v3.csv"),
        Path("data/curated/curated_oligo_tox_template_v3.csv"),
        Path("downloads/curated_oligo_tox_master_v3.csv"),
        Path("curated_oligo_tox_master_v3.csv"),
    ]
    for csv_path in candidate_paths:
        if csv_path.exists():
            return pd.read_csv(csv_path)
    raise FileNotFoundError("No curated CSV found.")'''
if old in t:
    t = t.replace(old, new)
else:
    if 'Path("data/curated/curated_oligo_tox_master_v3_expanded.csv")' not in t:
        raise SystemExit("Could not patch load_data() safely; inspect scripts/train_v3.py manually.")
p.write_text(t)
print("patched scripts/train_v3.py for expanded training")
PY
bash scripts/train_from_v03.sh
