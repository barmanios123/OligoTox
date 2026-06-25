#!/usr/bin/env bash
set -euo pipefail

cd ~/oligotox

mkdir -p data/curated
mkdir -p reports

MASTER="data/curated/curated_oligo_tox_master_v3.csv"
PAPARGYRI_XLSX="papargyri_2020_LNA_512gapmers_TableS1.xlsx"
PAPARGYRI_OUT="data/curated/papargyri_512gapmers_labeled.csv"
PLATELET_OUT="data/curated/fda_platelet_rows_v3.csv"
MERGED_OUT="data/curated/curated_oligo_tox_master_v3_expanded.csv"

echo "== Pending source integration =="

if [ -f "scripts/papargyri_extractor.py" ] && [ -f "$PAPARGYRI_XLSX" ]; then
  echo "Running Papargyri extractor..."
  python scripts/papargyri_extractor.py --input "$PAPARGYRI_XLSX" --output "$PAPARGYRI_OUT"
else
  echo "Papargyri step skipped: missing extractor or Table S1 workbook."
fi

if [ -f "scripts/fda_nda_platelet_extractor.py" ]; then
  echo "FDA platelet extractor present."
  echo "Skipping execution unless populated rows are intentionally generated."
else
  echo "FDA platelet extractor missing."
fi

python - <<'PY'
from pathlib import Path
import pandas as pd

master = Path("data/curated/curated_oligo_tox_master_v3.csv")
pap = Path("data/curated/papargyri_512gapmers_labeled.csv")
plt = Path("data/curated/fda_platelet_rows_v3.csv")
out = Path("data/curated/curated_oligo_tox_master_v3_expanded.csv")
report = Path("reports/pending_source_merge_report.txt")

if not master.exists():
    raise FileNotFoundError(f"Missing master dataset: {master}")

frames = []
notes = []

df_master = pd.read_csv(master)
frames.append(df_master)
notes.append(f"master_rows={len(df_master)}")

if pap.exists() and pap.stat().st_size > 0:
    df_pap = pd.read_csv(pap)
    if len(df_pap) > 0:
        frames.append(df_pap)
        notes.append(f"papargyri_rows={len(df_pap)}")
    else:
        notes.append("papargyri_rows=0 (empty dataframe)")
else:
    notes.append("papargyri_rows=0 (file not present)")

if plt.exists() and plt.stat().st_size > 0:
    try:
        df_plt = pd.read_csv(plt)
        if len(df_plt) > 0:
            frames.append(df_plt)
            notes.append(f"fda_platelet_rows={len(df_plt)}")
        else:
            notes.append("fda_platelet_rows=0 (empty dataframe)")
    except Exception as e:
        notes.append(f"fda_platelet_rows=0 (unreadable: {e})")
else:
    notes.append("fda_platelet_rows=0 (file not present or empty)")

all_cols = sorted(set().union(*[set(df.columns) for df in frames]))
aligned = [df.reindex(columns=all_cols) for df in frames]
merged = pd.concat(aligned, ignore_index=True)

dedupe_keys = [c for c in ["candidate_id", "compound", "name", "sequence"] if c in merged.columns]
before = len(merged)
if dedupe_keys:
    merged = merged.drop_duplicates(subset=dedupe_keys, keep="last")
after = len(merged)

merged.to_csv(out, index=False)

notes.append(f"merged_rows_before_dedup={before}")
notes.append(f"merged_rows_after_dedup={after}")
notes.append(f"columns={len(merged.columns)}")
notes.append(f"dedupe_keys={dedupe_keys}")

endpoint_cols = [
    "hepatocellular_injury",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "hypersensitivity",
    "cns_neurotoxicity",
    "offtarget_gene_silencing",
    "reproductive_dev_tox",
    "injection_site_reaction",
]
for col in endpoint_cols:
    if col in merged.columns:
        notes.append(f"{col}_positives={int(pd.to_numeric(merged[col], errors='coerce').fillna(0).sum())}")

report.write_text("\n".join(notes) + "\n")
print(f"wrote {out}")
print(f"wrote {report}")
PY

echo "== Done =="
echo "Expanded dataset: $MERGED_OUT"
echo "Merge report: reports/pending_source_merge_report.txt"
