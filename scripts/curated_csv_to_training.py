from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "data/processed/curated_oligo_tox_master_v2.csv"
LEGACY = ROOT / "data/curated/curated_oligo_tox_master.csv"
OUT = ROOT / "data/processed/training_matrix.csv"

source = V2 if V2.exists() else LEGACY
df = pd.read_csv(source)

endpoint_cols = [
    "hepatocellular_injury",
    "renal_tubular_injury",
    "hypersensitivity",
    "injection_site_reaction",
    "innate_immune_activation",
    "offtarget_gene_silencing",
    "thrombocytopenia",
    "complement_activation",
    "cns_neurotoxicity",
    "reproductive_dev_tox",
]

for col in endpoint_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

df.to_csv(OUT, index=False)
print(f"Using source: {source}")
print(f"Wrote: {OUT}")
print(df.shape)