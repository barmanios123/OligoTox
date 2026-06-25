from pathlib import Path
import pandas as pd

BASE_PATH = Path("data/processed/training_dataset.csv")
OUT_PATH = Path("data/metadata/candidate_to_oligo_map.csv")

CANDIDATE_HINT_COLS = [
    "candidate_id",
    "drug_name",
    "brand",
    "compound",
    "gapmer_id",
    "target_gene",
    "modality",
    "chemistry",
    "delivery",
    "route",
    "primary_source",
    "source_doi",
]

def normalize_text(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
        .str.strip()
        .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
    )

df = pd.read_csv(BASE_PATH, low_memory=False)

if "sequence" not in df.columns:
    raise ValueError("training_dataset.csv must contain a sequence column")

df["sequence"] = normalize_text(df["sequence"])
missing = df[df["sequence"].isna()].copy()

keep_cols = [c for c in CANDIDATE_HINT_COLS if c in missing.columns]
template = missing[keep_cols].copy()

template["oligo_id"] = pd.NA
template["registry_id"] = pd.NA
template["manual_sequence"] = pd.NA
template["manual_modification"] = pd.NA
template["mapping_notes"] = pd.NA
template["curation_status"] = "needs_mapping"

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
template.to_csv(OUT_PATH, index=False)

print("rows needing mapping:", len(template))
print("wrote:", OUT_PATH)
print("columns:", list(template.columns))
