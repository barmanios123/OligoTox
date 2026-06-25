from pathlib import Path
import pandas as pd

MAP_PATH = Path("data/metadata/candidate_to_oligo_map.csv")

def normalize_text(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
        .str.strip()
        .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
    )

df = pd.read_csv(MAP_PATH, low_memory=False)

for col in [
    "candidate_id",
    "oligo_id",
    "registry_id",
    "manual_sequence",
    "manual_modification",
    "curation_status",
]:
    if col in df.columns:
        df[col] = normalize_text(df[col])

empty = pd.Series([pd.NA] * len(df), index=df.index)

has_any_mapping = (
    df.get("oligo_id", empty).notna()
    | df.get("registry_id", empty).notna()
    | df.get("manual_sequence", empty).notna()
)

print("rows total:", len(df))
print("rows with any mapping info:", int(has_any_mapping.sum()))
print("rows still unresolved:", int((~has_any_mapping).sum()))

if "manual_sequence" in df.columns:
    seq_len = df["manual_sequence"].dropna().astype(str).str.len()
    if not seq_len.empty:
        print("manual_sequence length min/max:", int(seq_len.min()), int(seq_len.max()))

dupes = df["candidate_id"].duplicated().sum() if "candidate_id" in df.columns else 0
print("duplicate candidate_id rows:", int(dupes))
