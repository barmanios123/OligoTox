from pathlib import Path
import pandas as pd

BASE_PATH = Path("data/processed/training_dataset_rebuilt.csv")
OUT_ENDPOINT = Path("data/processed/audit_sequence_coverage_by_endpoint.csv")
OUT_MODALITY = Path("data/processed/audit_sequence_coverage_by_modality.csv")

LABEL_CANDIDATES = [
    "hepatocellular_injury",
    "hepatic_vacuolation",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "kupffer_cell_activation",
    "cns_toxicity",
    "offtarget_gene_silencing",
    "injection_site_reaction",
]

MODALITY_CANDIDATES = ["modality", "modality_pmo", "modality_pap", "modality_omics"]


def normalize_text(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
        .str.strip()
        .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
    )


def first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


df = pd.read_csv(BASE_PATH, low_memory=False)

if "sequence" not in df.columns:
    raise ValueError("training_dataset.csv must contain a sequence column")

df["sequence"] = normalize_text(df["sequence"])
df["has_sequence"] = df["sequence"].notna()

mod_col = first_present(df, MODALITY_CANDIDATES)
if mod_col is None:
    df["__modality_for_audit"] = "UNKNOWN"
    mod_col = "__modality_for_audit"
else:
    df[mod_col] = normalize_text(df[mod_col]).fillna("UNKNOWN")

label_cols = [c for c in LABEL_CANDIDATES if c in df.columns]
if not label_cols:
    raise ValueError("No expected label columns found in training_dataset.csv")

endpoint_rows = []
for label in label_cols:
    labeled_mask = df[label].notna()
    positive_mask = df[label] == 1
    negative_mask = df[label] == 0

    endpoint_rows.append({
        "endpoint": label,
        "rows_total": int(len(df)),
        "rows_labeled": int(labeled_mask.sum()),
        "rows_with_sequence": int(df["has_sequence"].sum()),
        "rows_labeled_with_sequence": int((labeled_mask & df["has_sequence"]).sum()),
        "rows_positive": int(positive_mask.sum()),
        "rows_positive_with_sequence": int((positive_mask & df["has_sequence"]).sum()),
        "rows_negative": int(negative_mask.sum()),
        "rows_negative_with_sequence": int((negative_mask & df["has_sequence"]).sum()),
        "pct_labeled_with_sequence": round(
            100 * ((labeled_mask & df["has_sequence"]).sum() / labeled_mask.sum()), 2
        ) if labeled_mask.sum() else 0.0,
    })

endpoint_audit = pd.DataFrame(endpoint_rows).sort_values(
    ["pct_labeled_with_sequence", "rows_labeled_with_sequence", "rows_labeled"],
    ascending=[True, True, False]
)

modality_rows = []
for modality, g in df.groupby(mod_col, dropna=False):
    row = {
        "modality": modality,
        "rows_total": int(len(g)),
        "rows_with_sequence": int(g["has_sequence"].sum()),
        "pct_rows_with_sequence": round(100 * g["has_sequence"].mean(), 2) if len(g) else 0.0,
    }
    for label in label_cols:
        labeled_mask = g[label].notna()
        row[f"{label}__labeled"] = int(labeled_mask.sum())
        row[f"{label}__labeled_with_sequence"] = int((labeled_mask & g["has_sequence"]).sum())
    modality_rows.append(row)

modality_audit = pd.DataFrame(modality_rows).sort_values(
    ["pct_rows_with_sequence", "rows_with_sequence", "rows_total"],
    ascending=[True, True, False]
)

endpoint_audit.to_csv(OUT_ENDPOINT, index=False)
modality_audit.to_csv(OUT_MODALITY, index=False)

print("=== SEQUENCE COVERAGE AUDIT ===")
print("rows total:", len(df))
print("rows with sequence:", int(df["has_sequence"].sum()))
print("rows missing sequence:", int((~df["has_sequence"]).sum()))
print()
print("Top endpoints needing sequence curation:")
print(endpoint_audit[[
    "endpoint",
    "rows_labeled",
    "rows_labeled_with_sequence",
    "pct_labeled_with_sequence",
    "rows_positive",
    "rows_positive_with_sequence"
]].head(10).to_string(index=False))
print()
print("Modalities with weakest sequence coverage:")
print(modality_audit[[
    "modality",
    "rows_total",
    "rows_with_sequence",
    "pct_rows_with_sequence"
]].head(10).to_string(index=False))
print()
print("wrote:", OUT_ENDPOINT)
print("wrote:", OUT_MODALITY)