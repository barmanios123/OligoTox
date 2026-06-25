import pandas as pd

CURATED = "data/raw/curated_oligo_tox_template_with_approved.csv"
INGEST = "data/raw/aso_tox_ingest_template.csv"
OUT = "data/raw/curated_oligo_tox_template_v2.csv"

LABEL_COLS = [
    "hepatocellular_injury","hepatic_vacuolation","renal_tubular_injury",
    "renal_accumulation","thrombocytopenia","coagulation_effect",
    "complement_activation","innate_immune_activation","cytokine_release",
    "kupffer_cell_activation","cns_toxicity","offtarget_gene_silencing",
    "injection_site_reaction",
]

def is_positive_like(x):
    if pd.isna(x):
        return False
    return str(x).strip().lower() in ("1", "1.0", "true", "yes")

curated = pd.read_csv(CURATED)
ingest = pd.read_csv(INGEST)

if "source_table" not in ingest.columns and "source_table" in curated.columns:
    ingest["source_table"] = ""

if "target_class" not in ingest.columns and "target_class" in curated.columns:
    ingest["target_class"] = ""

keep_mask = ingest[LABEL_COLS].apply(lambda r: any(is_positive_like(v) for v in r), axis=1)
ingest_keep = ingest.loc[keep_mask].copy()

common = [c for c in curated.columns if c in ingest_keep.columns]
ingest_keep = ingest_keep[common]

key_cols = [c for c in ["oligo_id", "candidate_id", "source_paper", "source"] if c in curated.columns and c in ingest_keep.columns]

merged = pd.concat([curated, ingest_keep], ignore_index=True)

if key_cols:
    merged = merged.drop_duplicates(subset=key_cols, keep="first")
else:
    merged = merged.drop_duplicates()

merged.to_csv(OUT, index=False)
print(f"Wrote {OUT}")
print(f"Curated rows: {len(curated)}")
print(f"Ingest positive rows kept: {len(ingest_keep)}")
print(f"Final merged rows: {len(merged)}")
