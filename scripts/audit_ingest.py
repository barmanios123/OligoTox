import pandas as pd

PATH = "data/raw/aso_tox_ingest_template.csv"

LABEL_COLS = [
    "hepatocellular_injury","hepatic_vacuolation","renal_tubular_injury",
    "renal_accumulation","thrombocytopenia","coagulation_effect",
    "complement_activation","innate_immune_activation","cytokine_release",
    "kupffer_cell_activation","cns_toxicity","offtarget_gene_silencing",
    "injection_site_reaction",
]

KEY_COLS = ["oligo_id","candidate_id","source_paper","source","series_id","sequence","modality","vehicle","route"]

def nonempty(x):
    if pd.isna(x):
        return False
    s = str(x).strip().lower()
    return s not in ("", "nan", "none")

def is_positive_like(x):
    if pd.isna(x):
        return False
    s = str(x).strip().lower()
    return s in ("1", "1.0", "true", "yes")

df = pd.read_csv(PATH)

print("=" * 80)
print("INGEST AUDIT")
print("=" * 80)
print("Rows:", len(df))
print()

print("Missingness in key columns")
for c in KEY_COLS:
    if c in df.columns:
        miss = (~df[c].apply(nonempty)).sum()
        print(f"{c:<18} missing={miss:<5} present={len(df)-miss}")

print()
print("Endpoint positives")
for c in LABEL_COLS:
    if c in df.columns:
        pos = df[c].apply(is_positive_like).sum()
        print(f"{c:<28} positives={pos}")

print()
has_any_label = df[LABEL_COLS].apply(lambda r: any(is_positive_like(v) for v in r), axis=1)
print("Rows with >=1 positive endpoint:", int(has_any_label.sum()))
print("Rows with no positive endpoints:", int((~has_any_label).sum()))

print()
print("Sample rows lacking sequence/modality but having notes:")
mask = (~df.get("sequence", pd.Series(index=df.index, dtype=object)).apply(nonempty)) & (
    df.get("notes", pd.Series(index=df.index, dtype=object)).apply(nonempty)
)
cols = [c for c in ["oligo_id","candidate_id","source","series_id","notes"] if c in df.columns]
print(df.loc[mask, cols].head(10).to_string(index=False))
