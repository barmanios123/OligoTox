import pandas as pd
import re

PATH = "data/raw/aso_tox_ingest_template.csv"
OUT = "data/raw/aso_tox_ingest_template.csv"

LABEL_MAP = {
    "hepatocellular_injury": [r"hepat", r"alt", r"ast", r"transaminase", r"liver"],
    "thrombocytopenia": [r"thrombocytopen", r"platelet"],
    "complement_activation": [r"complement", r"c3a", r"c5a", r"sc5b-?9", r"carpa"],
    "coagulation_effect": [r"coagulation", r"coagulopathy", r"pt\b", r"aptt", r"inr"],
    "innate_immune_activation": [r"innate immune", r"tlr", r"immune activation"],
    "cytokine_release": [r"cytokine", r"il-?6", r"tnf", r"ifn", r"interferon"],
    "renal_tubular_injury": [r"renal", r"kidney", r"creatinine", r"bun", r"tubular"],
    "injection_site_reaction": [r"injection site", r"isr", r"local reaction", r"injection reaction"],
    "cns_toxicity": [r"cns", r"neurolog", r"neurotox"],
    "offtarget_gene_silencing": [r"off-?target", r"hybridization", r"rnase h1"],
    "kupffer_cell_activation": [r"kupffer"],
    "hepatic_vacuolation": [r"vacuolat"],
    "renal_accumulation": [r"renal accumulation", r"kidney accumulation"],
}

MODALITY_MAP = {
    "sirna": "siRNA",
    "aso": "ASO",
    "gapmer": "ASO_gapmer",
    "splice": "ASO_splice",
    "pmo": "PMO",
}

VEHICLE_MAP = {
    "lnp": "LNP",
    "galnac": "GalNAc",
    "saline": "saline",
    "none": "none",
}

def nonempty(x):
    if pd.isna(x):
        return False
    return str(x).strip().lower() not in ("", "nan", "none")

df = pd.read_csv(PATH)

text_cols = [c for c in ["key_tox_signals_from_label", "notes", "source_paper", "source"] if c in df.columns]

for idx, row in df.iterrows():
    blob = " | ".join(str(row[c]) for c in text_cols if c in row and pd.notna(row[c])).lower()

    for label, patterns in LABEL_MAP.items():
        if label in df.columns and not nonempty(row.get(label, "")):
            if any(re.search(p, blob) for p in patterns):
                df.at[idx, label] = 1

    chem = str(row.get("notes", "")).lower() + " | " + str(row.get("source_paper", "")).lower()
    if "modality" in df.columns and not nonempty(row.get("modality", "")):
        for k, v in MODALITY_MAP.items():
            if k in chem:
                df.at[idx, "modality"] = v
                break

    veh_blob = blob + " | " + chem
    if "vehicle" in df.columns and not nonempty(row.get("vehicle", "")):
        for k, v in VEHICLE_MAP.items():
            if k in veh_blob:
                df.at[idx, "vehicle"] = v
                break

    if "route" in df.columns and not nonempty(row.get("route", "")):
        if re.search(r"\biv\b|intraven", blob):
            df.at[idx, "route"] = "iv"
        elif re.search(r"\bsc\b|subcut", blob):
            df.at[idx, "route"] = "sc"
        elif re.search(r"\bit\b|intrathecal", blob):
            df.at[idx, "route"] = "it"

df.to_csv(OUT, index=False)
print(f"Wrote {OUT} with auto-labeled approved-drug rows")
