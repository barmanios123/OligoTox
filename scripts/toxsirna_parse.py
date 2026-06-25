from __future__ import annotations
from pathlib import Path
import pandas as pd

RAW = Path("data/raw/sirna_modified/toxsirna_2026_scores.xlsx")
OUT = Path("data/raw/toxsirna_curated_append.csv")

df = pd.read_excel(RAW, sheet_name=0)

# EXPECTED columns: 'siRNA_id','sequence','tox_score','innate_score','offtarget_score', etc.
def label_innate(row):
    return int(row.get("innate_score", 0) >= 0.7)

def label_offtarget(row):
    return int(row.get("offtarget_score", 0) >= 0.7)

cur = pd.DataFrame(
    {
        "oligo_id": df["siRNA_id"],
        "sequence": df["sequence"],
        "modality": "siRNA",
        "length": df["sequence"].str.len(),
        "gc_content": df["sequence"].str.count("[GCgc]") / df["sequence"].str.len(),
        "cpg_density": df["sequence"].str.count("CG") / df["sequence"].str.len(),
        "ps_fraction": 0.0,
        "lna_fraction": 0.0,
        "moe_fraction": 0.0,
        "sugar_mod": "2F+2OMe",
        "conjugate": "none",
        "delivery_vehicle": "naked",
        "route": "in_vitro",
        "dose_mg_kg": 0.0,
        "repeat_dose": "single",
        "target_gene": df.get("target", "unspecified"),
        "target_tissue": "cell_line",
        "species": "human_cell",
        "study_type": "in_vitro_toxsirna",
        "source_paper": "TOXsiRNA_2026",
        "source_table": "supplementary",
        "hepatocellular_injury": 0,
        "hepatic_vacuolation": 0,
        "thrombocytopenia": 0,
        "complement_activation": 0,
        "innate_immune_activation": df.apply(label_innate, axis=1),
        "cytokine_release": 0,
        "renal_tubular_injury": 0,
        "renal_accumulation": 0,
        "offtarget_gene_silencing": df.apply(label_offtarget, axis=1),
        "coagulation_effect": 0,
        "cns_toxicity": 0,
        "injection_site_reaction": 0,
        "kupffer_cell_activation": 0,
        "curation_notes": "TOXsiRNA 2026; labels from tox score thresholds",
        "candidate_id": df["siRNA_id"],
        "source": "TOXsiRNA_2026",
        "notes": "auto-imported from toxsirna_2026_scores.xlsx",
        "vehicle": "naked",
    }
)

OUT.parent.mkdir(parents=True, exist_ok=True)
cur.to_csv(OUT, index=False)
print(f"Wrote {OUT} with {len(cur)} rows")
