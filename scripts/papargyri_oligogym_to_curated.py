from __future__ import annotations
from pathlib import Path
import pandas as pd

SRC = Path("data/raw/lna_gapmer/oligogym/Papargyri_2020_1_processed.csv.gz")
OUT = Path("data/raw/papargyri_oligogym_curated_append.csv")

df = pd.read_csv(SRC)

def label_offtarget(row):
    casp = row.get("y_raw", 1.0)
    return int(casp >= 300.0)

def label_hep(row):
    casp = row.get("y_raw", 1.0)
    return int(casp >= 300.0)

seq = df["fasta"].astype(str)

cur = pd.DataFrame(
    {
        "oligo_id": df.index.to_series().map(lambda i: f"OTX-PAP-OG-{i:04d}"),
        "candidate_id": df.index.to_series().map(lambda i: f"OTX-PAP-OG-{i:04d}"),
        "source": "Papargyri_2020_OligoGym",
        "notes": "Imported from OligoGym Papargyri_2020_1_processed; labels from y_raw threshold",
        "sequence": seq,
        "modality": "ASO_gapmer",
        "length": seq.str.len(),
        "gc_content": seq.str.count("[GCgc]") / seq.str.len(),
        "cpg_density": seq.str.count("CG") / seq.str.len(),
        "ps_fraction": 1.0,
        "lna_fraction": 0.3,
        "moe_fraction": 0.0,
        "sugar_mod": "LNA",
        "conjugate": "none",
        "delivery_vehicle": "naked",
        "route": "in_vitro",
        "dose_mg_kg": 0.0,
        "repeat_dose": "single",
        "target_gene": df["targets"],
        "target_tissue": "liver",
        "species": "human_cell",
        "study_type": "in_vitro_caspase",
        "source_paper": "Papargyri_2020",
        "source_table": "OligoGym_processed",
        "platelet_count_baseline": None,
        "platelet_count_treated": None,
        "pf4_ps_elisa_od": None,
        "platelet_c3d_c4d_positive": None,
        "alt_fold": None,
        "ast_fold": None,
        "hepatocellular_injury": df.apply(label_hep, axis=1),
        "hepatic_vacuolation": 0,
        "renal_tubular_injury": 0,
        "renal_accumulation": 0,
        "thrombocytopenia": 0,
        "coagulation_effect": 0,
        "complement_activation": 0,
        "innate_immune_activation": 0,
        "cytokine_release": 0,
        "kupffer_cell_activation": 0,
        "cns_toxicity": 0,
        "offtarget_gene_silencing": df.apply(label_offtarget, axis=1),
        "injection_site_reaction": 0,
    }
)

OUT.parent.mkdir(parents=True, exist_ok=True)
cur.to_csv(OUT, index=False)
print(f"Wrote {OUT} with {len(cur)} rows")
