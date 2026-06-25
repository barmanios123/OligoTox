from __future__ import annotations
from pathlib import Path
import pandas as pd

SRC = Path("data/raw/cns_neurotox/oligogym/MOE_Neurotox_1_processed.csv.gz")
OUT = Path("data/raw/moe_neurotox_oligogym_curated_append.csv")

df = pd.read_csv(SRC)

def label_cns(row):
    score = row.get("y", 0)
    return int(score >= 3)

seq = df["fasta"].astype(str)

cur = pd.DataFrame(
    {
        "oligo_id": df.index.to_series().map(lambda i: f"OTX-MOE-NTX-{i:04d}"),
        "candidate_id": df.index.to_series().map(lambda i: f"OTX-MOE-NTX-{i:04d}"),
        "source": "MOE_Neurotox_OligoGym",
        "notes": "Imported from OligoGym MOE_Neurotox_1_processed; labels from neurotox score threshold",
        "sequence": seq,
        "modality": "ASO_gapmer",
        "length": seq.str.len(),
        "gc_content": seq.str.count("[GCgc]") / seq.str.len(),
        "cpg_density": seq.str.count("CG") / seq.str.len(),
        "ps_fraction": 1.0,
        "lna_fraction": 0.0,
        "moe_fraction": 0.4,
        "sugar_mod": "MOE",
        "conjugate": "none",
        "delivery_vehicle": "naked",
        "route": "intrathecal",
        "dose_mg_kg": 0.0,
        "repeat_dose": "repeat",
        "target_gene": df["targets"],
        "target_tissue": "CNS",
        "species": "mouse",
        "study_type": "in_vivo_cns_neurotox",
        "source_paper": "MOE_Neurotox",
        "source_table": "OligoGym_processed",
        "platelet_count_baseline": None,
        "platelet_count_treated": None,
        "pf4_ps_elisa_od": None,
        "platelet_c3d_c4d_positive": None,
        "alt_fold": None,
        "ast_fold": None,
        "hepatocellular_injury": 0,
        "hepatic_vacuolation": 0,
        "renal_tubular_injury": 0,
        "renal_accumulation": 0,
        "thrombocytopenia": 0,
        "coagulation_effect": 0,
        "complement_activation": 0,
        "innate_immune_activation": 0,
        "cytokine_release": 0,
        "kupffer_cell_activation": 0,
        "cns_toxicity": df.apply(label_cns, axis=1),
        "offtarget_gene_silencing": 0,
        "injection_site_reaction": 0,
    }
)

OUT.parent.mkdir(parents=True, exist_ok=True)
cur.to_csv(OUT, index=False)
print(f"Wrote {OUT} with {len(cur)} rows")
