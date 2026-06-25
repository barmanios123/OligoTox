from __future__ import annotations
from pathlib import Path
import pandas as pd

RAW_CURATED = Path("data/raw/curated_oligo_tox_template_v3.csv")
OUT = Path("data/raw/thrombocytopenia_curation_template.csv")

df = pd.read_csv(RAW_CURATED)

# keep only columns relevant for platelet endpoints + IDs/notes
cols_keep = [
    "oligo_id","candidate_id","source","notes",
    "modality","delivery_vehicle","route","species","dose_mg_kg","repeat_dose",
    "platelet_count_baseline","platelet_count_treated",
    "pf4_ps_elisa_od","platelet_c3d_c4d_positive",
    "thrombocytopenia",
]

missing = [c for c in cols_keep if c not in df.columns]
for c in missing:
    df[c] = None

tpl = df[cols_keep].copy()

# filter to likely platelet-relevant sources
mask = tpl["source"].fillna("").str.contains(
    "inotersen|eplontersen|mipomersen|sider|uniox|unitox|nda_labels|sider_chembl",
    case=False,
    regex=True,
)
tpl = tpl[mask].copy()

# add a couple of empty rows for new NDA-based entries
for drug in ["inotersen", "eplontersen", "mipomersen", "vutrisiran", "givosiran"]:
    tpl = pd.concat(
        [
            tpl,
            pd.DataFrame(
                {
                    "oligo_id": [f"OTX-PLT-{drug.upper()}-001"],
                    "candidate_id": [f"OTX-PLT-{drug.upper()}-001"],
                    "source": [f"{drug}_NDA_label"],
                    "notes": [f"{drug} NDA platelet signal row — fill counts and PF4/C3d"],
                }
            ),
        ],
        ignore_index=True,
    )

tpl.to_csv(OUT, index=False)
print(f"Wrote {OUT} with {len(tpl)} rows")
