from __future__ import annotations
from pathlib import Path
import pandas as pd

RAW_CURATED = Path("data/raw/curated_oligo_tox_template_v3.csv")
OUT = Path("data/raw/cns_toxicity_curation_template.csv")

df = pd.read_csv(RAW_CURATED)

cols_keep = [
    "oligo_id","candidate_id","source","notes",
    "modality","delivery_vehicle","route","species","dose_mg_kg","repeat_dose",
    "cns_toxicity",
]

missing = [c for c in cols_keep if c not in df.columns]
for c in missing:
    df[c] = None

tpl = df[cols_keep].copy()

mask = tpl["source"].fillna("").str.contains(
    "nusinersen|scaso|nhp_deepbrain_aso|it|intrathecal",
    case=False,
    regex=True,
)
tpl = tpl[mask].copy()

# skeleton rows to fill from scASO / NHP deep brain ASO / nusinersen NDA
for tag in ["scASO", "nhp_deepbrain_aso", "nusinersen"]:
    tpl = pd.concat(
        [
            tpl,
            pd.DataFrame(
                {
                    "oligo_id": [f"OTX-CNS-{tag.upper()}-001"],
                    "candidate_id": [f"OTX-CNS-{tag.upper()}-001"],
                    "source": [f"{tag}_cns_data"],
                    "notes": [f"{tag} CNS tox row — fill NfL/GFAP/histopath evidence"],
                }
            ),
        ],
        ignore_index=True,
    )

tpl.to_csv(OUT, index=False)
print(f"Wrote {OUT} with {len(tpl)} rows")
