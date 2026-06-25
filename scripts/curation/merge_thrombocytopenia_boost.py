from __future__ import annotations
from pathlib import Path
import pandas as pd

BASE = Path("data/raw/curated_oligo_tox_template_v3.csv")
BOOST = Path("data/raw/thrombocytopenia_boost_template.csv")
OUT = Path("data/raw/curated_oligo_tox_template_v4.csv")

base = pd.read_csv(BASE)
boost = pd.read_csv(BOOST)

merged = pd.concat([base, boost], ignore_index=True)
merged.to_csv(OUT, index=False)

print(f"Wrote {OUT} with {len(merged)} rows")
print("Added sources:")
print(boost['source'].value_counts().to_string())
