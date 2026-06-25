from __future__ import annotations
from pathlib import Path
import pandas as pd

BASE = Path("data/raw/curated_oligo_tox_template_v3.csv")
PAP = Path("data/raw/papargyri_oligogym_curated_append.csv")
MOE = Path("data/raw/moe_neurotox_oligogym_curated_append.csv")
OUT = Path("data/raw/curated_oligo_tox_template_v3_plus_oligogym.csv")

dfs = [pd.read_csv(BASE)]

for p in [PAP, MOE]:
    if p.exists():
        dfs.append(pd.read_csv(p))
    else:
        print(f"WARNING: {p} does not exist, skipping")

merged = pd.concat(dfs, ignore_index=True)

merged.to_csv(OUT, index=False)
print(f"Wrote {OUT} with {len(merged)} rows")
