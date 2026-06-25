from pathlib import Path
import pandas as pd

BASE_PATH = Path("data/processed/training_dataset.csv")
OUT_PATH = Path("data/processed/training_dataset_rebuilt.csv")
AUDIT_PATH = Path("data/processed/training_dataset_rebuilt_audit.csv")

CANDIDATE_SOURCES = [
    Path("data/processed/pmo_training_ready.csv"),
    Path("data/processed/pmo_label_training_rows.csv"),
    Path("data/curated/curated_oligo_tox_master_v2.csv"),
    Path("data/metadata/candidate_to_oligo_map.csv"),
]

JOIN_KEYS = ["candidate_id", "oligo_id", "gapmer_id", "compound", "registry_id"]
SEQ_CANDIDATES = ["sequence", "seq", "oligo_sequence", "manual_sequence"]
MOD_CANDIDATES = ["modification", "chemistry", "chemistry_class", "manual_modification"]


def normalize_text(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
        .str.strip()
        .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
    )


def first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def prepare_source_table(src_path: Path) -> pd.DataFrame | None:
    if not src_path.exists():
        return None

    src = pd.read_csv(src_path, low_memory=False)
    seq_col = first_present(src, SEQ_CANDIDATES)
    mod_col = first_present(src, MOD_CANDIDATES)

    if seq_col is None:
        return None

    src = src.copy()
    src[seq_col] = normalize_text(src[seq_col])

    keep_cols = []
    for key in JOIN_KEYS:
        if key in src.columns:
            src[key] = normalize_text(src[key])
            keep_cols.append(key)

    if not keep_cols:
        return None

    cols = keep_cols + [seq_col]
    if mod_col is not None:
        src[mod_col] = normalize_text(src[mod_col])
        cols.append(mod_col)

    slim = src[cols].copy()
    slim = slim.rename(columns={seq_col: "__src_sequence"})
    if mod_col is not None:
        slim = slim.rename(columns={mod_col: "__src_modification"})
    slim["__source_file"] = str(src_path)
    return slim


def build_lookup(sources: list[pd.DataFrame], key: str) -> pd.DataFrame | None:
    frames = []
    for src in sources:
        if src is None or key not in src.columns:
            continue
        cols = [key, "__src_sequence", "__source_file"]
        if "__src_modification" in src.columns:
            cols.append("__src_modification")
        tmp = src[cols].dropna(subset=[key]).drop_duplicates(subset=[key], keep="first").copy()
        tmp = tmp.rename(columns={key: "__join_key"})
        frames.append(tmp)

    if not frames:
        return None

    out = pd.concat(frames, ignore_index=True)
    out = out.drop_duplicates(subset=["__join_key"], keep="first")
    return out


base = pd.read_csv(BASE_PATH, low_memory=False).copy()
base["__row_id"] = range(len(base))

if "sequence" in base.columns:
    base["sequence"] = normalize_text(base["sequence"])
else:
    base["sequence"] = pd.NA

if "modification" in base.columns:
    base["modification"] = normalize_text(base["modification"])

for key in JOIN_KEYS:
    if key in base.columns:
        base[key] = normalize_text(base[key])

before_valid = int(base["sequence"].notna().sum())

source_tables = [prepare_source_table(p) for p in CANDIDATE_SOURCES]
source_tables = [s for s in source_tables if s is not None]
if not source_tables:
    raise ValueError("No usable source tables with sequence were found")

lookups = {key: build_lookup(source_tables, key) for key in JOIN_KEYS}

merged = base.copy()
merged["sequence_before"] = merged["sequence"]
merged["sequence_source_key"] = pd.NA
merged["__source_file"] = pd.NA

if "modification" not in merged.columns:
    merged["modification"] = pd.NA

for key in JOIN_KEYS:
    lookup = lookups.get(key)
    if lookup is None or key not in merged.columns:
        continue

    step = merged[[key]].merge(lookup, left_on=key, right_on="__join_key", how="left")

    need_seq = merged["sequence"].isna()
    got_seq = step["__src_sequence"].notna()
    fill_mask = need_seq & got_seq

    merged.loc[fill_mask, "sequence"] = step.loc[fill_mask, "__src_sequence"].values
    merged.loc[fill_mask, "sequence_source_key"] = key
    merged.loc[fill_mask, "__source_file"] = step.loc[fill_mask, "__source_file"].values

    if "__src_modification" in step.columns:
        need_mod = merged["modification"].isna()
        got_mod = step["__src_modification"].notna()
        mod_mask = need_mod & got_mod
        merged.loc[mod_mask, "modification"] = step.loc[mod_mask, "__src_modification"].values

after_valid = int(merged["sequence"].notna().sum())

audit_cols = ["__row_id", "sequence_before", "sequence", "sequence_source_key", "__source_file"]
for key in JOIN_KEYS:
    if key in merged.columns:
        audit_cols.append(key)

audit = merged[audit_cols].copy()
audit["recovered_sequence"] = audit["sequence_before"].isna() & audit["sequence"].notna()

merged = merged.sort_values("__row_id").drop(columns=["__row_id"])
audit = audit.sort_values("__row_id")

merged.to_csv(OUT_PATH, index=False)
audit.to_csv(AUDIT_PATH, index=False)

print("rows:", len(merged))
print("valid sequence before:", before_valid)
print("valid sequence after:", after_valid)
print("recovered:", after_valid - before_valid)
print("wrote:", OUT_PATH)
print("audit:", AUDIT_PATH)
print("recovered by key:")
print(audit.loc[audit["recovered_sequence"], "sequence_source_key"].value_counts(dropna=False))