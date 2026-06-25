from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
NEW_BATCH = ROOT / "data/processed/pmo_training_ready.csv"
MASTER = ROOT / "data/curated/curated_oligo_tox_master.csv"

def normalize_key(s):
    return (
        s.astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"\s+", "_", regex=True)
    )

def ensure_oligo_id(df):
    df = df.copy()
    if "oligo_id" in df.columns:
        df["oligo_id"] = normalize_key(df["oligo_id"])
        return df
    if "candidate_id" in df.columns:
        df["oligo_id"] = normalize_key(df["candidate_id"])
        return df
    if "drug_name" in df.columns:
        df["oligo_id"] = normalize_key(df["drug_name"])
        return df
    raise ValueError("Need one of: oligo_id, candidate_id, drug_name")

def main():
    if not NEW_BATCH.exists():
        raise FileNotFoundError(f"Missing batch file: {NEW_BATCH}")

    new_df = pd.read_csv(NEW_BATCH)
    new_df = ensure_oligo_id(new_df)

    if MASTER.exists():
        master_df = pd.read_csv(MASTER)
        master_df = ensure_oligo_id(master_df)
    else:
        master_df = pd.DataFrame(columns=new_df.columns)

    all_cols = sorted(set(master_df.columns).union(new_df.columns))
    master_df = master_df.reindex(columns=all_cols)
    new_df = new_df.reindex(columns=all_cols)

    if "source_type" in new_df.columns:
        new_df["__source_priority"] = (new_df["source_type"].fillna("").astype(str).str.lower() == "label").astype(int)
    else:
        new_df["__source_priority"] = 1

    if "__source_priority" not in master_df.columns:
        master_df["__source_priority"] = 0

    combined = pd.concat([master_df, new_df], ignore_index=True)
    combined = combined.sort_values(["oligo_id", "__source_priority"], ascending=[True, False])
    combined = combined.drop_duplicates(subset=["oligo_id"], keep="first")
    combined = combined.drop(columns=["__source_priority"], errors="ignore")
    combined = combined.sort_values("oligo_id").reset_index(drop=True)

    MASTER.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(MASTER, index=False)

    print(f"Wrote master file: {MASTER}")
    print(f"Rows: {len(combined)}")
    print(f"Columns: {len(combined.columns)}")

if __name__ == '__main__':
    main()
