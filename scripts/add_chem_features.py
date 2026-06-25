import pandas as pd
from oligotox.features.aso_optimizer_v03 import batch_aso_features

PATH_IN = "data/processed/training_dataset_rebuilt.csv"
PATH_OUT = "data/processed/training_dataset_with_chem.csv"

def main() -> None:
    df = pd.read_csv(PATH_IN)

    if "sequence" not in df.columns:
        raise SystemExit("ERROR: training_dataset_rebuilt.csv does not have a 'sequence' column")

    seq_raw = df["sequence"]
    valid_mask = (
        seq_raw.notna()
        & seq_raw.astype(str).str.strip().ne("")
        & seq_raw.astype(str).str.lower().ne("nan")
    )
    dropped = int((~valid_mask).sum())
    if dropped:
        print("Dropping rows with missing/blank sequence:", dropped)

    df = df.loc[valid_mask].copy().reset_index(drop=True)
    sequences = df["sequence"].astype(str).str.strip().tolist()

    # Safely build mods list: use None for missing/nan modification strings
    if "modification" in df.columns:
        def safe_mod(v):
            s = str(v).strip()
            return None if s.lower() in ("nan", "", "none", "<na>") else s
        mods = [safe_mod(v) for v in df["modification"]]
    else:
        mods = None

    chem_list = batch_aso_features(sequences, mods=mods, ps_fractions=None)
    chem_df = pd.DataFrame(chem_list)

    existing_chem = [c for c in df.columns if c.startswith("chem_")]
    if existing_chem:
        df = df.drop(columns=existing_chem)

    df_chem = pd.concat([df.reset_index(drop=True), chem_df.reset_index(drop=True)], axis=1)
    df_chem.to_csv(PATH_OUT, index=False)

    print("chem_df shape:", chem_df.shape)
    print("Augmented df shape:", df_chem.shape)
    print("Wrote:", PATH_OUT)

if __name__ == "__main__":
    main()
