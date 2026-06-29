from pathlib import Path
import pandas as pd

KC_DEFAULT = Path("data/processed/kupffer_activation_30plus_examples.csv")
TRAIN_DEFAULT = Path("data/processed/training_dataset.csv")
OUT_DEFAULT = Path("data/processed/training_dataset_v3.csv")

def normalize_seq(s):
    return str(s).strip().upper().replace(" ", "").replace("U", "T")

def merge_kc_examples(kc_path=KC_DEFAULT, training_path=TRAIN_DEFAULT, out_path=OUT_DEFAULT):
    print(f"Loading KC examples from {kc_path}")
    kc = pd.read_csv(kc_path)
    train = pd.read_csv(training_path)

    kc["sequence_norm"] = kc["sequence"].map(normalize_seq)
    train["sequence_norm"] = train["sequence"].map(normalize_seq)

    kc = kc.rename(columns={
        "label_hepatic_vacuolation": "endpoint_hepatic_vacuolation",
        "label_kupffer_activation": "endpoint_kupffer_cell_activation",
    })

    needed = ["sequence_norm", "endpoint_hepatic_vacuolation", "endpoint_kupffer_cell_activation"]
    missing = [c for c in needed if c not in kc.columns]
    if missing:
        raise ValueError(f"KC CSV missing columns: {missing}")

    kc_small = kc[needed].drop_duplicates(subset=["sequence_norm"]).copy()
    merged = train.merge(kc_small, on="sequence_norm", how="left", suffixes=("", "_kc"))

    for col in ["endpoint_hepatic_vacuolation", "endpoint_kupffer_cell_activation"]:
        kc_col = f"{col}_kc"
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0)
        merged[kc_col] = pd.to_numeric(merged[kc_col], errors="coerce").fillna(0)
        merged[col] = merged[[col, kc_col]].max(axis=1)
        merged.drop(columns=[kc_col], inplace=True)

    merged.drop(columns=["sequence_norm"], inplace=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)

    print(f"Saved merged dataset to {out_path}")
    for col in ["endpoint_hepatic_vacuolation", "endpoint_kupffer_cell_activation"]:
        s = pd.to_numeric(merged[col], errors="coerce").fillna(0)
        print(f"{col}: positives = {int((s > 0).sum())}")

if __name__ == "__main__":
    merge_kc_examples()
