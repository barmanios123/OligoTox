from pathlib import Path
import pandas as pd
import sys

COMMON_COLUMNS = [
    "sequence",
    "modality",
    "chemistry",
    "conjugate",
    "endpoint",
    "label",
    "label_type",
    "source_dataset",
    "target_gene",
    "sequence_family",
    "notes",
]

def normalize_dataset(input_csv: Path, output_csv: Path, dataset_name: str, endpoint: str, label_type: str, modality: str):
    df = pd.read_csv(input_csv)
    lower = {c.lower(): c for c in df.columns}

    seq_col = next((lower[c] for c in ["sequence", "oligo_sequence", "sirna", "aso_sequence"] if c in lower), None)
    label_col = next((lower[c] for c in ["label", "tox", "toxicity", "response", "activity", "score"] if c in lower), None)

    if seq_col is None or label_col is None:
        raise ValueError(f"{dataset_name}: could not infer sequence/label columns from {list(df.columns)}")

    out = pd.DataFrame()
    out["sequence"] = df[seq_col].astype(str).str.upper().str.replace(" ", "", regex=False)
    out["modality"] = modality
    out["chemistry"] = df[lower["chemistry"]] if "chemistry" in lower else "unknown"
    out["conjugate"] = df[lower["conjugate"]] if "conjugate" in lower else "none"
    out["endpoint"] = endpoint
    out["label"] = df[label_col]
    out["label_type"] = label_type
    out["source_dataset"] = dataset_name
    out["target_gene"] = df[lower["target_gene"]] if "target_gene" in lower else ""
    out["sequence_family"] = df[lower["sequence_family"]] if "sequence_family" in lower else out["sequence"].str[:8]
    out["notes"] = ""

    for c in COMMON_COLUMNS:
        if c not in out.columns:
            out[c] = ""

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out[COMMON_COLUMNS].to_csv(output_csv, index=False)
    print(f"wrote {output_csv} rows={len(out)}")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        raise SystemExit("usage: python normalize_benchmark_dataset.py <input_csv> <output_csv> <dataset_name> <endpoint> <label_type> <modality>")
    normalize_dataset(
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        sys.argv[3],
        sys.argv[4],
        sys.argv[5],
        sys.argv[6],
    )
