from pathlib import Path
import json
import pandas as pd

LABEL_COLUMNS = [
    "hepatocellular_injury",
    "hepatic_vacuolation",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "kupffer_cell_activation",
    "cns_toxicity",
    "offtarget_gene_silencing",
    "injection_site_reaction",
]

ID_COLUMNS = ["candidate_id", "source", "notes", "series_id"]


def build_feature_and_label_matrices(training_csv, feature_out, label_out, meta_out):
    training_csv = Path(training_csv)
    feature_out = Path(feature_out)
    label_out = Path(label_out)
    meta_out = Path(meta_out)

    df = pd.read_csv(training_csv)

    label_cols = [c for c in LABEL_COLUMNS if c in df.columns]
    if not label_cols:
        raise ValueError("No label columns found in training dataset.")

    feature_cols = [c for c in df.columns if c not in set(ID_COLUMNS) | set(label_cols)]
    if not feature_cols:
        raise ValueError("No feature columns found after excluding ids and labels.")

    X = df[feature_cols].copy()
    y = df[label_cols].copy()

    feature_out.parent.mkdir(parents=True, exist_ok=True)
    label_out.parent.mkdir(parents=True, exist_ok=True)
    meta_out.parent.mkdir(parents=True, exist_ok=True)

    X.to_csv(feature_out, index=False)
    y.to_csv(label_out, index=False)

    metadata = {
        "training_csv": str(training_csv),
        "rows": int(df.shape[0]),
        "n_features": int(X.shape[1]),
        "n_labels": int(y.shape[1]),
        "feature_columns": feature_cols,
        "label_columns": label_cols,
        "id_columns": [c for c in ID_COLUMNS if c in df.columns],
    }

    meta_out.write_text(json.dumps(metadata, indent=2))
    return metadata