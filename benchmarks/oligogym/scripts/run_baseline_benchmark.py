from pathlib import Path
import argparse
import json
import pandas as pd

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    f1_score,
    balanced_accuracy_score,
    confusion_matrix,
)

def load_data(path):
    df = pd.read_csv(path)
    required = {"sequence", "label", "split"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df

def make_kmer_features(train_seqs, test_seqs, k=3):
    vec = CountVectorizer(analyzer="char", ngram_range=(k, k))
    X_train = vec.fit_transform(train_seqs)
    X_test = vec.transform(test_seqs)
    return X_train, X_test, vec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Split benchmark CSV")
    ap.add_argument("--output-prefix", required=True, help="Prefix for result files")
    args = ap.parse_args()

    df = load_data(args.input)

    train = df[df["split"] == "train"].copy()
    val = df[df["split"] == "val"].copy()
    test = df[df["split"] == "test"].copy()

    if len(train) == 0 or len(test) == 0:
        raise ValueError("Need non-empty train and test splits")

    X_train, X_val, _ = make_kmer_features(train["sequence"], val["sequence"] if len(val) else train["sequence"])
    _, X_test, vec = make_kmer_features(train["sequence"], test["sequence"])

    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    model.fit(vec.transform(train["sequence"]), train["label"])

    eval_frames = []
    for split_name, part in [("val", val), ("test", test)]:
        if len(part) == 0:
            continue
        X = vec.transform(part["sequence"])
        y_true = part["label"].astype(float)
        y_prob = model.predict_proba(X)[:, 1]
        y_pred = (y_prob >= 0.5).astype(int)

        metrics = {
            "split": split_name,
            "n": int(len(part)),
            "roc_auc": float(roc_auc_score(y_true, y_prob)),
            "pr_auc": float(average_precision_score(y_true, y_prob)),
            "f1": float(f1_score(y_true, y_pred)),
            "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
            "tn": int(confusion_matrix(y_true, y_pred).ravel()[0]),
            "fp": int(confusion_matrix(y_true, y_pred).ravel()[1]),
            "fn": int(confusion_matrix(y_true, y_pred).ravel()[2]),
            "tp": int(confusion_matrix(y_true, y_pred).ravel()[3]),
        }
        eval_frames.append(metrics)

    metrics_df = pd.DataFrame(eval_frames)
    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    metrics_csv = out_prefix.with_suffix(".metrics.csv")
    metrics_json = out_prefix.with_suffix(".metrics.json")
    preds_csv = out_prefix.with_suffix(".predictions.csv")

    metrics_df.to_csv(metrics_csv, index=False)

    with open(metrics_json, "w") as f:
        json.dump(metrics_df.to_dict(orient="records"), f, indent=2)

    test_probs = model.predict_proba(vec.transform(test["sequence"]))[:, 1]
    test_preds = (test_probs >= 0.5).astype(int)
    pred_df = test.copy()
    pred_df["pred_prob"] = test_probs
    pred_df["pred_label"] = test_preds
    pred_df.to_csv(preds_csv, index=False)

    print("Wrote:")
    print(" -", metrics_csv)
    print(" -", metrics_json)
    print(" -", preds_csv)
    print(metrics_df.to_string(index=False))

if __name__ == "__main__":
    main()
