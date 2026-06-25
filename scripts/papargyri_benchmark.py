from pathlib import Path
import re
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score

DATA_PATH = Path("data/raw/lna_gapmer/aso_data.csv")

def infer_design_col(df):
    candidates = ["design", "aso_design", "ASO", "pattern", "chemistry"]
    for c in candidates:
        if c in df.columns:
            return c
    return None

def infer_seq_col(df):
    candidates = ["sequence", "seq", "Sequence", "oligo_sequence"]
    for c in candidates:
        if c in df.columns:
            return c
    return None

def infer_target_col(df):
    candidates = ["target", "region", "Target"]
    for c in candidates:
        if c in df.columns:
            return c
    return None

def infer_casp_col(df):
    candidates = ["cas_avg", "caspase_avg", "CaspasePercent", "cas", "caspase"]
    for c in candidates:
        if c in df.columns:
            return c
    return None

def parse_design(design, seq=None):
    s = "" if pd.isna(design) else str(design)

    lna_tokens = len(re.findall(r"LNA|\\+|locked", s, flags=re.I))
    dna_tokens = len(re.findall(r"DNA|deoxy", s, flags=re.I))

    left = s[:max(0, min(8, len(s)))]
    right = s[max(0, len(s)-8):]

    lna_5p = int(bool(re.search(r"LNA|\\+|locked", left, flags=re.I)))
    lna_3p = int(bool(re.search(r"LNA|\\+|locked", right, flags=re.I)))

    if seq is not None and not pd.isna(seq):
        aso_len = len(str(seq).strip())
    else:
        aso_len = len(re.findall(r"[ACGTU]", s, flags=re.I))
        if aso_len == 0:
            aso_len = len(s)

    if dna_tokens == 0 and aso_len >= lna_tokens:
        dna_tokens = max(aso_len - lna_tokens, 0)

    return {
        "lna_5p": lna_5p,
        "lna_3p": lna_3p,
        "lna_count": int(lna_tokens),
        "dna_count": int(dna_tokens),
        "aso_len": int(aso_len),
    }

def main():
    if not DATA_PATH.exists():
        raise SystemExit(f"Missing {DATA_PATH}. Put aso_data.csv there first.")

    df = pd.read_csv(DATA_PATH)

    design_col = infer_design_col(df)
    seq_col = infer_seq_col(df)
    target_col = infer_target_col(df)
    casp_col = infer_casp_col(df)

    print("Detected columns:")
    print(" design_col =", design_col)
    print(" seq_col    =", seq_col)
    print(" target_col =", target_col)
    print(" casp_col   =", casp_col)

    if casp_col is None:
        raise SystemExit("Could not detect caspase column. Rename it to cas_avg or edit the script.")
    if design_col is None and seq_col is None:
        raise SystemExit("Could not detect a design or sequence column.")

    feat_rows = []
    for _, r in df.iterrows():
        feat_rows.append(
            parse_design(
                r[design_col] if design_col else "",
                r[seq_col] if seq_col else None
            )
        )

    feat_df = pd.DataFrame(feat_rows)
    work = pd.concat([df.reset_index(drop=True), feat_df], axis=1)

    work["too_toxic"] = (pd.to_numeric(work[casp_col], errors="coerce") > 300).astype(int)

    feature_cols = ["lna_5p", "lna_3p", "lna_count", "dna_count", "aso_len"]
    work = work.dropna(subset=feature_cols + ["too_toxic"]).copy()

    if target_col and work[target_col].astype(str).str.contains("A").any() and work[target_col].astype(str).str.contains("B").any():
        train = work[work[target_col].astype(str).str.contains("A", na=False)].copy()
        test = work[work[target_col].astype(str).str.contains("B", na=False)].copy()
        split_name = "target A train / target B test"
    else:
        train = work.sample(frac=0.7, random_state=42)
        test = work.drop(train.index).copy()
        split_name = "random 70/30 fallback"

    X_train = train[feature_cols]
    y_train = train["too_toxic"]
    X_test = test[feature_cols]
    y_test = test["too_toxic"]

    print("\\nRows:", len(work))
    print("Train:", len(train), "Test:", len(test))
    print("Split:", split_name)
    print("Positive rate train:", round(y_train.mean(), 4))
    print("Positive rate test :", round(y_test.mean(), 4))

    rf = RandomForestClassifier(
        n_estimators=300,
        random_state=42,
        class_weight="balanced"
    )
    rf.fit(X_train, y_train)
    rf_probs = rf.predict_proba(X_test)[:, 1]

    if y_test.nunique() > 1:
        print("\\nRandomForest benchmark")
        print(" ROC-AUC:", round(roc_auc_score(y_test, rf_probs), 4))
        print(" PR-AUC :", round(average_precision_score(y_test, rf_probs), 4))
    else:
        print("\\nRandomForest benchmark skipped: test set has one class.")

    try:
        import feyn
        ql = feyn.QLattice()
        train_df = X_train.copy()
        train_df["too_toxic"] = y_train.values
        models = ql.auto_run(
            train_df,
            output_name="too_toxic",
            kind="classification",
            max_complexity=7
        )
        best = models[0]
        preds = best.predict(X_test)
        if hasattr(preds, "values"):
            preds = preds.values
        preds = pd.Series(preds).astype(float)

        print("\\nFeyn QLattice benchmark")
        if y_test.nunique() > 1:
            print(" ROC-AUC:", round(roc_auc_score(y_test, preds), 4))
            print(" PR-AUC :", round(average_precision_score(y_test, preds), 4))
        else:
            print(" skipped metrics: test set has one class.")
        print(" Best model:", best)
    except Exception as e:
        print("\\nFeyn benchmark not run:", repr(e))

    out_path = Path("data/processed/papargyri_engineered.csv")
    work.to_csv(out_path, index=False)
    print("\\nWrote engineered dataset to", out_path)

if __name__ == "__main__":
    main()
