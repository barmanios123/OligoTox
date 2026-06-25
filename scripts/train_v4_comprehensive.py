#!/usr/bin/env python3
"""
OligoTox train_v4_comprehensive.py

Upgrades over v3:
  - Reads data/processed/training_dataset_with_chem.csv (796 rows, 207 cols)
  - Adds k-mer (2-mer, 3-mer) and position features from sequence
  - Per-endpoint classifiers: LR / RF / GBM auto-selected by n_pos
  - Stratified 5-fold CV: AUROC, AP, Brier score
  - Feature importance (RF/GBM) or coefficient (LR)
  - Exports reports/eval_v4_summary.csv + reports/eval_v4_feature_importance.csv

Run from repo root:
  PYTHONPATH="$HOME:$PWD:$PWD/src" python3 scripts/train_v4_comprehensive.py
"""
import os, json, warnings, pickle
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
warnings.filterwarnings("ignore")

DATA_PATH = Path("data/processed/training_dataset_with_chem.csv")
MODEL_DIR = Path("artifacts/models_v4")
REPORT_DIR = Path("reports")
MODEL_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

FEATURE_COLS = [
    "length_nt","gc_content","ps_fraction","lna_fraction","moe_fraction",
    "has_GalNAc","lnp_formulated","black_box_warning",
    "warning_hepatic_monitoring","warning_renal_monitoring",
    "warning_thrombocytopenia","warning_neurological",
]

NUCLEOTIDES = list("ACGTU")
DINUCS     = ["".join(p) for p in product(NUCLEOTIDES, repeat=2)]
TRINUCS    = ["".join(p) for p in product(NUCLEOTIDES, repeat=3)]

def kmer_features(seq):
    """Return dict of normalised 2-mer and 3-mer frequencies from a sequence."""
    seq = str(seq).upper().replace(" ","").replace("-","")
    n = len(seq)
    feats = {}
    for di in DINUCS:
        feats[f"di_{di}"] = seq.count(di) / max(n-1, 1)
    for tri in TRINUCS:
        feats[f"tri_{tri}"] = seq.count(tri) / max(n-2, 1)
    # positional GC in thirds
    thirds = max(n//3, 1)
    for i, label in enumerate(["gc_5p","gc_mid","gc_3p"]):
        seg = seq[i*thirds:(i+1)*thirds]
        feats[label] = (seg.count("G")+seg.count("C")) / max(len(seg),1)
    return feats

def load_data():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Cannot find {DATA_PATH}")
    print(f"Loading: {DATA_PATH}")
    df = pd.read_csv(DATA_PATH)
    print(f"  Shape: {df.shape}")
    return df

def build_features(df):
    # Base chemistry features
    X = df[[c for c in FEATURE_COLS if c in df.columns]].copy().fillna(0).astype(float)
    # Modality dummies
    for col in ["ASO_gapmer","siRNA","SSO","PMO","mRNA_LNP"]:
        X[f"mod_{col}"] = (df.get("modality","") == col).astype(float)
    for col in ["2MOE_PS","LNA_PS","GalNAc_siRNA","LNA_siRNA","PMO"]:
        X[f"chem_{col}"] = (df.get("chemistry","") == col).astype(float)
    for col in ["naked_sc","naked_iv","naked_it","GalNAc_sc","LNP_iv","LNP_im"]:
        X[f"del_{col}"] = (df.get("delivery","") == col).astype(float)
    # K-mer features from sequence
    seq_col = None
    for c in ["sequence","Sequence","seq","nt_sequence"]:
        if c in df.columns:
            seq_col = c; break
    if seq_col:
        kmer_df = pd.DataFrame([kmer_features(s) for s in df[seq_col]]).fillna(0)
        X = pd.concat([X.reset_index(drop=True), kmer_df.reset_index(drop=True)], axis=1)
        print(f"  Added {kmer_df.shape[1]} k-mer features")
    else:
        print("  WARNING: No sequence column found; skipping k-mer features")
    return X

def select_model(n_pos):
    if n_pos < 4: return None
    if n_pos < 10:
        return Pipeline([
            ("scaler", StandardScaler()),
            ("clf", CalibratedClassifierCV(
                LogisticRegression(C=0.01, class_weight="balanced", max_iter=1000, solver="liblinear"),
                cv=2))
        ])
    if n_pos < 30:
        return CalibratedClassifierCV(
            RandomForestClassifier(n_estimators=200, max_depth=4,
                class_weight="balanced", min_samples_leaf=3, random_state=42),
            cv=3)
    return CalibratedClassifierCV(
        GradientBoostingClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
            subsample=0.8, min_samples_leaf=5, random_state=42),
        cv=3)

def get_feature_importance(model, feature_names):
    """Extract feature importance or coef from the model."""
    try:
        if hasattr(model, "named_steps"):
            est = model.named_steps.get("clf", model)
        else:
            est = model
        # CalibratedClassifierCV wraps the base estimator
        if hasattr(est, "calibrated_classifiers_"):
            base = est.calibrated_classifiers_[0].estimator
        elif hasattr(est, "base_estimator"):
            base = est.base_estimator
        else:
            base = est
        if hasattr(base, "feature_importances_"):
            return dict(zip(feature_names, base.feature_importances_))
        if hasattr(base, "coef_"):
            return dict(zip(feature_names, np.abs(base.coef_[0])))
    except Exception:
        pass
    return {}

def main():
    df = load_data()
    X = build_features(df)
    feature_names = X.columns.tolist()
    print(f"  Feature matrix: {X.shape}")

    endpoint_cols = [c for c in df.columns if c.startswith("endpoint_")]
    print(f"  {len(endpoint_cols)} endpoint columns found")

    results = []
    feat_importance_rows = []
    skv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    for ep in endpoint_cols:
        y_raw = df[ep].fillna(0)
        y = y_raw.astype(int)
        n_pos = int(y.sum())
        n_neg = int((y==0).sum())
        print(f"\n--- {ep}: pos={n_pos}, neg={n_neg} ---")

        model = select_model(n_pos)
        if model is None:
            print(f"  Skipping: too few positives ({n_pos})")
            results.append({
                "endpoint": ep, "n": len(y), "n_pos": n_pos, "n_neg": n_neg,
                "auroc": None, "ap": None, "brier": None,
                "model": "skipped", "cv_folds": 0
            })
            continue

        # Determine CV folds available
        max_cv = min(5, n_pos)
        if max_cv < 2:
            print(f"  Skipping CV: not enough positive folds")
            continue
        skv_ep = StratifiedKFold(n_splits=max_cv, shuffle=True, random_state=42)

        try:
            probas = cross_val_predict(model, X, y, cv=skv_ep, method="predict_proba")[:,1]
            auroc = roc_auc_score(y, probas) if n_pos > 1 else None
            ap    = average_precision_score(y, probas) if n_pos > 1 else None
            brier = brier_score_loss(y, probas)
            print(f"  AUROC={auroc:.3f}  AP={ap:.3f}  Brier={brier:.3f}")
        except Exception as e:
            print(f"  CV failed: {e}")
            auroc = ap = brier = None

        # Fit full model for saving + feature importance
        try:
            model.fit(X, y)
            fi = get_feature_importance(model, feature_names)
            top = sorted(fi.items(), key=lambda x: -x[1])[:10]
            for rank,(feat,imp) in enumerate(top,1):
                feat_importance_rows.append({
                    "endpoint": ep, "rank": rank, "feature": feat, "importance": round(imp,5)
                })
        except Exception as e:
            print(f"  Model fit failed: {e}")

        model_name = type(model).__name__
        results.append({
            "endpoint": ep, "n": len(y), "n_pos": n_pos, "n_neg": n_neg,
            "auroc": round(auroc,3) if auroc else None,
            "ap": round(ap,3) if ap else None,
            "brier": round(brier,3) if brier else None,
            "model": model_name, "cv_folds": max_cv
        })

        # Save model
        ep_short = ep.replace("endpoint_","")
        with open(MODEL_DIR / f"{ep_short}.pkl", "wb") as f:
            pickle.dump(model, f)

    # Save reports
    res_df = pd.DataFrame(results).sort_values("auroc", ascending=False, na_position="last")
    res_df.to_csv(REPORT_DIR / "eval_v4_summary.csv", index=False)
    print(f"\nSaved: {REPORT_DIR}/eval_v4_summary.csv")

    if feat_importance_rows:
        fi_df = pd.DataFrame(feat_importance_rows)
        fi_df.to_csv(REPORT_DIR / "eval_v4_feature_importance.csv", index=False)
        print(f"Saved: {REPORT_DIR}/eval_v4_feature_importance.csv")

    print("\n=== EVALUATION SUMMARY ===")
    print(res_df[["endpoint","n_pos","auroc","ap","brier","model"]].to_string(index=False))

    valid = res_df.dropna(subset=["auroc"])
    if len(valid):
        print(f"\nMean AUROC across {len(valid)} endpoints: {valid['auroc'].mean():.3f}")
        print(f"Median AUROC: {valid['auroc'].median():.3f}")
        print(f"Best:  {valid.iloc[0]['endpoint']} = {valid.iloc[0]['auroc']:.3f}")
        print(f"Worst: {valid.iloc[-1]['endpoint']} = {valid.iloc[-1]['auroc']:.3f}")

if __name__ == "__main__":
    main()
