from pathlib import Path
#!/usr/bin/env python3
"""OligoTox train_v3.py — per-endpoint RF/LR/GBM classifiers with calibration.
Run: PYTHONPATH=$HOME python oligotox/train_v3.py
"""
import os, json, warnings, pickle
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, cross_validate
warnings.filterwarnings("ignore")

ENDPOINT_COLS = [
    "endpoint_hepatocellular_injury","endpoint_renal_tubular_injury","endpoint_renal_accumulation",
    "endpoint_thrombocytopenia","endpoint_coagulation_effect","endpoint_complement_activation",
    "endpoint_innate_immune_activation","endpoint_cytokine_release","endpoint_hypersensitivity",
    "endpoint_cns_neurotoxicity","endpoint_offtarget_gene_silencing","endpoint_reproductive_dev_tox",
    "endpoint_injection_site_reaction",
]
FEATURE_COLS = [
    "length_nt","gc_content","ps_fraction","lna_fraction","moe_fraction",
    "has_GalNAc","lnp_formulated","black_box_warning","warning_hepatic_monitoring",
    "warning_renal_monitoring","warning_thrombocytopenia","warning_neurological",
]
MODEL_DIR = os.path.expanduser("~/oligotox/artifacts/models_v3")
os.makedirs(MODEL_DIR, exist_ok=True)

def load_data():
    for p in [
        os.path.expanduser("~/oligotox_datasets/processed/curated_oligo_tox_master_v3.csv"),
        os.path.expanduser("~/oligotox/data/raw/curated_oligo_tox_master_v2.csv"),
        os.path.expanduser("~/oligotox/data/raw/curated_oligo_tox_template.csv"),
    ]:
        if os.path.exists(p):
            print(f"Loading: {p}"); return pd.read_csv(p)
        candidate_paths = [
        Path("data/curated/curated_oligo_tox_master_v3.csv"),
        Path("data/curated/curated_oligo_tox_template_v3.csv"),
        Path("downloads/curated_oligo_tox_master_v3.csv"),
        Path("curated_oligo_tox_master_v3.csv"),
    ]
    for csv_path in candidate_paths:
        if csv_path.exists():
            return pd.read_csv(csv_path)
    raise FileNotFoundError("No curated CSV found.")

def build_features(df):
    X = df[[c for c in FEATURE_COLS if c in df.columns]].copy().fillna(0).astype(float)
    for col in ["ASO_gapmer","siRNA","SSO","PMO","mRNA_LNP"]:
        X[f"mod_{col}"] = (df.get("modality","") == col).astype(float)
    for col in ["2MOE_PS","LNA_PS","GalNAc_siRNA","LNP_siRNA","PMO"]:
        X[f"chem_{col}"] = (df.get("chemistry","") == col).astype(float)
    for col in ["naked_sc","naked_iv","naked_it","GalNAc_sc","LNP_iv","LNP_im"]:
        X[f"del_{col}"] = (df.get("delivery","") == col).astype(float)
    return X

def oversample_minority(X, y, seed=42):
    counts = {v: (y==v).sum() for v in [0,1]}
    if counts.get(1,0) < 2 or counts.get(1,0) >= counts.get(0,0)//3:
        return X, y
    target = min(counts[1]*2, counts[0]//3, 20)
    if target <= counts[1]: return X, y
    rng = np.random.default_rng(seed)
    pos_idx = np.where(y==1)[0]
    extras = rng.choice(pos_idx, target - counts[1])
    X = pd.concat([X, X.iloc[extras]], ignore_index=True)
    y = np.concatenate([y, y[extras]])
    return X, y

def select_model(n_pos):
    if n_pos < 4: return None
    if n_pos < 8:
        return CalibratedClassifierCV(
            LogisticRegression(C=0.1, class_weight="balanced", max_iter=500), cv=2)
    if n_pos < 15:
        return CalibratedClassifierCV(
            RandomForestClassifier(n_estimators=300, max_depth=4,
                class_weight="balanced", min_samples_leaf=2, random_state=42), cv=3)
    return CalibratedClassifierCV(
        GradientBoostingClassifier(n_estimators=200, max_depth=3,
            learning_rate=0.05, subsample=0.8, random_state=42), cv=3)

def train():
    df = load_data(); X = build_features(df); report = []
    for ep in ENDPOINT_COLS:
        if ep not in df.columns: continue
        y = df[ep].fillna(0).astype(int).values
        n_pos = int(y.sum())
        estimator = select_model(n_pos)
        row = dict(endpoint=ep, n=len(y), n_positive=n_pos,
                   positive_rate=round(n_pos/len(y),3))
        if estimator is None:
            row.update(status="RULES_ONLY", pr_auc=None, roc_auc=None)
            print(f"  {ep}: RULES_ONLY ({n_pos} pos)")
            report.append(row); continue
        X_ep, y_ep = oversample_minority(X, y)
        cv = StratifiedKFold(n_splits=min(5, n_pos), shuffle=True, random_state=42)
        try:
            res = cross_validate(estimator, X_ep, y_ep, cv=cv,
                scoring=["average_precision","roc_auc","neg_brier_score"], n_jobs=-1)
            pr = float(np.mean(res["test_average_precision"]))
            roc = float(np.mean(res["test_roc_auc"]))
            brier = float(np.mean(-res["test_neg_brier_score"]))
            status = "STRONG" if pr >= 0.65 else ("MODERATE" if pr >= 0.40 else "WEAK")
            estimator.fit(X_ep, y_ep)
            with open(os.path.join(MODEL_DIR, f"{ep}.pkl"), "wb") as f:
                pickle.dump(estimator, f)
            row.update(status=status, pr_auc=round(pr,3), roc_auc=round(roc,3), brier=round(brier,3))
            print(f"  {ep}: {status} PR={pr:.3f} ROC={roc:.3f}")
        except Exception as e:
            row.update(status="ERROR", pr_auc=None, roc_auc=None, note=str(e))
        report.append(row)
    with open(os.path.join(MODEL_DIR, "training_report_v3.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n✓ Models saved to {MODEL_DIR}")

if __name__ == "__main__":
    train()
