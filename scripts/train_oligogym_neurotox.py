#!/usr/bin/env python3
"""
train_oligogym_neurotox.py

Trains a classifier on MOE_Neurotox (OligoGym) in-vivo ASO safety data.
- 2437 ASOs tested in mice (FOB score 0-7)
- Binary label: y >= 3 = neurotoxic (based on Hagedorn lab threshold)
- Features: sequence k-mers + length + GC content
- Uses RF with 5-fold stratified CV
- Output: artifacts/models_neurotox/ + reports/neurotox_eval.csv

Run: PYTHONPATH="$HOME:$PWD:$PWD/src" python3 scripts/train_oligogym_neurotox.py
"""
import pickle, warnings
from pathlib import Path
from itertools import product
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings('ignore')

DIR = Path('data/oligogym')
MODEL_DIR = Path('artifacts/models_neurotox')
REPORT_DIR = Path('reports')
MODEL_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(exist_ok=True)

TOX_THRESHOLD = 3  # FOB score >= 3 = toxic

def extract_seq(fasta):
    if pd.isna(fasta): return ''
    lines = str(fasta).strip().split('\n')
    return ''.join(l for l in lines if not l.startswith('>')).upper()

def kmer_features_fast(seq, k=2):
    """Fast k-mer frequency dict for k=2 only for speed."""
    nts = list('ACGTU')
    kmers = [''.join(p) for p in product(nts, repeat=k)]
    n = len(seq)
    if n < k: return {f'k{k}_{km}': 0.0 for km in kmers}
    return {f'k{k}_{km}': seq.count(km)/max(n-k+1,1) for km in kmers}

def build_features(df):
    records = []
    for seq in df['sequence']:
        seq = str(seq).upper()
        n = len(seq)
        gc = (seq.count('G') + seq.count('C')) / max(n, 1)
        feats = {'length_nt': n, 'gc_content': gc}
        feats.update(kmer_features_fast(seq, k=2))
        feats.update(kmer_features_fast(seq, k=3))
        # Positional thirds
        t = max(n//3, 1)
        for i, lbl in enumerate(['gc_5p','gc_mid','gc_3p']):
            seg = seq[i*t:(i+1)*t]
            feats[lbl] = (seg.count('G')+seg.count('C')) / max(len(seg),1)
        records.append(feats)
    return pd.DataFrame(records).fillna(0)

def main():
    print('Loading MOE_Neurotox dataset...')
    df = pd.read_csv(DIR / 'MOE_Neurotox_1_processed.csv.gz')
    print(f'  Shape: {df.shape}')
    
    df['sequence'] = df['fasta'].apply(extract_seq)
    df['y_bin'] = (df['y_raw'] >= TOX_THRESHOLD).astype(int)
    
    n_pos = df['y_bin'].sum()
    n_neg = (df['y_bin']==0).sum()
    print(f'  Labels: {n_pos} toxic / {n_neg} safe (threshold={TOX_THRESHOLD})')
    print(f'  Positive rate: {n_pos/len(df)*100:.1f}%')
    
    print('\nBuilding k-mer features...')
    X = build_features(df)
    y = df['y_bin'].values
    print(f'  Feature matrix: {X.shape}')
    
    print('\nTraining RandomForest (5-fold CV)...')
    rf = RandomForestClassifier(
        n_estimators=300,
        max_depth=5,
        min_samples_leaf=5,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    skv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    probas_rf = cross_val_predict(rf, X, y, cv=skv, method='predict_proba')[:, 1]
    auroc_rf = roc_auc_score(y, probas_rf)
    ap_rf = average_precision_score(y, probas_rf)
    brier_rf = brier_score_loss(y, probas_rf)
    print(f'  RF   AUROC={auroc_rf:.3f}  AP={ap_rf:.3f}  Brier={brier_rf:.3f}')
    
    print('\nTraining Logistic Regression (5-fold CV)...')
    lr = Pipeline([
        ('scaler', StandardScaler()),
        ('clf', LogisticRegression(C=0.1, class_weight='balanced', max_iter=1000, solver='saga'))
    ])
    probas_lr = cross_val_predict(lr, X, y, cv=skv, method='predict_proba')[:, 1]
    auroc_lr = roc_auc_score(y, probas_lr)
    ap_lr = average_precision_score(y, probas_lr)
    brier_lr = brier_score_loss(y, probas_lr)
    print(f'  LR   AUROC={auroc_lr:.3f}  AP={ap_lr:.3f}  Brier={brier_lr:.3f}')
    
    # Fit final RF + save
    rf.fit(X, y)
    with open(MODEL_DIR / 'moe_neurotox_rf.pkl', 'wb') as f:
        pickle.dump({'model': rf, 'feature_names': X.columns.tolist()}, f)
    print(f'\nSaved model to {MODEL_DIR}/moe_neurotox_rf.pkl')
    
    # Feature importance
    fi = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)
    print('\nTop 20 features:')
    print(fi.head(20).to_string())
    fi.to_csv(REPORT_DIR / 'neurotox_feature_importance.csv', header=['importance'])
    
    # Save eval
    results = pd.DataFrame([
        {'endpoint':'moe_neurotox_invivo','n':len(y),'n_pos':int(n_pos),'n_neg':int(n_neg),
         'threshold':TOX_THRESHOLD,'model':'RandomForest',
         'auroc':round(auroc_rf,3),'ap':round(ap_rf,3),'brier':round(brier_rf,3)},
        {'endpoint':'moe_neurotox_invivo','n':len(y),'n_pos':int(n_pos),'n_neg':int(n_neg),
         'threshold':TOX_THRESHOLD,'model':'LogisticRegression',
         'auroc':round(auroc_lr,3),'ap':round(ap_lr,3),'brier':round(brier_lr,3)},
    ])
    results.to_csv(REPORT_DIR / 'neurotox_eval.csv', index=False)
    print(f'\nSaved eval to {REPORT_DIR}/neurotox_eval.csv')
    print('\n=== SUMMARY ===')
    print(results.to_string(index=False))

if __name__ == '__main__':
    main()
