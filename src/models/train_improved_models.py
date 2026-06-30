#!/usr/bin/env python3
"""
Improved Model Training Pipeline
Uses XGBoost + SMOTE to handle class imbalance and improve performance
"""

import pandas as pd
import numpy as np
import joblib
import json
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, precision_recall_fscore_support
from imblearn.over_sampling import SMOTE
import xgboost as xgb
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_PATH = Path('data/processed')
MODEL_PATH = Path('artifacts/models/latest_train_baseline')
MODEL_PATH.mkdir(exist_ok=True, parents=True)

# Endpoints to train
ENDPOINTS = [
    'hepatocellular_injury', 'hepatic_vacuolation', 'renal_tubular_injury',
    'renal_accumulation', 'thrombocytopenia', 'coagulation_effect',
    'complement_activation', 'innate_immune_activation', 'cytokine_release',
    'kupffer_cell_activation', 'cns_toxicity', 'offtarget_gene_silencing',
    'injection_site_reaction', 'body_weight_loss', 'mortality',
    'proinflammatory_cytokines', 'tlr_activation', 'hybridization_related',
    'developmental_toxicity'
]

def main():
    print("\n" + "="*70)
    print("  OLIGOTOX IMPROVED MODEL TRAINING PIPELINE")
    print("  Algorithm: XGBoost + SMOTE")
    print("="*70 + "\n")
    
    # Load data
    print("Loading data...")
    X = pd.read_csv(DATA_PATH / 'training_dataset.csv')
    if 'oligo_id' in X.columns:
        X = X.drop('oligo_id', axis=1)
    
    y = pd.read_csv(DATA_PATH / 'label_matrix.csv')
    if 'oligo_id' in y.columns:
        y = y.drop('oligo_id', axis=1)
    
    print(f"Loaded {len(X)} samples with {len(X.columns)} features")
    print(f"Training {len(ENDPOINTS)} toxicity endpoints\n")
    
    results = []
    
    for endpoint in ENDPOINTS:
        if endpoint not in y.columns:
            continue
        
        print(f"\n{'='*60}")
        print(f"Training: {endpoint}")
        print(f"{'='*60}")
        
        y_endpoint = y[endpoint]
        base_rate = y_endpoint.mean()
        n_positive = y_endpoint.sum()
        
        print(f"Positive cases: {n_positive}/{len(y)} ({base_rate:.2%})")
        
        if base_rate < 0.05:
            severity = 'High'
        elif base_rate < 0.25:
            severity = 'Moderate'
        else:
            severity = 'Low'
        
        if n_positive < 10:
            print(f"⚠️  Skipping - insufficient positive cases")
            results.append({'endpoint': endpoint, 'status': 'skipped', 'n_positive': int(n_positive)})
            continue
        
        # Train/test split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_endpoint, test_size=0.2, stratify=y_endpoint, random_state=42
        )
        
        print(f"Train: {len(X_train)} ({y_train.sum()} pos) | Test: {len(X_test)} ({y_test.sum()} pos)")
        
        # Apply SMOTE for extreme imbalance
        X_train_resampled = X_train
        y_train_resampled = y_train
        smote_applied = False
        
        if base_rate < 0.15 and n_positive >= 20:
            try:
                print("Applying SMOTE...")
                smote = SMOTE(sampling_strategy=0.3, random_state=42, k_neighbors=min(5, n_positive-1))
                X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
                print(f"After SMOTE: {len(X_train_resampled)} samples ({y_train_resampled.sum()} positive)")
                smote_applied = True
            except Exception as e:
                print(f"SMOTE failed: {e}")
        
        # Calculate scale_pos_weight
        n_neg = (y_train_resampled == 0).sum()
        n_pos = (y_train_resampled == 1).sum()
        scale_pos_weight = n_neg / n_pos if n_pos > 0 else 1
        
        # Train XGBoost
        print("Training XGBoost...")
        model = xgb.XGBClassifier(
            n_estimators=200,
            max_depth=6,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            scale_pos_weight=scale_pos_weight,
            eval_metric='auc',
            random_state=42,
            n_jobs=-1
        )
        model.fit(X_train_resampled, y_train_resampled, verbose=False)
        
        # Evaluate
        y_pred_proba = model.predict_proba(X_test)[:, 1]
        y_pred = model.predict(X_test)
        
        try:
            auroc = roc_auc_score(y_test, y_pred_proba)
        except:
            auroc = 0.5
        
        precision, recall, f1, _ = precision_recall_fscore_support(
            y_test, y_pred, average='binary', zero_division=0
        )
        
        print(f"\n📊 Results: AUROC={auroc:.4f} | Precision={precision:.4f} | Recall={recall:.4f} | F1={f1:.4f}")
        
        # Save model
        model_path = MODEL_PATH / f'{endpoint}_xgboost.pkl'
        joblib.dump(model, model_path)
        print(f"✅ Saved: {model_path.name}")
        
        # Feature importance
        feature_importance = dict(zip(X.columns, model.feature_importances_))
        top_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)[:5]
        
        results.append({
            'endpoint': endpoint,
            'status': 'success',
            'base_rate': float(base_rate),
            'severity': severity,
            'n_positive': int(n_positive),
            'smote_applied': smote_applied,
            'auroc': float(auroc),
            'precision': float(precision),
            'recall': float(recall),
            'f1_score': float(f1),
            'top_features': [(feat, float(imp)) for feat, imp in top_features]
        })
    
    # Save results
    results_path = MODEL_PATH / 'training_results.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Summary
    successful = [r for r in results if r['status'] == 'success']
    if successful:
        avg_auroc = np.mean([r['auroc'] for r in successful])
        print(f"\n\n{'='*70}")
        print("  TRAINING COMPLETE")
        print(f"{'='*70}")
        print(f"\n📊 Summary:")
        print(f"  Successfully trained: {len(successful)}/{len(results)} endpoints")
        print(f"  Average AUROC: {avg_auroc:.4f}")
        print(f"  SMOTE applied: {sum([r.get('smote_applied', False) for r in successful])} endpoints")
        
        top_5 = sorted(successful, key=lambda x: x['auroc'], reverse=True)[:5]
        print(f"\n🏆 Top 5 Models:")
        for i, r in enumerate(top_5, 1):
            print(f"  {i}. {r['endpoint']}: {r['auroc']:.4f}")
        
        print(f"\n✅ Results: {results_path}")
        print(f"✅ Models: {MODEL_PATH}/\n")

if __name__ == '__main__':
    main()
