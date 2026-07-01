import pandas as pd
import json
import os
from pathlib import Path
from datetime import datetime
import joblib

print("="*80)
print("OLIGOTOX COMPREHENSIVE PLATFORM AUDIT")
print(f"Run time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# 1. DATA FILES AUDIT
print("\n" + "="*80)
print("1. DATA FILES VERIFICATION")
print("="*80)

data_files = {
    'training_dataset': 'data/processed/training_dataset.csv',
    'label_matrix': 'data/processed/label_matrix.csv',
    'feature_matrix': 'data/processed/feature_matrix.csv'
}

for name, path in data_files.items():
    if os.path.exists(path):
        df = pd.read_csv(path)
        mtime = datetime.fromtimestamp(os.path.getmtime(path))
        print(f"\n✓ {name}:")
        print(f"  Shape: {df.shape}")
        print(f"  Last modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')}")
        if name == 'label_matrix':
            endpoints = df.columns[1:].tolist()
            print(f"  Endpoints ({len(endpoints)}): {endpoints}")
            print(f"\n  Class distribution:")
            for col in endpoints:
                pos = int(df[col].sum())
                neg = len(df) - pos
                print(f"    {col:30s} | Pos: {pos:4d} | Neg: {neg:4d} | Rate: {pos/len(df):.3f}")
    else:
        print(f"\n✗ {name}: FILE NOT FOUND at {path}")

# 2. TRAINED MODELS AUDIT
print("\n" + "="*80)
print("2. TRAINED MODELS VERIFICATION")
print("="*80)

model_dir = Path('artifacts/models')
if model_dir.exists():
    model_files = list(model_dir.glob('*.joblib'))
    print(f"\nFound {len(model_files)} model files:")
    for i, model_file in enumerate(sorted(model_files), 1):
        try:
            model = joblib.load(model_file)
            mtime = datetime.fromtimestamp(os.path.getmtime(model_file))
            print(f"{i:2d}. {model_file.stem:30s} | Modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')} | Type: {type(model).__name__}")
        except Exception as e:
            print(f"{i:2d}. {model_file.stem:30s} | ERROR: {str(e)}")
else:
    print("\n✗ Models directory not found")

# 3. TRAINING METRICS AUDIT
print("\n" + "="*80)
print("3. TRAINING METRICS VERIFICATION")
print("="*80)

metrics_file = 'artifacts/metrics/training_summary.json'
if os.path.exists(metrics_file):
    with open(metrics_file) as f:
        metrics = json.load(f)
    print(f"\nTraining metrics for {len(metrics)} endpoints:")
    print(f"\n{'Endpoint':<30s} {'Samples':>8s} {'Pos Rate':>9s} {'PR-AUC':>8s} {'Brier':>8s} {'Status':<20s}")
    print("-"*100)
    for ep, data in sorted(metrics.items()):
        pos_rate = data.get('positive_rate', 0)
        pr_auc = data.get('pr_auc', 'N/A')
        brier = data.get('brier', 'N/A')
        status = data.get('status', 'unknown')
        samples = data.get('samples', 0)
        
        pr_auc_str = f"{pr_auc:.3f}" if isinstance(pr_auc, float) else str(pr_auc)
        brier_str = f"{brier:.3f}" if isinstance(brier, float) else str(brier)
        
        print(f"{ep:<30s} {samples:>8d} {pos_rate:>9.3f} {pr_auc_str:>8s} {brier_str:>8s} {status:<20s}")
else:
    print("\n✗ Training metrics file not found")

# 4. DOCUMENTATION CHECK
print("\n" + "="*80)
print("4. DOCUMENTATION VERIFICATION")
print("="*80)

docs = [
    'docs/TRAINING_PIPELINE_ARCHITECTURE.md',
    'README.md'
]

for doc in docs:
    if os.path.exists(doc):
        mtime = datetime.fromtimestamp(os.path.getmtime(doc))
        with open(doc) as f:
            content = f.read()
        # Check for "19 endpoints" mention
        has_19 = '19' in content and 'endpoint' in content.lower()
        print(f"\n✓ {doc}:")
        print(f"  Last modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')}")
        if has_19:
            print(f"  ⚠️  Contains '19 endpoints' - may need update to reflect actual count")
    else:
        print(f"\n✗ {doc}: NOT FOUND")

# 5. FEATURE ENGINEERING CHECK
print("\n" + "="*80)
print("5. FEATURE ENGINEERING VERIFICATION")
print("="*80)

if os.path.exists('data/processed/feature_matrix.csv'):
    df_features = pd.read_csv('data/processed/feature_matrix.csv')
    print(f"\nFeature matrix shape: {df_features.shape}")
    print(f"Feature columns: {df_features.shape[1] - 1} (excluding oligo_id)")
    
    # Count feature types
    feature_types = {
        'chemical': 0,
        'sequence': 0,
        'structure': 0,
        'embedding': 0,
        'mechanistic': 0,
        'other': 0
    }
    
    for col in df_features.columns[1:]:
        col_lower = col.lower()
        if any(x in col_lower for x in ['gc', 'length', 'mod', 'chemical']):
            feature_types['chemical'] += 1
        elif any(x in col_lower for x in ['purine', 'pyrimidine', 'cpg', 'motif', 'sequence']):
            feature_types['sequence'] += 1
        elif any(x in col_lower for x in ['structure', 'fold', 'mfe']):
            feature_types['structure'] += 1
        elif any(x in col_lower for x in ['embed', 'rna', 'lm']):
            feature_types['embedding'] += 1
        elif any(x in col_lower for x in ['immune', 'tlr', 'seed', 'offtarget']):
            feature_types['mechanistic'] += 1
        else:
            feature_types['other'] += 1
    
    print(f"\nFeature breakdown:")
    for ftype, count in feature_types.items():
        print(f"  {ftype:15s}: {count:4d}")

# 6. DATA CONSISTENCY CHECK
print("\n" + "="*80)
print("6. DATA CONSISTENCY VERIFICATION")
print("="*80)

if all(os.path.exists(data_files[k]) for k in ['training_dataset', 'label_matrix']):
    df_train = pd.read_csv('data/processed/training_dataset.csv')
    df_labels = pd.read_csv('data/processed/label_matrix.csv')
    
    checks = [
        ("Same number of rows", df_train.shape[0] == df_labels.shape[0]),
        ("Expected 1614 oligonucleotides", df_train.shape[0] == 1614),
        ("No missing values in labels", not df_labels.isnull().any().any()),
        ("Labels are binary (0/1)", set(df_labels.iloc[:, 1:].values.flatten()) <= {0, 1, 0.0, 1.0})
    ]
    
    print("\nConsistency checks:")
    for check_name, result in checks:
        status = "✓" if result else "✗"
        print(f"  {status} {check_name}")

# 7. SUMMARY
print("\n" + "="*80)
print("7. SUMMARY & RECOMMENDATIONS")
print("="*80)

if os.path.exists('data/processed/label_matrix.csv'):
    df_labels = pd.read_csv('data/processed/label_matrix.csv')
    actual_endpoints = len(df_labels.columns) - 1
    
    print(f"\n📊 Current Platform State:")
    print(f"   - Oligonucleotides: {df_labels.shape[0]}")
    print(f"   - Endpoints: {actual_endpoints}")
    print(f"   - Features: {df_features.shape[1] - 1 if 'df_features' in locals() else 'N/A'}")
    
    if os.path.exists(metrics_file):
        with open(metrics_file) as f:
            metrics = json.load(f)
        print(f"   - Trained models: {len(metrics)}")
        
        # Calculate average performance
        valid_pr_aucs = [m['pr_auc'] for m in metrics.values() if isinstance(m.get('pr_auc'), float)]
        if valid_pr_aucs:
            avg_pr_auc = sum(valid_pr_aucs) / len(valid_pr_aucs)
            print(f"   - Average PR-AUC: {avg_pr_auc:.3f}")
    
    print(f"\n⚠️  Action Items:")
    # Check if docs mention 19 endpoints
    if os.path.exists('docs/TRAINING_PIPELINE_ARCHITECTURE.md'):
        with open('docs/TRAINING_PIPELINE_ARCHITECTURE.md') as f:
            if '19' in f.read():
                print(f"   1. Update TRAINING_PIPELINE_ARCHITECTURE.md to reflect {actual_endpoints} endpoints")
    
    print(f"\n✅ Platform Status: OPERATIONAL")

print("\n" + "="*80)
print("AUDIT COMPLETE")
print("="*80)
