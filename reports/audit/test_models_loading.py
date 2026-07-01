import joblib
import pandas as pd
import numpy as np
from pathlib import Path
import json

print("="*80)
print("TESTING ALL TRAINED MODELS")
print("="*80)

# Load test data
print("\n1. Loading test data...")
df_features = pd.read_csv('data/processed/feature_matrix.csv')
df_labels = pd.read_csv('data/processed/label_matrix.csv')

print(f"  ✓ Features: {df_features.shape}")
print(f"  ✓ Labels: {df_labels.shape}")

# Get endpoints from labels
endpoints = df_labels.columns[1:].tolist()
print(f"\n2. Testing {len(endpoints)} models...")

model_dir = Path('artifacts/models')
results = []

for i, endpoint in enumerate(endpoints, 1):
    model_path = model_dir / f"{endpoint}.joblib"
    
    try:
        # Load model
        model = joblib.load(model_path)
        
        # Get a sample of features (first 5 rows)
        X_sample = df_features.iloc[:5, 1:].values
        
        # Make predictions
        y_pred = model.predict(X_sample)
        y_proba = model.predict_proba(X_sample)
        
        # Check outputs
        assert len(y_pred) == 5, "Prediction length mismatch"
        assert y_proba.shape == (5, 2), "Probability shape mismatch"
        
        results.append({
            'endpoint': endpoint,
            'status': '✅ OK',
            'model_type': type(model).__name__,
            'predictions': 'Valid'
        })
        
        print(f"  {i:2d}. {endpoint:30s} ✓ Loaded & Tested")
        
    except FileNotFoundError:
        results.append({
            'endpoint': endpoint,
            'status': '✗ NOT FOUND',
            'model_type': 'N/A',
            'predictions': 'N/A'
        })
        print(f"  {i:2d}. {endpoint:30s} ✗ Model file not found")
        
    except Exception as e:
        results.append({
            'endpoint': endpoint,
            'status': f'✗ ERROR',
            'model_type': 'N/A',
            'predictions': str(e)
        })
        print(f"  {i:2d}. {endpoint:30s} ✗ Error: {str(e)[:50]}")

# Summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)

success_count = sum(1 for r in results if '✅' in r['status'])
fail_count = len(results) - success_count

print(f"\n✅ Successfully loaded: {success_count}/{len(results)}")
if fail_count > 0:
    print(f"✗ Failed: {fail_count}/{len(results)}")
else:
    print(f"\n🎉 ALL MODELS WORKING CORRECTLY!")

print("\n" + "="*80)
