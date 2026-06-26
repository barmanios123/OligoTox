#!/usr/bin/env python3
"""
OligoTox Model Training Script
Trains multi-endpoint toxicity prediction models
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import joblib
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("🧠 OLIGOTOX MODEL TRAINING")
print("="*70)

# Load training dataset
print("\n📊 Loading training dataset...")
df = pd.read_csv('data/processed/training_dataset.csv')
print(f"✅ Loaded {len(df)} sequences with {len(df.columns)} features")

# Define endpoints to train
endpoints = [
    'endpoint_cns_delayed_toxicity',
    'endpoint_cns_neurotoxicity', 
    'endpoint_coagulation_effect',
    'endpoint_complement_activation',
    'endpoint_cytokine_release',
    'endpoint_hepatocellular_injury',
    'endpoint_hepatic_vacuolation',
    'endpoint_injection_site_reaction',
    'endpoint_innate_immune_activation',
    'endpoint_kupffer_cell_activation',
    'endpoint_renal_accumulation',
    'endpoint_renal_tubular_injury',
    'endpoint_thrombocytopenia'
]

# Prepare features - use numeric features only
feature_cols = ['gc_content', 'cas_avg', 'kd_avg', 'ps_fraction', 'lna_fraction']
X = df[feature_cols].fillna(0)

print(f"\n🔧 Training features: {feature_cols}")
print(f"📈 Training {len(endpoints)} endpoint models...\n")

# Create models directory
model_dir = Path('src/oligotox/models/trained')
model_dir.mkdir(parents=True, exist_ok=True)

# Train a model for each endpoint
for i, endpoint in enumerate(endpoints, 1):
    if endpoint not in df.columns:
        print(f"⚠️  {i}/{len(endpoints)}: {endpoint} - NOT FOUND IN DATA")
        continue
    
    # Prepare target
    y = df[endpoint].fillna(0).astype(int)
    
    # Skip if no positive examples
    if y.sum() == 0:
        print(f"⚠️  {i}/{len(endpoints)}: {endpoint} - NO POSITIVE EXAMPLES")
        continue
    
    # Train model
    model = RandomForestClassifier(
        n_estimators=100,
        max_depth=10,
        random_state=42,
        n_jobs=-1do e
    )
    
    model.fit(X, y)
    
    # Save model
    model_path = model_dir / f"{endpoint}.joblib"
    joblib.dump(model, model_path)
    
    # Calculate accuracy
    train_acc = model.score(X, y)
    n_pos = y.sum()
    
    print(f"✅ {i:2d}/{len(endpoints)}: {endpoint:40s} (acc={train_acc:.3f}, n_pos={n_pos:4d})")

print("\n" + "="*70)
print("✨ MODEL TRAINING COMPLETE!")
print(f"📁 Models saved to: {model_dir}")

# List trained models
trained_models = list(model_dir.glob('*.joblib'))
print(f"\n🎯 Total models trained: {len(trained_models)}")

for model_file in trained_models:
    size_kb = model_file.stat().st_size / 1024
    print(f"   - {model_file.name} ({size_kb:.1f} KB)")

print("\n🚀 Platform ready for predictions!")
