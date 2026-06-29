#!/usr/bin/env python3
"""
Load literature-derived hepatic toxicity examples and integrate into training data.

Usage:
    python scripts/load_literature_examples.py data/literature_examples/your_file.csv
"""

import sys
import pandas as pd
from pathlib import Path

def load_literature_examples(csv_path):
    """
    Load literature examples and prepare for training.
    
    Expected columns:
    - sequence: oligonucleotide sequence
    - modality: ASO, siRNA, etc.
    - chemistry: PS, 2OME, MOE, etc.
    - conjugate: GalNAc, None, etc.
    - hepatic_vacuolation: 0 or 1
    - kupffer_cell_activation: 0 or 1
    - source_citation: literature reference
    """
    df = pd.read_csv(csv_path)
    
    required_cols = ['sequence', 'modality', 'chemistry', 'hepatic_vacuolation', 'kupffer_cell_activation']
    missing = [c for c in required_cols if c not in df.columns]
    
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    print(f"Loaded {len(df)} examples from {csv_path}")
    print(f"  Hepatic vacuolation: {df['hepatic_vacuolation'].sum()} positive")
    print(f"  Kupffer activation: {df['kupffer_cell_activation'].sum()} positive")
    
    return df

def featurize_examples(df):
    """Convert sequences to features matching score_candidate template."""
    from oligotox.features.featurizer import score_candidate
    
    features_list = []
    
    for idx, row in df.iterrows():
        try:
            # Create input matching score_candidate signature
            inp = {
                'sequence': row['sequence'],
                'modality': row['modality'],
                'chemistry': row['chemistry'],
                'conjugate': row.get('conjugate', 'None'),
                'length': len(row['sequence']),
                'dose_mg_kg': row.get('dose_mg_kg', 3.0)
            }
            
            features = score_candidate(inp)
            features['hepatic_vacuolation'] = row['hepatic_vacuolation']
            features['kupffer_cell_activation'] = row['kupffer_cell_activation']
            features['source'] = row.get('source_citation', 'literature')
            
            features_list.append(features)
            
        except Exception as e:
            print(f"Warning: Failed to featurize row {idx}: {e}")
            continue
    
    return pd.DataFrame(features_list)

def append_to_training_data(new_df, output_path='data/training_data_with_literature.csv'):
    """Append new examples to existing training data."""
    
    # Load existing training data if it exists
    if Path(output_path).exists():
        existing = pd.read_csv(output_path)
        combined = pd.concat([existing, new_df], ignore_index=True)
        print(f"Appended to existing data. Total rows: {len(combined)}")
    else:
        combined = new_df
        print(f"Created new training file with {len(combined)} rows")
    
    combined.to_csv(output_path, index=False)
    print(f"Saved to: {output_path}")
    
    return output_path

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python scripts/load_literature_examples.py <csv_file>")
        print("Example: python scripts/load_literature_examples.py data/literature_examples/hepatic_toxicity_template.csv")
        sys.exit(1)
    
    csv_path = sys.argv[1]
    
    # Load
    df = load_literature_examples(csv_path)
    
    # Featurize
    print("\nFeaturizing sequences...")
    features_df = featurize_examples(df)
    print(f"Created feature matrix: {features_df.shape}")
    
    # Append to training data
    print("\nIntegrating with training data...")
    output_path = append_to_training_data(features_df)
    
    print("\n✅ Done! Next steps:")
    print("1. Review the combined data:", output_path)
    print("2. Retrain models: python train_models.py")
    print("3. Validate: pytest tests/test_training.py")
