import pandas as pd
from pathlib import Path

path = Path('data/processed/training_dataset.csv')
df = pd.read_csv(path)
print('shape:', df.shape)
print('columns:', list(df.columns))
print('\nmissing values:')
print(df.isna().sum().sort_values(ascending=False).head(20))
