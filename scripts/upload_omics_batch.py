from pathlib import Path
import pandas as pd

omics_dir = Path("data/raw/omics")
omics_dir.mkdir(parents=True, exist_ok=True)

for csv_path in omics_dir.glob("*.csv"):
    df = pd.read_csv(csv_path)
    print(csv_path.name, df.shape)
