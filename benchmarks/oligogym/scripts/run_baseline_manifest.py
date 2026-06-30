from pathlib import Path
import pandas as pd

manifest = pd.read_csv("benchmarks/oligogym/manifests/gapmer_hepatotox_manifest.csv")
print("Benchmark manifest loaded:")
print(manifest[["dataset_name", "task_type", "endpoint", "status"]].to_string(index=False))

missing = []
for _, row in manifest.iterrows():
    raw_path = Path(row["raw_path"])
    if not raw_path.exists():
        missing.append(str(raw_path))

if missing:
    print("\nMissing raw datasets:")
    for p in missing:
        print(" -", p)
else:
    print("\nAll raw datasets present.")
