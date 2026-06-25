from pathlib import Path
import json
from datetime import datetime

registry = Path("data/versions/dataset_registry.json")
payload = json.loads(registry.read_text()) if registry.exists() else []
payload.append({
    "timestamp": datetime.utcnow().isoformat(),
    "training_dataset": "data/processed/training_dataset.csv",
    "feature_matrix": "data/processed/feature_matrix.csv",
    "label_matrix": "data/processed/label_matrix.csv",
})
registry.write_text(json.dumps(payload, indent=2))
print("Registered dataset version")
