from pathlib import Path
from oligotox.integration.external_registry import DATASET_PATHS

def benchmark_status():
    return {
        name: {
            "path": str(path),
            "exists": Path(path).exists(),
        }
        for name, path in DATASET_PATHS.items()
    }
