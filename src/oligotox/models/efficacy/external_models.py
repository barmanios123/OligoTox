from pathlib import Path
from oligotox.integration.external_registry import MODEL_PATHS

def external_model_status():
    return {
        name: {
            "path": str(path),
            "exists": Path(path).exists(),
        }
        for name, path in MODEL_PATHS.items()
    }
