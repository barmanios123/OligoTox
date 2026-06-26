from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXTERNAL = ROOT / "external"

MODEL_PATHS = {
    "pfred": EXTERNAL / "models" / "PFRED",
    "oligograph": EXTERNAL / "models" / "OligoGraph",
    "sirnadiscovery": EXTERNAL / "models" / "siRNADiscovery",
    "asoptimizer": EXTERNAL / "models" / "ASOptimizer",
}

DATASET_PATHS = {
    "oligogym": EXTERNAL / "datasets" / "OligoGym",
}
