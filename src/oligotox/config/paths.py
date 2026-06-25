# config/paths.py
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent.parent

DATA_DIR = ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"

ARTIFACTS_DIR = ROOT / "artifacts"
MODEL_ARTIFACTS_DIR = ARTIFACTS_DIR / "models"
METRICS_ARTIFACTS_DIR = ARTIFACTS_DIR / "metrics"
CONFIG_ARTIFACTS_DIR = ARTIFACTS_DIR / "config"

for p in [
    DATA_DIR,
    RAW_DATA_DIR,
    PROCESSED_DATA_DIR,
    ARTIFACTS_DIR,
    MODEL_ARTIFACTS_DIR,
    METRICS_ARTIFACTS_DIR,
    CONFIG_ARTIFACTS_DIR,
]:
    p.mkdir(parents=True, exist_ok=True)