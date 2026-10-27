from __future__ import annotations

from pathlib import Path
import json
import joblib
import pandas as pd


def _read_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return json.loads(path.read_text())


def load_dataset_metadata(processed_dir: str | Path) -> dict:
    processed_dir = Path(processed_dir)
    return _read_json(processed_dir / "dataset_metadata.json")


def load_feature_columns(processed_dir: str | Path) -> list[str]:
    processed_dir = Path(processed_dir)
    feature_csv = processed_dir / "feature_matrix.csv"
    if not feature_csv.exists():
        raise FileNotFoundError(f"Missing feature matrix: {feature_csv}")
    X = pd.read_csv(feature_csv, nrows=1)
    return list(X.columns)


def load_feature_template(processed_dir: str | Path) -> pd.DataFrame:
    cols = load_feature_columns(processed_dir)
    return pd.DataFrame([[0.0] * len(cols)], columns=cols)


def load_thresholds(config_dir: str | Path) -> dict:
    config_dir = Path(config_dir)
    return _read_json(config_dir / "endpoint_thresholds.json")


def load_composite_config(config_dir: str | Path) -> dict:
    config_dir = Path(config_dir)
    return _read_json(config_dir / "composite_config.json")


def load_endpoint_features(config_dir: str | Path) -> dict:
    """Per-endpoint feature columns, for models not using the shared list.

    Missing file means every endpoint uses the shared trained columns.
    """
    config_dir = Path(config_dir)
    path = config_dir / "endpoint_features.json"
    return json.loads(path.read_text()) if path.exists() else {}


def load_endpoint_evidence(config_dir: str | Path) -> dict:
    """Evidence tier per endpoint, from scripts/derive_evidence_tiers.py."""
    config_dir = Path(config_dir)
    return _read_json(config_dir / "endpoint_evidence.json")


def load_endpoint_models(model_dir: str | Path, endpoints: list[str]) -> dict[str, object]:
    model_dir = Path(model_dir)
    out = {}
    for endpoint in endpoints:
        path = model_dir / f"{endpoint}.joblib"
        if path.exists():
            out[endpoint] = joblib.load(path)
    if not out:
        raise RuntimeError(f"No endpoint models found in {model_dir}")
    return out
