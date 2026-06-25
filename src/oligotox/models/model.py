from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import joblib
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

DEFAULT_ENDPOINTS = [
    "hepatocellular_injury",
    "hepatic_vacuolation",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "kupffer_cell_activation",
    "cns_toxicity",
    "offtarget_gene_silencing",
    "injection_site_reaction",
]


@dataclass
class EndpointModelResult:
    endpoint: str
    auc_cv_mean: float | None
    n_samples: int
    positive_rate: float
    fitted_model: Pipeline


class EndpointToxModel:
    def __init__(
        self,
        endpoints: List[str] | None = None,
        n_estimators: int = 200,
        calibration: str = "sigmoid",
        random_state: int = 42,
    ) -> None:
        self.endpoints = endpoints or DEFAULT_ENDPOINTS
        self.n_estimators = n_estimators
        self.calibration = calibration
        self.random_state = random_state
        self.models: Dict[str, Pipeline] = {}
        self.feature_columns: List[str] = []
        self.training_summary: Dict[str, Dict] = {}

    def _build_pipeline(self) -> Pipeline:
        base_rf = RandomForestClassifier(
            n_estimators=self.n_estimators,
            random_state=self.random_state,
            class_weight="balanced",
            min_samples_leaf=1,
            n_jobs=-1,
        )
        return Pipeline([
            ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
            ("model", base_rf),
        ])

    def fit(self, df: pd.DataFrame, feature_columns: List[str]) -> Dict[str, EndpointModelResult]:
        self.feature_columns = feature_columns[:]
        X = df[self.feature_columns].copy()
        results: Dict[str, EndpointModelResult] = {}

        for endpoint in self.endpoints:
            if endpoint not in df.columns:
                continue

            y = df[endpoint].astype(int)
            if y.nunique() < 2:
                continue

            pipe = self._build_pipeline()
            pipe.fit(X, y)
            self.models[endpoint] = pipe

            result = EndpointModelResult(
                endpoint=endpoint,
                auc_cv_mean=None,
                n_samples=int(len(y)),
                positive_rate=float(y.mean()),
                fitted_model=pipe,
            )
            results[endpoint] = result
            self.training_summary[endpoint] = {
                "auc_cv_mean": None,
                "n_samples": int(len(y)),
                "positive_rate": float(y.mean()),
                "use_calibration": False,
            }

        return results

    def predict_proba(self, X: pd.DataFrame) -> pd.DataFrame:
        out = pd.DataFrame(index=X.index)
        for endpoint, model in self.models.items():
            out[endpoint] = model.predict_proba(X[self.feature_columns])[:, 1]
        return out

    def save(self, output_dir: str | Path) -> None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        joblib.dump(self.models, output_dir / "endpoint_models.joblib")
        joblib.dump(self.feature_columns, output_dir / "feature_columns.joblib")
        joblib.dump(self.training_summary, output_dir / "training_summary.joblib")

    def load(self, output_dir: str | Path) -> None:
        output_dir = Path(output_dir)
        self.models = joblib.load(output_dir / "endpoint_models.joblib")
        self.feature_columns = joblib.load(output_dir / "feature_columns.joblib")
        self.training_summary = joblib.load(output_dir / "training_summary.joblib")
