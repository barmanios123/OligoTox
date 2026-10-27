"""Train the shipped siRNA potency model.

Validated by scripts/evaluate_efficacy.py: Spearman 0.53 under gene-grouped
CV on 2,431 siRNAs, and 0.44 / 0.48 transferring between two independent
studies. The target-gene-mean baseline carries none of it.

The ASO model from McQuisten is trained too but tiered `weak`: 0.25 is real
but modest, and it is a 2007 activity scale that does not map onto a modern
potency readout.

Writes artifacts/models/efficacy_<name>.joblib and
artifacts/config/efficacy_models.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from evaluate_efficacy import ASO, SIRNA, _feature_cols, _model, build, guide_features  # noqa: F401

ROOT = Path(__file__).resolve().parent.parent
MODELS = ROOT / "artifacts" / "models"
CONFIG = ROOT / "artifacts" / "config"
RESULTS = ROOT / "reports" / "results"


def main() -> None:
    results = json.loads((RESULTS / "efficacy.json").read_text())["results"]

    def stat(dataset, model, evaluation):
        for r in results:
            if r["dataset"] == dataset and r["model"] == model and r["evaluation"] == evaluation:
                return r
        return None

    # siRNA: pool both Ichihara studies, which share no sequences.
    sirna = pd.concat([build(s) for s in SIRNA], ignore_index=True)
    cols = _feature_cols(sirna)
    sirna_model = _model().fit(sirna[cols].astype(float), sirna["y_raw"].values)

    aso = build(ASO)
    aso_cols = _feature_cols(aso)
    aso_model = _model().fit(aso[aso_cols].astype(float), aso["y_raw"].values)

    MODELS.mkdir(parents=True, exist_ok=True)
    joblib.dump(sirna_model, MODELS / "efficacy_sirna_potency.joblib")
    joblib.dump(aso_model, MODELS / "efficacy_aso_activity.joblib")

    transfer = [r for r in results if r["evaluation"] == "cross-study transfer"
                and r["model"] == "All guide features"]
    worst = min(transfer, key=lambda r: r["spearman"])

    config = {
        "sirna_potency": {
            "model_file": "efficacy_sirna_potency.joblib",
            "modality": "siRNA",
            "readout": "percent target inhibition",
            "units": "percent knockdown, higher is more potent",
            "feature_columns": cols,
            "tier": "validated",
            "training_sequences": int(len(sirna)),
            "training_genes": int(sirna.target_gene.nunique()),
            "grouped_cv_spearman": stat(SIRNA[0], "All guide features", "gene-grouped CV")["spearman"],
            "cross_study_spearman_low": worst["spearman"],
            "evidence": (
                f"Spearman {stat(SIRNA[0], 'All guide features', 'gene-grouped CV')['spearman']:.2f} "
                f"under target-gene-grouped CV, and {worst['spearman']:.2f} or better transferring "
                "between two independent studies with no shared sequences. It beats GC-only and "
                "asymmetry-only baselines, and a target-gene-mean baseline carries none of it."
            ),
            "caveat": "Trained on 2007 siRNA chemistry with no modification annotation, so it "
                      "scores sequence design only and knows nothing about chemistry.",
        },
        "aso_activity": {
            "model_file": "efficacy_aso_activity.joblib",
            "modality": "ASO",
            "readout": "antisense activity, 0-1",
            "units": "relative activity, higher is more active",
            "feature_columns": aso_cols,
            "tier": "weak",
            "training_sequences": int(len(aso)),
            "training_genes": int(aso.target_gene.nunique()),
            "grouped_cv_spearman": stat(ASO, "All guide features", "gene-grouped CV")["spearman"],
            "evidence": (
                f"Spearman {stat(ASO, 'All guide features', 'gene-grouped CV')['spearman']:.2f} "
                "under target-gene-grouped CV. Above the GC-only baseline but modest, and it has "
                "no second study to transfer to."
            ),
            "caveat": "Shown as a ranking aid only, never as a predicted knockdown.",
        },
    }
    (CONFIG / "efficacy_models.json").write_text(json.dumps(config, indent=2) + "\n")

    for name, spec in config.items():
        oof = spearmanr(
            (sirna if name == "sirna_potency" else aso)["y_raw"].values,
            (sirna_model if name == "sirna_potency" else aso_model).predict(
                (sirna[cols] if name == "sirna_potency" else aso[aso_cols]).astype(float)),
        ).statistic
        print(f"{name:16s} tier={spec['tier']:9s} n={spec['training_sequences']:5d} "
              f"genes={spec['training_genes']:3d} in-sample rho={oof:.3f}")


if __name__ == "__main__":
    main()
