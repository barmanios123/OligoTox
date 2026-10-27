from __future__ import annotations

from pathlib import Path
from typing import Dict, Any
import pandas as pd

from oligotox.models.loader import (
    load_dataset_metadata,
    load_feature_template,
    load_endpoint_models,
    load_endpoint_evidence,
    load_endpoint_features,
    load_thresholds,
    load_composite_config,
)
from oligotox.pipelines.composite import (
    weighted_composite_score,
    top_risk_endpoint,
    summarize_endpoint,
    probability_to_band,
)


def build_feature_row(
    feature_template: pd.DataFrame,
    candidate_features: Dict[str, Any],
    trained_feature_columns: list[str] | None = None,
) -> pd.DataFrame:
    row = feature_template.copy().astype(object)

    for col in row.columns:
        val = candidate_features.get(col, 0.0)
        if val is None:
            val = 0.0
        row.loc[0, col] = val

    for col in row.columns:
        row[col] = pd.to_numeric(row[col], errors="coerce").fillna(0.0)

    if trained_feature_columns is not None:
        for col in trained_feature_columns:
            if col not in row.columns:
                row[col] = 0.0
        row = row.reindex(columns=trained_feature_columns, fill_value=0.0)

    return row


def score_candidate(
    candidate_features: Dict[str, Any],
    processed_dir: str | Path,
    model_dir: str | Path,
    config_dir: str | Path,
    endpoints: list[str] | None = None,
) -> Dict[str, Any]:
    processed_dir = Path(processed_dir)
    model_dir = Path(model_dir)
    config_dir = Path(config_dir)

    thresholds = load_thresholds(config_dir)
    composite_cfg = load_composite_config(config_dir)
    evidence = load_endpoint_evidence(config_dir)
    per_endpoint_features = load_endpoint_features(config_dir)

    # The shipped config carries the trained columns and endpoints, so scoring
    # works from a fresh clone without the processed training data.
    trained_feature_columns = composite_cfg.get("trained_feature_columns")
    label_cols = endpoints or composite_cfg.get("label_columns")
    if trained_feature_columns is None or label_cols is None:
        meta = load_dataset_metadata(processed_dir)
        trained_feature_columns = trained_feature_columns or list(load_feature_template(processed_dir).columns)
        label_cols = label_cols or meta["label_columns"]
    feature_template = pd.DataFrame([[0.0] * len(trained_feature_columns)], columns=trained_feature_columns)

    models = load_endpoint_models(model_dir, label_cols)

    X_row = build_feature_row(
        feature_template=feature_template,
        candidate_features=candidate_features,
        trained_feature_columns=trained_feature_columns,
    )

    endpoint_evidence = evidence["endpoints"]
    endpoint_probs = {}
    per_endpoint = {}

    for endpoint, model in models.items():
        # cns_toxicity is trained on external sequence-only data, so it has its
        # own column list rather than the shared 55.
        own = per_endpoint_features.get(endpoint, {}).get("feature_columns")
        if own:
            missing = [c for c in own if c not in candidate_features]
            if missing:
                raise KeyError(
                    f"{endpoint} needs features the caller did not supply: {missing}. "
                    "Build them with oligotox.features.candidate.candidate_features()."
                )
            row = build_feature_row(
                feature_template=pd.DataFrame([[0.0] * len(own)], columns=own),
                candidate_features=candidate_features,
                trained_feature_columns=own,
            )
        else:
            row = X_row
        proba = model.predict_proba(row); prob = float(proba[0, 1] if proba.shape[1] > 1 else proba[0, 0])
        endpoint_probs[endpoint] = prob
        thr = thresholds.get(endpoint, {}).get("threshold")
        summary = summarize_endpoint(
            endpoint=endpoint,
            prob=prob,
            threshold=thr,
            bands=composite_cfg["risk_bands"],
        )
        # An endpoint the evaluation doesn't support is reported, but never as a
        # prediction: see scripts/derive_evidence_tiers.py.
        ev = endpoint_evidence.get(endpoint, {})
        tier = ev.get("tier", "insufficient_data")
        summary["evidence_tier"] = tier
        summary["evidence_note"] = ev.get("reason")
        summary["cv_pr_auc"] = ev.get("pr_auc")
        summary["base_rate"] = ev.get("base_rate")
        summary["positive_sequences"] = ev.get("positive_sequences")
        if tier != "validated":
            summary["predicted_positive"] = None
        per_endpoint[endpoint] = summary

    # The composite only averages endpoints that survived the baselines, so it
    # can't dilute the one supported endpoint with eleven unsupported ones.
    basis = [e for e in endpoint_probs if per_endpoint[e]["evidence_tier"] == "validated"]
    validated_probs = {e: endpoint_probs[e] for e in basis}

    composite_score = (
        weighted_composite_score(
            endpoint_probs=validated_probs,
            endpoint_weights=composite_cfg["endpoint_weights"],
        )
        if validated_probs
        else None
    )
    dominant = top_risk_endpoint(validated_probs)

    tiers: Dict[str, list[str]] = {}
    for endpoint, info in per_endpoint.items():
        tiers.setdefault(info["evidence_tier"], []).append(endpoint)

    return {
        "composite_score": composite_score,
        "composite_band": (
            probability_to_band(composite_score, composite_cfg["risk_bands"])
            if composite_score is not None
            else None
        ),
        "composite_basis": basis,
        "dominant_endpoint": dominant,
        "per_endpoint": per_endpoint,
        "endpoints_by_tier": {tier: sorted(v) for tier, v in sorted(tiers.items())},
        "feature_columns": trained_feature_columns,
    }
