"""Assign each endpoint an evidence tier from the baseline results.

The app ships 13 endpoint models, but the baselines show that only one of them
predicts anything a simple baseline doesn't. This script turns that evaluation
into a config the app and the report read at runtime, so the console can't
present an unsupported endpoint as a prediction.

Tiers, applied in order:

  insufficient_data  Fewer than MIN_POSITIVE_SEQUENCES positive sequences, or
                     not evaluated at all. PR-AUC on 4-7 positives is noise.
  source_confounded  The data-source-only baseline matches or beats the model,
                     so the model is recognising which dataset a row came from.
  weak               Clears the source baseline, but the model's 95% CI does
                     not clear the best simple baseline (GC-only, length-only),
                     so it is not distinguishable from one crude number.
  validated          Clears all of the above, and beats its own base rate
                     inside a single source dataset.

Reads reports/results/baselines.csv, within_source.csv and
training_summary_after.json. Writes artifacts/config/endpoint_evidence.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "reports" / "results"
CONFIG = ROOT / "artifacts" / "config"

MODEL = "Logistic (55 features)"
SIMPLE_BASELINES = ("GC content only", "Length only")
SOURCE_BASELINE = "Data source only"
MIN_POSITIVE_SEQUENCES = 10

TIER_ORDER = ("validated", "weak", "source_confounded", "insufficient_data")


def classify(model, source_pr_auc, best_simple, base_rate, n_sequences, within_beats_base):
    """Return (tier, reason, checks) for one endpoint."""
    checks = {
        "has_enough_sequences": n_sequences >= MIN_POSITIVE_SEQUENCES,
        "beats_source_baseline": model["pr_auc"] > source_pr_auc,
        "ci_clears_base_rate": model["ci_low"] > base_rate,
        "ci_clears_simple_baselines": model["ci_low"] > best_simple,
        "beats_base_rate_within_source": within_beats_base,
    }

    if not checks["has_enough_sequences"]:
        return (
            "insufficient_data",
            f"{n_sequences} positive sequences; PR-AUC is not a usable estimate below "
            f"{MIN_POSITIVE_SEQUENCES}",
            checks,
        )
    if not checks["beats_source_baseline"]:
        return (
            "source_confounded",
            f"the data-source-only baseline scores {source_pr_auc:.2f} against the model's "
            f"{model['pr_auc']:.2f}, so the model is recognising the source dataset, not the design",
            checks,
        )
    if not checks["ci_clears_simple_baselines"]:
        return (
            "weak",
            f"the model's 95% CI reaches {model['ci_low']:.2f}, below the best simple baseline "
            f"({best_simple:.2f}), so it is not distinguishable from one crude sequence number",
            checks,
        )
    if not checks["ci_clears_base_rate"]:
        return (
            "weak",
            f"the model's 95% CI reaches {model['ci_low']:.2f}, below the {base_rate:.3f} base rate",
            checks,
        )
    if not within_beats_base:
        return (
            "weak",
            "the signal does not survive inside a single source dataset",
            checks,
        )
    return (
        "validated",
        f"PR-AUC {model['pr_auc']:.2f} [{model['ci_low']:.2f}-{model['ci_high']:.2f}] against a "
        f"{base_rate:.3f} base rate, clearing the source, GC and length baselines and holding up "
        "within a single source",
        checks,
    )


def main() -> None:
    baselines = pd.read_csv(RESULTS / "baselines.csv")
    within = pd.read_csv(RESULTS / "within_source.csv")
    summary = json.loads((RESULTS / "training_summary_after.json").read_text())

    endpoints = {}
    for endpoint, sub in baselines.groupby("endpoint"):
        by_model = sub.set_index("model")
        model = by_model.loc[MODEL]
        source_pr_auc = float(by_model.loc[SOURCE_BASELINE, "pr_auc"])
        best_simple = max(float(by_model.loc[name, "pr_auc"]) for name in SIMPLE_BASELINES)
        base_rate = float(model["base_rate"])
        n_sequences = int(model["positive_sequences"])

        rows = within[(within.endpoint == endpoint) & (within.model == MODEL)]
        within_beats_base = bool((rows.pr_auc > rows.base_rate).any()) if len(rows) else False

        model_metrics = {
            "pr_auc": float(model["pr_auc"]),
            "ci_low": float(model["ci_low"]),
            "ci_high": float(model["ci_high"]),
        }
        tier, reason, checks = classify(
            model_metrics, source_pr_auc, best_simple, base_rate, n_sequences, within_beats_base
        )
        endpoints[endpoint] = {
            "tier": tier,
            "reason": reason,
            "pr_auc": model_metrics["pr_auc"],
            "ci_low": model_metrics["ci_low"],
            "ci_high": model_metrics["ci_high"],
            "base_rate": base_rate,
            "positive_sequences": n_sequences,
            "source_only_pr_auc": source_pr_auc,
            "best_simple_baseline_pr_auc": best_simple,
            "checks": checks,
        }

    # Endpoints the training run could not evaluate never reach baselines.csv.
    for endpoint, stats in summary.items():
        if endpoint in endpoints:
            continue
        n_sequences = int(stats.get("positive_sequences") or 0)
        endpoints[endpoint] = {
            "tier": "insufficient_data",
            "reason": f"{n_sequences} positive sequences; not evaluated",
            "pr_auc": None,
            "positive_sequences": n_sequences,
        }

    validated = sorted(e for e, v in endpoints.items() if v["tier"] == "validated")
    payload = {
        "description": (
            "Evidence tier per endpoint, derived from reports/results/baselines.csv by "
            "scripts/derive_evidence_tiers.py. Only 'validated' endpoints are presented as "
            "predictions and only they enter the composite score."
        ),
        "min_positive_sequences": MIN_POSITIVE_SEQUENCES,
        "tiers": list(TIER_ORDER),
        "composite_endpoints": validated,
        "endpoints": dict(
            sorted(endpoints.items(), key=lambda kv: (TIER_ORDER.index(kv[1]["tier"]), kv[0]))
        ),
    }

    CONFIG.mkdir(parents=True, exist_ok=True)
    out = CONFIG / "endpoint_evidence.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")

    counts = {t: sum(1 for v in endpoints.values() if v["tier"] == t) for t in TIER_ORDER}
    print(f"wrote {out.relative_to(ROOT)}")
    for tier in TIER_ORDER:
        names = sorted(e for e, v in endpoints.items() if v["tier"] == tier)
        print(f"  {tier:18s} {counts[tier]:2d}  {', '.join(names)}")


if __name__ == "__main__":
    main()
