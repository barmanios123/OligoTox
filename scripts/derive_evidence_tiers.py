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

MODEL_PREFIX = "Logistic ("
SIMPLE_BASELINES = ("GC content only", "Length only")
SOURCE_BASELINE = "Data source only"
MIN_POSITIVE_SEQUENCES = 10
# Calibration has to clear a margin, not merely beat zero. See the note where
# this is applied: a gate at exactly zero flips on noise.
MIN_BRIER_SKILL = 0.02

TIER_ORDER = ("validated", "weak", "source_confounded", "insufficient_data")


def _model_label(df) -> str:
    """The model label carries the feature count, so match on the prefix."""
    return next(m for m in df.model.unique() if m.startswith(MODEL_PREFIX))

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
        model = by_model.loc[_model_label(sub)]
        source_pr_auc = float(by_model.loc[SOURCE_BASELINE, "pr_auc"])
        best_simple = max(float(by_model.loc[name, "pr_auc"]) for name in SIMPLE_BASELINES)
        base_rate = float(model["base_rate"])
        n_sequences = int(model["positive_sequences"])

        rows = within[(within.endpoint == endpoint) & (within.model.str.startswith(MODEL_PREFIX))]
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

    # Endpoints validated against an independent external assay instead of the
    # internal baselines. A tier here needs the transfer CI to clear chance.
    external = RESULTS / "cns_external.json"
    if external.exists():
        ext = json.loads(external.read_text())
        transfers = [r for r in ext["results"]
                     if r["evaluation"] == "cross-assay transfer"
                     and r["features"].startswith("All")]
        if transfers and all(r["ci_low"] > 0.5 for r in transfers):
            worst = min(transfers, key=lambda r: r["roc_auc"])
            best = max(transfers, key=lambda r: r["roc_auc"])
            total = sum(d["sequences"] for d in ext["datasets"].values())
            rule = [r for r in ext["results"]
                    if r["evaluation"] == "cross-assay transfer" and r["features"].startswith("G ")]
            beats_rule = all(
                t["roc_auc"] > g["roc_auc"]
                for t, g in zip(sorted(transfers, key=lambda r: r["assay"]),
                                sorted(rule, key=lambda r: r["assay"]))
            ) if rule else True
            caveat = ("" if beats_rule else
                      " It does not improve on the published G-count rule, which scores the same "
                      "or better, so treat it as that rule rather than a new finding.")
            endpoints["cns_toxicity"] = {
                "tier": "validated",
                "reason": (
                    f"trained on {total} sequences from 2 independent assays and validated by "
                    f"training on one and testing on the other, with no shared sequences "
                    f"(ROC-AUC {worst['roc_auc']:.2f}-{best['roc_auc']:.2f}, "
                    f"CI low {worst['ci_low']:.2f})." + caveat
                ),
                "pr_auc": float(best["pr_auc"]),
                "ci_low": float(worst["ci_low"]),
                "ci_high": float(best["ci_high"]),
                "base_rate": float(best["base_rate"]),
                "positive_sequences": int(sum(d["positives"] for d in ext["datasets"].values())),
                "validation": "cross-assay external",
                "checks": {"external_transfer_above_chance": True,
                           "beats_published_rule": bool(beats_rule)},
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

    # Calibration is separate from ranking: an endpoint can beat its baselines
    # on PR-AUC while its probabilities are worse than quoting the base rate.
    calibration_path = RESULTS / "calibration.csv"
    if calibration_path.exists():
        import csv

        with calibration_path.open() as fh:
            for row in csv.DictReader(fh):
                entry = endpoints.get(row["endpoint"])
                if not entry:
                    continue
                # Rounded: BLAS differences between platforms move these in the
                # 5th decimal, which would otherwise rewrite the config on every
                # run and bury real changes in noise.
                skill = round(float(row["brier_skill"]), 4)
                entry["brier_skill"] = skill
                entry["calibration_mae"] = round(float(row["calibration_mae"]), 4)
                # A bare skill > 0 is not a gate. Hepatocellular injury came in
                # at +0.0017 after an unrelated feature change, which is a 0.2%
                # improvement on quoting the base rate, and that was enough to
                # flip its label from "ranking only" to "probability". A label
                # that turns over in the third decimal is not a label, so a
                # margin is required.
                entry["probabilities_calibrated"] = skill >= MIN_BRIER_SKILL
                if skill < MIN_BRIER_SKILL:
                    entry["reason"] += (
                        " Its probabilities are not calibrated: Brier skill "
                        f"{skill:+.4f} against the base rate, below the "
                        f"{MIN_BRIER_SKILL} needed to be worth reading as a "
                        "percentage, so treat the score as a ranking."
                    )

    # A directional external check, too small to change a tier but not to report.
    hepatotox = RESULTS / "hepatotox_external.json"
    if hepatotox.exists():
        ext = json.loads(hepatotox.read_text())
        entry = endpoints.get("hepatocellular_injury")
        if entry:
            entry["external_check"] = {
                "source": ext["source"], "sequences": ext["sequences"],
                "positives": ext["positives"], "roc_auc": ext["roc_auc"],
                "roc_ci": ext["roc_ci"], "status": ext["status"],
            }
            if ext["roc_auc"] < 0.5:
                entry["reason"] += (
                    f" On an independent screen of {ext['sequences']} LNA gapmers it ranks "
                    f"the {ext['positives']} hepatotoxic ones no better than chance "
                    f"(ROC-AUC {ext['roc_auc']:.2f}, CI {ext['roc_ci'][0]:.2f}-"
                    f"{ext['roc_ci'][1]:.2f}). Too few positives to be conclusive, but "
                    "there is no evidence this endpoint transfers off its own source."
                )

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
