"""What would each failing endpoint need in order to pass?

"Insufficient data" is true but useless on its own. For each endpoint that
does not reach `validated`, this works out the binding constraint and what it
would take to clear it, so the gap becomes a shopping list rather than a
verdict.

Four constraints, in the order the tier rules apply them:

  too few positives      needs MIN_POSITIVE_SEQUENCES positive sequences;
                         the shortfall is exact arithmetic
  source confounded      positives sit in one source, so the source predicts
                         them; needs positives from an independent source
  CI does not clear      the estimate is real but too noisy; a bootstrap CI
  a baseline             narrows roughly as 1/sqrt(n), so the multiple needed
                         is an order-of-magnitude estimate, not a promise
  uncalibrated           ranks without usable probabilities; reports the gap
                         between predicted and observed rather than inventing
                         a sample size

Estimates are labelled as estimates. The arithmetic ones are not.

Writes reports/results/endpoint_requirements.json and prints a summary.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "artifacts" / "config"
RESULTS = ROOT / "reports" / "results"
PROCESSED = ROOT / "data" / "processed"
MIN_POSITIVE_SEQUENCES = 10


def _source_breakdown() -> dict[str, dict[str, int]]:
    """Positive sequences per source, per endpoint."""
    labels = pd.read_csv(PROCESSED / "label_matrix.csv")
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    seq = raw["sequence"].astype(str).str.upper().str.replace("U", "T", regex=False)
    main = ("oligotox_main", "papargyri_2020")
    source = raw["source"].where(raw["source"].isin(main), "curated_other")

    out = {}
    for endpoint in labels.columns:
        # This counts positive sequences per source, so an unmeasured row
        # cannot affect it either way. Written as a mask rather than a
        # zero-fill so that nobody copies the zero-fill out of here into a
        # place where it does matter, which is how it reached the trainer.
        y = pd.to_numeric(labels[endpoint], errors="coerce") == 1
        positives = pd.DataFrame({"seq": seq, "source": source})[y == 1]
        out[endpoint] = positives.groupby("source")["seq"].nunique().to_dict()
    return out


def _needed_for_ci(current_positives: int, ci_low: float, target: float) -> int | None:
    """Rough positives needed for the CI lower bound to clear `target`.

    A bootstrap interval's half-width shrinks about as 1/sqrt(n). This scales
    the current gap by that relationship. It is an order-of-magnitude guide
    and nothing more; it assumes the point estimate holds as data is added,
    which is exactly what more data might disprove.
    """
    if ci_low >= target or current_positives <= 0:
        return None
    gap = target - ci_low
    if gap <= 0:
        return None
    # Needed multiple on the sample size to close the gap, capped so the
    # output stays an estimate rather than a fantasy.
    multiple = min((gap / max(ci_low, 0.01) + 1) ** 2, 50.0)
    return int(math.ceil(current_positives * multiple))


def requirements() -> list[dict]:
    evidence = json.loads((CONFIG / "endpoint_evidence.json").read_text())["endpoints"]
    sources = _source_breakdown()
    baselines = pd.read_csv(RESULTS / "baselines.csv")
    logistic = next(m for m in baselines.model.unique() if m.startswith("Logistic ("))

    rows = []
    for endpoint, info in evidence.items():
        tier = info["tier"]
        entry = {
            "endpoint": endpoint,
            "tier": tier,
            "positive_sequences": info.get("positive_sequences"),
            "blocking_constraint": None,
            "needs": [],
            "estimate": False,
        }

        if tier == "validated" and info.get("probabilities_calibrated") is not False:
            entry["blocking_constraint"] = "none"
            entry["needs"].append("Nothing. It is presented as a prediction.")
            rows.append(entry)
            continue

        positives = info.get("positive_sequences") or 0

        if tier == "insufficient_data":
            entry["blocking_constraint"] = "too few positive sequences"
            shortfall = max(MIN_POSITIVE_SEQUENCES - positives, 0)
            entry["needs"].append(
                f"{shortfall} more positive sequences to reach the minimum of "
                f"{MIN_POSITIVE_SEQUENCES}, which is where a PR-AUC becomes a usable "
                f"estimate. It has {positives}."
            )
            breakdown = sources.get(endpoint, {})
            if len(breakdown) == 1:
                only = next(iter(breakdown))
                entry["needs"].append(
                    f"They should not all come from {only}: positives from a second "
                    "independent source are what separate biology from curation."
                )

        elif tier == "source_confounded":
            entry["blocking_constraint"] = "positives come from one source"
            breakdown = sources.get(endpoint, {})
            entry["source_breakdown"] = breakdown
            biggest = max(breakdown, key=breakdown.get) if breakdown else "one source"
            entry["needs"].append(
                f"Positive sequences from a source other than {biggest}, which holds "
                f"{breakdown.get(biggest, 0)} of its {positives}. While every positive "
                "sits in one bucket, a model that recognises the bucket scores better "
                "than one that reads the molecule."
            )
            entry["needs"].append(
                "A systematic screen rather than hand-curated literature examples, "
                "since curated rows exist because a finding was published about them."
            )

        elif tier == "weak":
            entry["blocking_constraint"] = "estimate too noisy to clear a baseline"
            sub = baselines[baselines.endpoint == endpoint].set_index("model")
            simple = max(float(sub.loc["GC content only", "pr_auc"]),
                         float(sub.loc["Length only", "pr_auc"]))
            ci_low = float(info.get("ci_low") or sub.loc[logistic, "ci_low"])
            needed = _needed_for_ci(positives, ci_low, simple)
            entry["estimate"] = True
            if needed:
                entry["needs"].append(
                    f"Roughly {needed} positive sequences, against {positives} now, for "
                    f"the CI lower bound ({ci_low:.2f}) to clear the best simple "
                    f"baseline ({simple:.2f}). An order-of-magnitude estimate assuming "
                    "the point estimate holds, which more data may disprove."
                )
            else:
                entry["needs"].append(
                    "A narrower confidence interval, which means more positive sequences."
                )

        if info.get("probabilities_calibrated") is False:
            skill = info.get("brier_skill")
            entry["needs"].append(
                f"Calibration: Brier skill {skill:+.2f} means the probabilities are no "
                "better than quoting the base rate. Ranking and calibration are "
                "separate problems, and this one is not fixed by a better model. It "
                "needs more labelled examples across the probability range."
            )
            if entry["blocking_constraint"] is None:
                entry["blocking_constraint"] = "uncalibrated probabilities"

        rows.append(entry)
    return rows


def main() -> None:
    rows = requirements()
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "endpoint_requirements.json").write_text(json.dumps({
        "min_positive_sequences": MIN_POSITIVE_SEQUENCES,
        "note": "Shortfalls against the minimum are exact. Sample sizes needed to "
                "narrow a confidence interval are order-of-magnitude estimates and "
                "are flagged with estimate: true.",
        "endpoints": rows,
    }, indent=2) + "\n")

    order = {"none": 0, "uncalibrated probabilities": 1,
             "estimate too noisy to clear a baseline": 2,
             "positives come from one source": 3, "too few positive sequences": 4}
    for row in sorted(rows, key=lambda r: order.get(r["blocking_constraint"], 9)):
        flag = "  [estimate]" if row["estimate"] else ""
        print(f"\n{row['endpoint']}  [{row['tier']}]  blocked by: "
              f"{row['blocking_constraint']}{flag}")
        for need in row["needs"]:
            print(f"  - {need}")
    print(f"\nwrote {(RESULTS / 'endpoint_requirements.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
