import json

import pandas as pd

BASELINES = {
    "Data source only",
    "GC content only",
    "Length only",
}



def test_baselines_cover_every_evaluated_endpoint(root):
    summary = json.loads((root / "reports/results/training_summary_after.json").read_text())
    evaluated = {e for e, s in summary.items() if s["status"] == "ok"}
    df = pd.read_csv(root / "reports/results/baselines.csv")
    assert set(df.endpoint) == evaluated
    for endpoint, sub in df.groupby("endpoint"):
        # The two model labels carry the feature count, which changes with
        # the matrix, so resolve them rather than hard-coding a number.
        named = {m for m in sub.model if m.startswith(("Logistic (", "Gradient boosting ("))}
        assert len(named) == 2, endpoint
        assert set(sub.model) - named == BASELINES, endpoint


def test_confidence_intervals_are_valid(root):
    for name in ("baselines.csv", "within_source.csv"):
        df = pd.read_csv(root / "reports/results" / name)
        assert (df.ci_low <= df.ci_high).all()
        assert df[["pr_auc", "ci_low", "ci_high"]].stack().between(0, 1).all()
