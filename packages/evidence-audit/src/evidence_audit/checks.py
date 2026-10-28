"""The checks themselves, independent of any particular dataset.

Each one corresponds to a failure mode in docs/NEGATIVE_RESULTS.md and
answers a question that a PR-AUC on its own cannot:

  leakage          is any feature a relabelled copy of the target?
  grouping         how many independent subjects are behind the row count?
  source_baseline  does a model that sees only the source do as well?
  within_source    does anything survive inside a single source?
  simple_baselines does one crude number do as well as the whole matrix?
  calibration      are the probabilities probabilities, or only ranks?

Nothing here imports the rest of OligoTox. The only inputs are a feature
frame, a binary label, a grouping key per row, and optionally a source per
row.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Sequence

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

LEAKAGE_CORRELATION = 0.9
MIN_POSITIVE_GROUPS = 10
N_BOOT = 1000
# A Brier skill of zero means the probabilities are worth no more than quoting
# the base rate. Testing against zero exactly puts the verdict where the noise
# is: the same dataset scored +0.002 on one machine and -0.009 on another, and
# flipped between "passes" and "ranks but uncalibrated". A calibration claim
# needs a margin it can survive, not a sign.
MIN_BRIER_SKILL = 0.02


def _model():
    base = LogisticRegression(max_iter=2000, class_weight="balanced", solver="liblinear")
    return CalibratedClassifierCV(base, method="sigmoid", cv=2)


def _normalise(values: Sequence[Any], nucleotide: bool = False) -> np.ndarray:
    """Group keys as comparable strings.

    `nucleotide` additionally folds U onto T, so RNA and DNA spellings of one
    sequence count as one group. It is off by default because the same
    substitution silently corrupts any other kind of key: a subject id of
    "USER-1" would become "TSER-1" and stop matching itself.
    """
    out = pd.Series(values).astype(str).str.strip().str.upper()
    if nucleotide:
        out = out.str.replace("U", "T", regex=False)
    return out.values


@dataclass
class EndpointAudit:
    endpoint: str
    rows: int
    unique_groups: int
    positives: int
    positive_groups: int
    base_rate: float
    leaky_features: list[str] = field(default_factory=list)
    grouped_pr_auc: float | None = None
    ungrouped_pr_auc: float | None = None
    source_only_pr_auc: float | None = None
    best_simple_pr_auc: float | None = None
    simple_baseline_used: str | None = None
    within_source: dict[str, float] = field(default_factory=dict)
    brier: float | None = None
    brier_skill: float | None = None
    verdict: str = "unknown"
    findings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class AuditReport:
    endpoints: list[EndpointAudit]

    def to_dict(self) -> dict:
        return {"endpoints": [e.to_dict() for e in self.endpoints]}

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame([
            {
                "endpoint": e.endpoint, "verdict": e.verdict,
                "rows": e.rows, "groups": e.unique_groups,
                "pos_groups": e.positive_groups, "base_rate": e.base_rate,
                "grouped": e.grouped_pr_auc, "ungrouped": e.ungrouped_pr_auc,
                "source_only": e.source_only_pr_auc, "simple": e.best_simple_pr_auc,
                "brier_skill": e.brier_skill,
            }
            for e in self.endpoints
        ])


def find_leaky_features(X: pd.DataFrame, y: np.ndarray) -> list[str]:
    """Features that are near-copies of the target.

    A curated annotation of the outcome correlates with it far beyond what a
    design feature can. This is what turned a 0.93 into a 0.29.
    """
    leaky = []
    for col in X.columns:
        values = pd.to_numeric(X[col], errors="coerce")
        if values.nunique(dropna=True) <= 1:
            continue
        corr = abs(np.corrcoef(values.fillna(values.mean()), y)[0, 1])
        if np.isfinite(corr) and corr >= LEAKAGE_CORRELATION:
            leaky.append(col)
    return leaky


def _oof(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray | None, n_splits: int):
    if groups is None:
        from sklearn.model_selection import StratifiedKFold

        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
        return cross_val_predict(_model(), X, y, cv=cv, method="predict_proba")[:, 1]
    cv = StratifiedGroupKFold(n_splits=n_splits)
    return cross_val_predict(_model(), X, y, cv=cv, groups=groups,
                             method="predict_proba")[:, 1]


def audit_endpoint(
    X: pd.DataFrame,
    y: Sequence[int],
    groups: Sequence[Any],
    sources: Sequence[Any] | None = None,
    endpoint: str = "endpoint",
    simple_features: Sequence[str] = ("gc_content", "seq_length"),
    nucleotide_groups: bool = False,
) -> EndpointAudit:
    """Run every check on one endpoint."""
    y = np.asarray(y).astype(int)
    groups = _normalise(groups, nucleotide_groups)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    positive_groups = len(set(groups[y == 1]))
    audit = EndpointAudit(
        endpoint=endpoint, rows=len(y), unique_groups=len(set(groups)),
        positives=int(y.sum()), positive_groups=positive_groups,
        base_rate=float(y.mean()),
    )

    audit.leaky_features = find_leaky_features(X, y)
    if audit.leaky_features:
        audit.findings.append(
            f"{len(audit.leaky_features)} feature(s) correlate >= {LEAKAGE_CORRELATION} "
            f"with the target and are probably annotations of it: "
            f"{audit.leaky_features[:5]}"
        )

    if audit.unique_groups < audit.rows:
        audit.findings.append(
            f"{audit.rows} rows cover {audit.unique_groups} unique groups; "
            "treat the smaller number as the sample size"
        )

    if positive_groups < 2 or (len(y) - y.sum()) < 2:
        audit.verdict = "not evaluable"
        audit.findings.append("too few positives to cross-validate")
        return audit

    n_splits = int(min(5, positive_groups))
    clean = X.drop(columns=audit.leaky_features, errors="ignore")
    if clean.shape[1] == 0:
        audit.verdict = "not evaluable"
        audit.findings.append("no features left after removing leaky ones")
        return audit

    prob = _oof(clean, y, groups, n_splits)
    audit.grouped_pr_auc = float(average_precision_score(y, prob))
    audit.brier = float(brier_score_loss(y, prob))
    base_brier = float(np.mean((y - y.mean()) ** 2))
    audit.brier_skill = float(1 - audit.brier / base_brier) if base_brier > 0 else None

    # The same thing without grouping, to show what the split is worth.
    ungrouped = _oof(clean, y, None, n_splits)
    audit.ungrouped_pr_auc = float(average_precision_score(y, ungrouped))

    if sources is not None:
        src = pd.Series(sources).astype(str).values
        dummies = pd.get_dummies(src, dtype=float)
        if dummies.shape[1] > 1:
            sp = _oof(dummies, y, groups, n_splits)
            audit.source_only_pr_auc = float(average_precision_score(y, sp))
            for name in sorted(set(src)):
                mask = src == name
                if y[mask].sum() >= 2 and (~y[mask].astype(bool)).sum() >= 2:
                    audit.within_source[name] = float(y[mask].mean())

    available = [c for c in simple_features if c in clean.columns]
    for feature in available:
        sp = _oof(clean[[feature]], y, groups, n_splits)
        score = float(average_precision_score(y, sp))
        if audit.best_simple_pr_auc is None or score > audit.best_simple_pr_auc:
            audit.best_simple_pr_auc = score
            audit.simple_baseline_used = feature

    audit.verdict = _verdict(audit)
    return audit


def _verdict(a: EndpointAudit) -> str:
    if a.positive_groups < MIN_POSITIVE_GROUPS:
        a.findings.append(
            f"{a.positive_groups} positive groups; below {MIN_POSITIVE_GROUPS} "
            "a PR-AUC is not a usable estimate"
        )
        return "insufficient data"
    if a.source_only_pr_auc is not None and a.source_only_pr_auc >= (a.grouped_pr_auc or 0):
        a.findings.append(
            f"a model seeing only the source scores {a.source_only_pr_auc:.2f} against "
            f"{a.grouped_pr_auc:.2f}; the model is recognising the dataset"
        )
        return "source confounded"
    if a.best_simple_pr_auc is not None and a.best_simple_pr_auc >= (a.grouped_pr_auc or 0):
        a.findings.append(
            f"{a.simple_baseline_used} alone scores {a.best_simple_pr_auc:.2f} against "
            f"{a.grouped_pr_auc:.2f}; the feature matrix adds nothing"
        )
        return "no better than a simple baseline"
    if a.brier_skill is not None and a.brier_skill < MIN_BRIER_SKILL:
        a.findings.append(
            f"Brier skill {a.brier_skill:+.3f}, under the {MIN_BRIER_SKILL} needed to "
            "claim calibration; it ranks but its probabilities are no better than "
            "quoting the base rate, so report it as a ranking"
        )
        return "ranks but uncalibrated"
    return "passes"


def audit_dataset(
    X: pd.DataFrame,
    Y: pd.DataFrame | pd.Series,
    groups: Sequence[Any],
    sources: Sequence[Any] | None = None,
    simple_features: Sequence[str] = ("gc_content", "seq_length"),
    nucleotide_groups: bool = False,
) -> AuditReport:
    """Audit every label column in Y."""
    if isinstance(Y, pd.Series):
        Y = Y.to_frame()
    reports = []
    for endpoint in Y.columns:
        y = pd.to_numeric(Y[endpoint], errors="coerce").fillna(0).astype(int).values
        reports.append(audit_endpoint(X, y, groups, sources, endpoint,
                                      simple_features, nucleotide_groups))
    return AuditReport(endpoints=reports)
