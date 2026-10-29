"""Does target-site accessibility improve the per-assay potency models?

`evaluate_accessibility.py` showed accessibility is worth something pooled
across assays: on 282 genes it adds +0.019 to a 0.229 baseline and scores
+0.062 alone. That is a measurement, not a decision. The shipped potency
models are per-context, and a feature earns a place in them by improving
them, inside their own data, under the same gates everything else faces.

**The trap this is built around.** Accessibility needs the target transcript.
Only some genes are cached, and inside a cached gene only some oligos locate
a binding site. Filling the rest with zeros would hand the model a column
that is really "was this site locatable", which tracks gene identity and
cache membership rather than biology. That is exactly the failure in
NEGATIVE_RESULTS entry 13, where zero-filling chemistry turned it into a
source indicator and manufactured a +0.290 gain out of nothing.

So nothing is zero-filled. Each context is restricted to the rows whose site
actually locates, and design-only and design-plus-accessibility are compared
on that identical subset. The comparison is therefore honest and the
population is narrower than the shipped model's: this answers "does
accessibility help where it is computable", which is the only question the
data can answer.

    python scripts/evaluate_context_accessibility.py

Writes reports/results/context_accessibility.{csv,json}.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold, cross_val_predict

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from oligotox.efficacy.potency import guide_features  # noqa: E402
from oligotox.features.accessibility import site_accessibility  # noqa: E402
from train_context_potency import _model, load_contexts  # noqa: E402

GENOMIC = ROOT / "data" / "external" / "genomic"
RESULTS = ROOT / "reports" / "results"
ACCESS_COLS = ["site_accessibility", "site_log_accessibility",
               "site_flank_accessibility", "site_relative_position"]
# Grouped CV needs several genes, and a gain measured on a handful of genes
# is the mistake entry 12 already records.
MIN_ROWS = 300
MIN_GENES = 5
N_SPLITS = 5
# A feature earns its place by clearing a margin, not by being positive.
MIN_GAIN = 0.01
# Each located row costs one RNAplfold fold of an 800nt window, so a full
# pass over every context runs for hours and an earlier attempt was killed
# before it wrote anything. A seeded sample per context keeps the run bounded.
# This costs precision on the gain, not direction, and the sample is large
# enough that a +0.01 effect is still visible.
MAX_ROWS_PER_CONTEXT = 1200


def load_sequences() -> dict[str, str]:
    out = {}
    for path in sorted(GENOMIC.glob("*.json")):
        if path.name.startswith("_"):
            continue
        try:
            record = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        seq = (record.get("sequence") or "").upper().replace("U", "T")
        if seq:
            out[record.get("symbol", path.stem)] = seq
    return out


def _rho(X: pd.DataFrame, y: np.ndarray, genes: np.ndarray, cv) -> float:
    pred = cross_val_predict(_model(), X.astype(float), y, cv=cv, groups=genes)
    value = spearmanr(y, pred).statistic
    return float(value) if np.isfinite(value) else 0.0


def main() -> None:
    transcripts = load_sequences()
    if not transcripts:
        raise SystemExit("no genomic spans cached; run scripts/fetch_transcripts.py --genomic")
    print(f"{len(transcripts)} genomic span(s) cached", flush=True)
    contexts = load_contexts()
    print(f"{len(contexts)} context(s) to test, "
          f"up to {MAX_ROWS_PER_CONTEXT} rows each\n", flush=True)

    rows = []
    for slug, frame, _meta in contexts:
        frame = frame[frame["target_gene"].isin(transcripts)]
        if frame.empty:
            continue
        if len(frame) > MAX_ROWS_PER_CONTEXT:
            frame = frame.sample(MAX_ROWS_PER_CONTEXT, random_state=42)

        # Keep only rows whose site actually locates. Nothing is zero-filled,
        # so no column can encode "this one was not computable".
        located, access = [], []
        for seq, gene in frame[["sequence", "target_gene"]].itertuples(index=False):
            site = site_accessibility(seq, transcripts[gene])
            if site.site_found:
                located.append(True)
                access.append({c: site.as_features()[c] for c in ACCESS_COLS})
            else:
                located.append(False)
        subset = frame[pd.Series(located, index=frame.index)]
        if len(subset) < MIN_ROWS or subset["target_gene"].nunique() < MIN_GENES:
            print(f"  {slug:32s} skipped, {len(subset)} located rows over "
                  f"{subset['target_gene'].nunique()} gene(s)")
            continue

        design = pd.DataFrame([guide_features(s) for s in subset["sequence"]],
                              index=subset.index)
        acc = pd.DataFrame(access, index=subset.index)
        y = subset["y"].values
        genes = subset["target_gene"].values
        cv = GroupKFold(n_splits=min(N_SPLITS, subset["target_gene"].nunique()))

        base = _rho(design, y, genes, cv)
        combined = _rho(pd.concat([design, acc], axis=1), y, genes, cv)
        alone = _rho(acc, y, genes, cv)
        gain = combined - base

        rows.append({
            "context": slug,
            "rows_located": int(len(subset)),
            "rows_sampled": int(len(frame)),
            "sampled": bool(len(frame) == MAX_ROWS_PER_CONTEXT),
            "coverage": round(len(subset) / max(len(frame), 1), 3),
            "genes": int(subset["target_gene"].nunique()),
            "design_only": round(base, 4),
            "design_plus_accessibility": round(combined, 4),
            "accessibility_only": round(alone, 4),
            "gain": round(gain, 4),
            "helps": bool(gain > MIN_GAIN),
        })
        print(f"  {slug:32s} n={len(subset):6,} genes={subset['target_gene'].nunique():3d}  "
              f"design={base:+.3f} +access={combined:+.3f} alone={alone:+.3f} "
              f"gain={gain:+.4f}", flush=True)

    if not rows:
        raise SystemExit("no context had enough located rows to evaluate")

    helped = [r for r in rows if r["helps"]]
    verdict = (f"accessibility clears +{MIN_GAIN} in {len(helped)} of {len(rows)} context(s)"
               if helped else
               "accessibility does not clear the margin in any context")
    print(f"\n  {verdict}")
    if helped:
        print("  Worth wiring into those contexts, for rows where a site locates. "
              "Rows without one must stay unscored rather than zero-filled.")

    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS / "context_accessibility.csv", index=False)
    (RESULTS / "context_accessibility.json").write_text(json.dumps({
        "question": "Does target-site accessibility improve the per-assay potency "
                    "models, inside their own data?",
        "method": "Restricted to rows whose binding site locates in a cached genomic "
                  "span. Nothing zero-filled, so no column encodes computability. "
                  "Target-gene-grouped folds, same model as the shipped contexts.",
        "caveat": "The located subset is narrower than the shipped model's "
                  "population. This answers whether accessibility helps where it "
                  "is computable, not whether it helps everywhere.",
        "min_gain": MIN_GAIN,
        "verdict": verdict,
        "results": rows,
    }, indent=2) + "\n")
    print(f"\nwrote {(RESULTS / 'context_accessibility.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
