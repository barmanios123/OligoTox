"""Normalize the external datasets into one labelled table, with the
direction of every readout written down.

The bug this exists to prevent: `data/oligogym/oligogym_tox_combined.csv` has
a column named `y_binary` that actually holds rounded raw readouts (calcium
scores up to 200, percent inhibition down to -27). The companion
`oligogym_tox_combined_fixed.csv` is genuinely binary but gets Hagedorn's
direction backwards, labelling 83% of those compounds neurotoxic. Neither
file records a threshold or a direction, so neither can be checked.

Turning a continuous readout into a toxic/not label is a modeling decision.
Every such decision lives in SOURCES below, with its direction, threshold and
citation, and a source missing from that table is an error rather than a
silent default.

Three of the six datasets measure potency, not toxicity. They are kept and
marked `efficacy`, and excluded from the toxicity output, because a
target-inhibition readout in a file named "tox_combined" is how a potency
number ends up trained as a safety label.

Writes data/processed/external_tox_normalized.csv,
data/processed/external_efficacy_normalized.csv and a manifest recording
every decision.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
COMBINED = ROOT / "data" / "oligogym" / "oligogym_tox_combined.csv"
OUT = ROOT / "data" / "processed"

HIGHER_IS_TOXIC = "higher_is_toxic"
LOWER_IS_TOXIC = "lower_is_toxic"


@dataclass(frozen=True)
class Source:
    readout: str
    kind: str              # "toxicity" or "efficacy"
    endpoint: str
    direction: str | None  # None for efficacy, which gets no toxic label
    threshold: float | None
    citation: str
    note: str = ""


SOURCES: dict[str, Source] = {
    "MOE_Neurotox_2024": Source(
        readout="FOB_score_3h",
        kind="toxicity",
        endpoint="cns_toxicity",
        direction=HIGHER_IS_TOXIC,
        threshold=3.0,
        citation="MOE neurotoxicity in vivo panel, functional observational battery",
        note="FOB scores functional abnormality, so a higher score is worse. Unambiguous.",
    ),
    "Hagedorn_2022": Source(
        readout="calcium_oscillation_score",
        kind="toxicity",
        endpoint="cns_toxicity",
        direction=LOWER_IS_TOXIC,
        threshold=50.0,
        citation="Hagedorn et al., Nucleic Acid Ther. 2022, doi:10.1089/nat.2021.0071",
        note="Acute tolerability correlates with REDUCED spontaneous calcium "
             "oscillation, so a low score is toxic. oligogym_tox_combined_fixed.csv "
             "has this inverted.",
    ),
    "Papargyri_2020_OligoGym": Source(
        readout="caspase_minmax_scaled",
        kind="toxicity",
        endpoint="cytotoxicity",
        direction=HIGHER_IS_TOXIC,
        threshold=0.5,
        citation="Papargyri et al. 2020, caspase activation",
        note="Caspase activation marks apoptosis, so higher is worse. Already "
             "min-max scaled to 0-1. Overlaps the internal training set, so it "
             "must not be used as external validation.",
    ),
    "Ichihara_2007_1": Source(
        readout="pct_inhibition", kind="efficacy", endpoint="target_inhibition",
        direction=None, threshold=None,
        citation="Ichihara et al. 2007",
        note="Potency, not toxicity. Excluded from the toxicity table.",
    ),
    "Ichihara_2007_2": Source(
        readout="pct_inhibition", kind="efficacy", endpoint="target_inhibition",
        direction=None, threshold=None,
        citation="Ichihara et al. 2007",
        note="Potency, not toxicity. Excluded from the toxicity table.",
    ),
    "McQuisten_2007": Source(
        readout="activity_range_0to1", kind="efficacy", endpoint="target_inhibition",
        direction=None, threshold=None,
        citation="McQuisten & Peek, BMC Bioinformatics 2007;8:184",
        note="Antisense activity, not toxicity. Excluded from the toxicity table.",
    ),
}


def _normalize_sequence(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().str.replace("U", "T", regex=False).str.strip()


def normalize() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    raw = pd.read_csv(COMBINED, low_memory=False)

    unknown = sorted(set(raw["source"]) - set(SOURCES))
    if unknown:
        raise SystemExit(
            f"No direction declared for: {unknown}. Add each to SOURCES with its "
            "readout direction and threshold before normalizing."
        )

    tox_rows, eff_rows, per_source = [], [], {}
    for name, spec in SOURCES.items():
        df = raw[raw.source == name].copy()
        if df.empty:
            continue
        df["sequence"] = _normalize_sequence(df["sequence"])
        df["y_raw"] = pd.to_numeric(df["y_raw"], errors="coerce")
        df = df.dropna(subset=["y_raw"])

        declared = set(df["readout"].dropna().unique())
        if declared and spec.readout not in declared:
            raise SystemExit(
                f"{name}: expected readout {spec.readout!r}, found {sorted(declared)}. "
                "The upstream file changed; re-check the direction before trusting it."
            )

        # One row per sequence. Replicates are averaged, never split across a label.
        agg = df.groupby("sequence", as_index=False).agg(
            y_raw=("y_raw", "mean"), n_replicates=("y_raw", "size"),
        )
        agg["source"] = name
        agg["readout"] = spec.readout
        agg["endpoint"] = spec.endpoint
        agg["modality"] = df["modality"].iloc[0]
        agg["assay_type"] = df["assay_type"].iloc[0]

        if spec.kind == "toxicity":
            agg["direction"] = spec.direction
            agg["threshold"] = spec.threshold
            agg["toxic"] = (
                (agg["y_raw"] >= spec.threshold) if spec.direction == HIGHER_IS_TOXIC
                else (agg["y_raw"] < spec.threshold)
            ).astype(int)
            tox_rows.append(agg)
            rate = float(agg["toxic"].mean())
        else:
            eff_rows.append(agg)
            rate = None

        per_source[name] = {
            **asdict(spec),
            "sequences": int(len(agg)),
            "rows_in": int(len(df)),
            "positive_rate": rate,
        }

    tox = pd.concat(tox_rows, ignore_index=True)
    eff = pd.concat(eff_rows, ignore_index=True)
    manifest = {
        "input": str(COMBINED.relative_to(ROOT)),
        "note": "Labels are produced only from a declared direction and threshold. "
                "Efficacy readouts are kept separate and never labelled toxic.",
        "toxicity_sequences": int(len(tox)),
        "efficacy_sequences": int(len(eff)),
        "sources": per_source,
    }
    return tox, eff, manifest


def main() -> None:
    tox, eff, manifest = normalize()
    OUT.mkdir(parents=True, exist_ok=True)
    tox.to_csv(OUT / "external_tox_normalized.csv", index=False)
    eff.to_csv(OUT / "external_efficacy_normalized.csv", index=False)
    (OUT / "external_normalization_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"toxicity  {len(tox):5d} sequences across {tox.source.nunique()} sources")
    for name, sub in tox.groupby("source"):
        spec = SOURCES[name]
        print(f"  {name:26s} {len(sub):5d} seq  {spec.direction:16s} "
              f"thr={spec.threshold:<6} pos={sub.toxic.mean():.3f}")
    print(f"efficacy  {len(eff):5d} sequences across {eff.source.nunique()} sources "
          f"(excluded from toxicity)")


if __name__ == "__main__":
    main()
