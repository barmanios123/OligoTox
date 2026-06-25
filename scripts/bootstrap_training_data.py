"""
Bootstrap a minimal training_dataset.csv from demo-like candidates.

This is a placeholder so train.py can run end-to-end.
Replace rows with real curated data as you ingest Dieckmann,
Papargyri, GalNAc siRNA, etc.
"""

from pathlib import Path

import pandas as pd

from oligotox.data_models import (
    OligoInput,
    DeliveryInput,
    Modality,
    DeliveryVehicle,
    PlateletComplement,
    OmicsReadout,
    OligoToxInput,
)
from oligotox.features.engineer import extract_features

OUT_PATH = Path("data/processed/training_dataset.csv")


def build_row(
    candidate_id: str,
    oligo: OligoInput,
    delivery: DeliveryInput,
    platelet: PlateletComplement | None = None,
    omics: OmicsReadout | None = None,
    # Placeholder labels – set manually after bootstrapping
    labels: dict[str, int] | None = None,
) -> dict:
    if platelet is None:
        platelet = PlateletComplement()
    if omics is None:
        omics = OmicsReadout()
    if labels is None:
        labels = {}

    inp = OligoToxInput(
        oligo=oligo,
        delivery=delivery,
        platelet_complement=platelet,
        omics=omics,
    )

    feats = extract_features(inp)
    row = {"candidate_id": candidate_id}
    row.update(feats)
    row.update(labels)
    return row


def main() -> None:
    rows: list[dict] = []

    # Example: LNP siRNA (low–moderate tox placeholder labels)
    oligo_lnp = OligoInput(
        sequence="AUGCUACGUACGUACGUACG",
        modality=Modality.SIRNA_LNP,
        target_gene="TTR",
    )
    delivery_lnp = DeliveryInput(
        vehicle=DeliveryVehicle.LNP,
        route="intravenous",
        dose_mg_kg=0.3,
        repeat_dose=True,
    )
    labels_lnp = {
        "hepatocellular_injury": 0,
        "hepatic_vacuolation": 0,
        "renal_tubular_injury": 0,
        "renal_accumulation": 0,
        "thrombocytopenia": 0,
        "coagulation_effect": 0,
        "complement_activation": 1,
        "innate_immune_activation": 1,
        "cytokine_release": 0,
        "kupffer_cell_activation": 0,
        "cns_toxicity": 0,
        "offtarget_gene_silencing": 0,
        "injection_site_reaction": 0,
    }
    rows.append(
        build_row("DEMO_LNP_siRNA", oligo_lnp, delivery_lnp, labels=labels_lnp)
    )

    # Add more demo-style rows (GalNAc siRNA, PS-ASO gapmer, IT lncRNA) similarly
    # with placeholder labels you will overwrite later.

    df = pd.DataFrame(rows)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(df)} rows to {OUT_PATH}")


if __name__ == "__main__":
    main()