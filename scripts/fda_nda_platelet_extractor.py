from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


SEEDED_NDA_DATA = {
    "inotersen": [],
    "mipomersen": [],
    "nusinersen": [],
    "eplontersen": [],
}


def label_thrombocytopenia(platelet_nadir, baseline_or_control, reviewer_flag=False):
    if platelet_nadir is None or baseline_or_control in (None, 0):
        return None
    ratio = platelet_nadir / baseline_or_control
    if reviewer_flag:
        return 1
    if platelet_nadir < 50 or ratio <= 0.50:
        return 1
    if platelet_nadir < 75 or ratio <= 0.75:
        return 1
    return 0


def thrombocytopenia_grade(platelet_nadir, baseline_or_control):
    if platelet_nadir is None or baseline_or_control in (None, 0):
        return None
    ratio = platelet_nadir / baseline_or_control
    if platelet_nadir < 50 or ratio <= 0.50:
        return 3
    if platelet_nadir < 75 or ratio <= 0.75:
        return 2
    if platelet_nadir < 100 or ratio <= 0.90:
        return 1
    return 0


def label_coagulation(pt_ratio=None, aptt_ratio=None, reviewer_flag=False):
    if reviewer_flag:
        return 1
    vals = [v for v in [pt_ratio, aptt_ratio] if v is not None]
    if not vals:
        return None
    if any(v >= 1.5 for v in vals):
        return 1
    return 0


def label_complement(fold_change=None, reviewer_flag=False):
    if reviewer_flag:
        return 1
    if fold_change is None:
        return None
    if fold_change >= 2.0:
        return 1
    return 0


def build_rows(drug: str):
    rows = []
    for i, arm in enumerate(SEEDED_NDA_DATA.get(drug, []), start=1):
        platelet_baseline = arm.get("platelet_baseline_or_control")
        platelet_nadir = arm.get("platelet_nadir")
        pt_ratio = arm.get("pt_ratio")
        aptt_ratio = arm.get("aptt_ratio")
        complement_fold = arm.get("complement_fold_change_max")

        row = {
            "candidate_id": arm.get("candidate_id", f"{drug}_nda_arm_{i}"),
            "compound": drug,
            "source_type": "fda_nda_pharmtox",
            "source_file": arm.get("source_file", ""),
            "study_id": arm.get("study_id", ""),
            "species": arm.get("species", ""),
            "route": arm.get("route", ""),
            "dose_mg_kg": arm.get("dose_mg_kg"),
            "duration_days": arm.get("duration_days"),
            "platelet_baseline_or_control": platelet_baseline,
            "platelet_nadir": platelet_nadir,
            "platelet_nadir_day": arm.get("platelet_nadir_day"),
            "pt_ratio": pt_ratio,
            "aptt_ratio": aptt_ratio,
            "complement_marker": arm.get("complement_marker", ""),
            "complement_fold_change_max": complement_fold,
            "reviewer_flag_thrombocytopenia": int(bool(arm.get("reviewer_flag_thrombocytopenia", False))),
            "reviewer_flag_coagulopathy": int(bool(arm.get("reviewer_flag_coagulopathy", False))),
            "reviewer_flag_complement": int(bool(arm.get("reviewer_flag_complement", False))),
            "thrombocytopenia": label_thrombocytopenia(
                platelet_nadir,
                platelet_baseline,
                reviewer_flag=bool(arm.get("reviewer_flag_thrombocytopenia", False)),
            ),
            "thrombocytopenia_grade": thrombocytopenia_grade(
                platelet_nadir,
                platelet_baseline,
            ),
            "coagulation_effect": label_coagulation(
                pt_ratio=pt_ratio,
                aptt_ratio=aptt_ratio,
                reviewer_flag=bool(arm.get("reviewer_flag_coagulopathy", False)),
            ),
            "complement_activation": label_complement(
                fold_change=complement_fold,
                reviewer_flag=bool(arm.get("reviewer_flag_complement", False)),
            ),
        }
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--drug", choices=["inotersen", "mipomersen", "nusinersen", "eplontersen"])
    parser.add_argument("--pdf", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    drugs = ["inotersen", "mipomersen", "nusinersen", "eplontersen"] if args.all else [args.drug]
    rows = []
    for drug in drugs:
        if drug is None:
            continue
        rows.extend(build_rows(drug))

    df = pd.DataFrame(rows)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} rows)")


if __name__ == "__main__":
    main()
