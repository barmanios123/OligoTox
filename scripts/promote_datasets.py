from pathlib import Path
import csv
import shutil

ROOT = Path.cwd()
REGISTRY = ROOT / "data/registry/dataset_registry.csv"

PROMOTION_RULES = [
    ("archive/_legacy/manifests/", "manifests/"),
    ("archive/_legacy/archive/download_helpers", "archive/download_helpers"),
    ("archive/app_mvp/", "app_mvp_legacy_snapshot/"),
]

def safe_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        if not dst.exists():
            shutil.copytree(src, dst)
    else:
        if not dst.exists():
            shutil.copy2(src, dst)

if not REGISTRY.exists():
    raise FileNotFoundError(f"Missing registry: {REGISTRY}")

rows = []
with REGISTRY.open() as fp:
    rows = list(csv.DictReader(fp))

for row in rows:
    src_rel = row["path"]
    src = ROOT / src_rel
    if not src.exists():
        continue

    promoted = False
    for old_prefix, new_prefix in PROMOTION_RULES:
        if src_rel.startswith(old_prefix):
            suffix = src_rel[len(old_prefix):]
            dst = ROOT / new_prefix / suffix
            safe_copy(src, dst)
            row["promote_to"] = str(dst.relative_to(ROOT))
            row["notes"] = "auto-promoted by rule"
            promoted = True
            break

    if not promoted and row["status"] == "candidate_for_ingestion":
        ext = row["extension"]
        if ext in {".csv", ".tsv", ".xlsx", ".json", ".parquet"}:
            dst = ROOT / "data/processed/unified" / Path(src_rel).name
            safe_copy(src, dst)
            row["promote_to"] = str(dst.relative_to(ROOT))
            row["notes"] = "staged for unified review"

with REGISTRY.open("w", newline="") as fp:
    writer = csv.DictWriter(fp, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("Promotion complete.")
