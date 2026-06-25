from pathlib import Path
import csv
import hashlib

ROOT = Path.cwd()

SEARCH_DIRS = [
    "data",
    "manifests",
    "models",
    "omics",
    "archive",
    "assays",
    "translation",
    "artifacts",
]

EXTS = {".csv", ".tsv", ".xlsx", ".json", ".parquet", ".pkl", ".joblib", ".pdf"}

def sha1_small(path: Path, chunk_size=1024 * 1024):
    h = hashlib.sha1()
    with path.open("rb") as f:
        chunk = f.read(chunk_size)
        h.update(chunk)
    return h.hexdigest()

def infer_role(path: Path) -> str:
    s = str(path).lower()
    if "manifest" in s:
        return "manifest"
    if "/models/" in s or s.startswith("models/"):
        return "model_artifact"
    if "processed" in s:
        return "processed_data"
    if "raw" in s:
        return "raw_data"
    if "archive" in s or "_legacy" in s:
        return "legacy_archive"
    if "omics" in s:
        return "omics"
    if path.suffix.lower() == ".pdf":
        return "reference_pdf"
    return "dataset"

def infer_status(path: Path) -> str:
    s = str(path).lower()
    if "archive/_legacy" in s or "/archive/" in s:
        return "preserved"
    if "data/raw" in s:
        return "candidate_for_ingestion"
    if "manifests/" in s or "data/processed/" in s or "models/" in s:
        return "active"
    return "review"

rows = []
for d in SEARCH_DIRS:
    p = ROOT / d
    if not p.exists():
        continue
    for f in p.rglob("*"):
        if f.is_file() and f.suffix.lower() in EXTS:
            rel = f.relative_to(ROOT)
            rows.append({
                "path": str(rel),
                "filename": f.name,
                "extension": f.suffix.lower(),
                "size_bytes": f.stat().st_size,
                "sha1_head": sha1_small(f),
                "role": infer_role(rel),
                "status": infer_status(rel),
                "source_bucket": str(rel).split("/")[0],
                "promote_to": "",
                "notes": "",
            })

out = ROOT / "data/registry/dataset_registry.csv"
out.parent.mkdir(parents=True, exist_ok=True)

with out.open("w", newline="") as fp:
    writer = csv.DictWriter(fp, fieldnames=list(rows[0].keys()) if rows else [
        "path","filename","extension","size_bytes","sha1_head","role","status","source_bucket","promote_to","notes"
    ])
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda r: r["path"]))

print(f"Wrote {len(rows)} rows to {out}")
