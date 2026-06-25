from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_CANDIDATES = [
    ROOT / "data/processed/training/training_manifest.csv",
    ROOT / "data/processed/training/oligotox_training_manifest.csv",
    ROOT / "data/processed/manifests/training_manifest.csv",
    ROOT / "data/processed/manifests/oligotox_training_manifest.csv",
    ROOT / "data/interim/training_manifest.csv",
    ROOT / "data/interim/oligotox_training_manifest.csv",
    ROOT / "data/training_manifest.csv",
]
OUT_TABLE = ROOT / "data/processed/unified/oligotox_unified_training_table.csv"
OUT_LOG = ROOT / "data/processed/unified/oligotox_unified_training_log.json"


def find_manifest() -> Path:
    for p in MANIFEST_CANDIDATES:
        if p.exists():
            return p
    hits = list(ROOT.rglob("*training*manifest*.csv"))
    if hits:
        return hits[0]
    raise FileNotFoundError("Could not find a training manifest CSV")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = []
    for c in df.columns:
        c = str(c).strip().lower()
        c = c.replace("\n", " ")
        c = "_".join(c.split())
        cols.append(c)
    df.columns = cols
    return df


def make_columns_unique(df: pd.DataFrame) -> pd.DataFrame:
    seen = {}
    out = []
    for c in map(str, df.columns):
        n = seen.get(c, 0)
        out.append(c if n == 0 else f"{c}__dup{n}")
        seen[c] = n + 1
    df.columns = out
    return df


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", low_memory=False)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported file type: {path}")


def resolve_source_path(rel: str) -> Path:
    rel = str(rel).strip()
    p = Path(rel)
    if p.is_absolute():
        return p
    candidates = [
        ROOT / rel,
        ROOT / "data" / rel,
        ROOT / "archive" / rel,
    ]
    for c in candidates:
        if c.exists():
            return c
    hits = list(ROOT.rglob(Path(rel).name))
    if hits:
        return hits[0]
    raise FileNotFoundError(f"Could not resolve source path: {rel}")


def load_manifest(manifest_path: Path) -> pd.DataFrame:
    df = pd.read_csv(manifest_path, low_memory=False)
    df.columns = [str(c).strip().lower() for c in df.columns]

    path_col = None
    for candidate in ["path", "source_path", "file_path", "filepath", "file", "dataset_path"]:
        if candidate in df.columns:
            path_col = candidate
            break
    if path_col is None:
        raise ValueError(f"Manifest missing path column. Found columns: {df.columns.tolist()}")

    if "status" not in df.columns:
        df["status"] = ""
    if "role" not in df.columns:
        df["role"] = ""
    if "source_bucket" not in df.columns:
        df["source_bucket"] = ""

    df = df.rename(columns={path_col: "path"})
    return df


def unify() -> None:
    manifest_path = find_manifest()
    manifest = load_manifest(manifest_path)

    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    if OUT_TABLE.exists():
        OUT_TABLE.unlink()

    meta = []
    first_write = True
    seen_rows = set()
    total_rows_written = 0

    print(f"Using manifest: {manifest_path}")

    for _, row in manifest.iterrows():
        rel = str(row.get("path", "")).strip()
        if not rel:
            meta.append({"path": rel, "status": "skipped", "reason": "empty path"})
            continue

        try:
            src = resolve_source_path(rel)
            df = read_table(src)
            df = normalize_columns(df)
            df = make_columns_unique(df)
            df = df.loc[:, ~df.columns.duplicated()]
            df = df.reset_index(drop=True)

            df["__source_path"] = rel
            df["__source_role"] = row.get("role", "")
            df["__source_status"] = row.get("status", "")
            df["__source_bucket"] = row.get("source_bucket", "")

            df = df.astype("string")

            row_keys = df.fillna("").astype(str).agg("||".join, axis=1)
            keep_mask = ~row_keys.isin(seen_rows)
            kept = df.loc[keep_mask].copy()

            for k in row_keys[keep_mask].tolist():
                seen_rows.add(k)

            if not kept.empty:
                kept.to_csv(OUT_TABLE, mode="w" if first_write else "a", index=False, header=first_write)
                first_write = False
                total_rows_written += len(kept)

            meta.append({
                "path": rel,
                "resolved_path": str(src),
                "status": "included",
                "n_rows_read": int(df.shape[0]),
                "n_rows_written": int(kept.shape[0]),
                "n_cols": int(df.shape[1]),
                "role": row.get("role", ""),
                "status_tag": row.get("status", ""),
                "source_bucket": row.get("source_bucket", ""),
            })
            print(f"[OK] {rel} | read={len(df)} | wrote={len(kept)} | cols={df.shape[1]}")

        except Exception as e:
            meta.append({
                "path": rel,
                "status": "error",
                "reason": str(e),
            })
            print(f"[ERR] {rel} | {e}")

    if first_write:
        raise RuntimeError("No tables included in unified dataset")

    OUT_LOG.write_text(json.dumps(meta, indent=2))
    print(f"\nWrote unified training table: {OUT_TABLE} (rows={total_rows_written})")
    print(f"Wrote log: {OUT_LOG}")

    included = [m for m in meta if m.get("status") == "included"]
    errored = [m for m in meta if m.get("status") == "error"]

    print("\n===== Inclusion summary =====")
    print(f"Included sources: {len(included)}")
    print(f"Errored sources: {len(errored)}")
    for m in included:
        print(f"[OK] {m['path']} | read={m['n_rows_read']} | wrote={m['n_rows_written']} | cols={m['n_cols']}")
    for m in errored:
        print(f"[ERR] {m['path']} | reason={m['reason']}")


if __name__ == "__main__":
    unify()
