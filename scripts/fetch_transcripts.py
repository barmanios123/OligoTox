"""Fetch canonical transcripts for the genes ASO Atlas targets.

Two features need real target sequence, and ASO Atlas gives only a gene
symbol: target-site accessibility, and the off-target complementarity burden
that is the published mechanism behind gapmer hepatotoxicity.

Asks Ensembl for the canonical transcript specifically. Requesting every cDNA
for a gene returns dozens of isoforms (APP alone has 63), which is slow
enough to time out and no more informative for locating a binding site.

Caching is per gene and resumable: a rerun fetches only what is missing, and
a gene Ensembl cannot place is recorded as a miss rather than retried
forever. Symbols in the corpus are untidy (`K-Ras`, `A-431`), so a few
normalisations are tried before giving up.

    python scripts/fetch_transcripts.py              # every gene in the corpus
    python scripts/fetch_transcripts.py --limit 50   # a sample

Writes data/external/transcripts/<SYMBOL>.json and a manifest. Ensembl asks
for no more than ~15 requests a second; this stays well under.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT / "data" / "external" / "aso_atlas.parquet"
CACHE = ROOT / "data" / "external" / "transcripts"
MANIFEST = CACHE / "_manifest.json"
ENSEMBL = "https://rest.ensembl.org"
PAUSE = 0.12
MAX_CDNA = 25_000

NETWORK_ERRORS = (urllib.error.HTTPError, urllib.error.URLError,
                  json.JSONDecodeError, TimeoutError, OSError)


def _get(path: str, timeout: int = 45):
    request = urllib.request.Request(
        f"{ENSEMBL}{path}",
        headers={"Content-Type": "application/json", "User-Agent": "oligotox"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read())


def _candidate_symbols(symbol: str) -> list[str]:
    raw = str(symbol).strip()
    seen, out = set(), []
    for candidate in (raw, raw.upper(), raw.replace("-", ""), raw.upper().replace("-", "")):
        if candidate and candidate not in seen:
            seen.add(candidate)
            out.append(candidate)
    return out


def resolve(symbol: str) -> dict | None:
    """Symbol -> canonical transcript cDNA, or None if Ensembl cannot place it."""
    for candidate in _candidate_symbols(symbol):
        try:
            gene = _get(f"/lookup/symbol/homo_sapiens/{urllib.parse.quote(candidate)}?expand=1")
        except NETWORK_ERRORS:
            continue
        time.sleep(PAUSE)

        transcripts = gene.get("Transcript") or []
        if not transcripts:
            continue
        canonical = next((x for x in transcripts if x.get("is_canonical")), None)
        chosen = canonical or max(transcripts, key=lambda x: x.get("length") or 0)

        try:
            payload = _get(f"/sequence/id/{chosen['id']}?type=cdna")
        except NETWORK_ERRORS:
            continue
        time.sleep(PAUSE)

        sequence = (payload.get("seq") or "").upper()
        if not sequence:
            continue
        return {
            "symbol": symbol, "resolved_as": candidate, "gene_id": gene.get("id"),
            "transcript_id": chosen["id"], "canonical": bool(canonical),
            "length": len(sequence), "sequence": sequence[:MAX_CDNA],
            "truncated": len(sequence) > MAX_CDNA,
        }
    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--retry-misses", action="store_true")
    args = parser.parse_args(argv)

    if not ATLAS.exists():
        raise SystemExit(f"missing {ATLAS}; run scripts/ingest_aso_atlas.py first")

    import pandas as pd

    genes = (
        pd.read_parquet(ATLAS, columns=["target_gene"])["target_gene"]
        .dropna().astype(str).value_counts()
    )
    symbols = list(genes.index)[: args.limit] if args.limit else list(genes.index)

    CACHE.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {"hits": {}, "misses": []}
    misses = set(manifest.get("misses", []))

    fetched = skipped = failed = 0
    for i, symbol in enumerate(symbols, 1):
        path = CACHE / f"{symbol.replace('/', '_')}.json"
        if path.exists() or (symbol in misses and not args.retry_misses):
            skipped += 1
            continue

        record = resolve(symbol)
        if record is None:
            misses.add(symbol)
            failed += 1
        else:
            path.write_text(json.dumps(record))
            manifest["hits"][symbol] = {"gene_id": record["gene_id"],
                                        "length": record["length"]}
            fetched += 1
        if i % 25 == 0:
            print(f"  {i}/{len(symbols)}  fetched={fetched} cached={skipped} missed={failed}")
            manifest["misses"] = sorted(misses)
            MANIFEST.write_text(json.dumps(manifest, indent=2))

    manifest["misses"] = sorted(misses)
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nfetched {fetched}, already cached {skipped}, unresolved {failed}")
    print(f"{len(manifest['hits'])} gene(s) have a transcript in {CACHE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
