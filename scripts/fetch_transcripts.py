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
GENOMIC_CACHE = ROOT / "data" / "external" / "genomic"
MANIFEST = CACHE / "_manifest.json"
ENSEMBL = "https://rest.ensembl.org"
PAUSE = 0.12
MAX_CDNA = 25_000
MAX_GENOMIC = 400_000
GENOMIC = [False]  # set by --genomic

NETWORK_ERRORS = (urllib.error.HTTPError, urllib.error.URLError,
                  json.JSONDecodeError, TimeoutError, OSError)

# Retries for a request Ensembl refused for a reason that will pass.
MAX_ATTEMPTS = 4
BACKOFF = 2.0


class Unknown(Exception):
    """Ensembl says there is no such thing. Permanent, so stop asking."""


class Transient(Exception):
    """Ensembl declined this time. Says nothing about the gene."""


def _get(path: str, timeout: int = 45):
    """Ensembl JSON, with the two failure kinds kept apart.

    Treating these as one error is how 151 of 350 genes ended up permanently
    blacklisted. A 404 means the symbol is not in Ensembl and asking again is
    pointless. A 429, a 5xx or a timeout means the request lost, and the gene
    is fine: APOC3 and BRCA1 were both recorded as unresolvable and both
    resolve on the first try when asked again. Only the first kind may be
    remembered as a miss.
    """
    request = urllib.request.Request(
        f"{ENSEMBL}{path}",
        headers={"Content-Type": "application/json", "User-Agent": "oligotox"},
    )
    last = "no attempt made"
    for attempt in range(MAX_ATTEMPTS):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read())
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):
                raise Unknown(f"HTTP {exc.code}") from exc
            if exc.code == 429:
                # Ensembl sends Retry-After and means it.
                wait = float(exc.headers.get("Retry-After") or BACKOFF * (attempt + 1))
                time.sleep(min(wait, 30.0))
                continue
            if exc.code >= 500:
                time.sleep(BACKOFF * (attempt + 1))
                continue
            raise Transient(f"HTTP {exc.code}") from exc
        except (urllib.error.URLError, json.JSONDecodeError,
                TimeoutError, OSError) as exc:
            time.sleep(BACKOFF * (attempt + 1))
            last = f"{type(exc).__name__}: {exc}"
    raise Transient(f"gave up after {MAX_ATTEMPTS} attempts ({last})")


def _candidate_symbols(symbol: str) -> list[str]:
    raw = str(symbol).strip()
    seen, out = set(), []
    for candidate in (raw, raw.upper(), raw.replace("-", ""), raw.upper().replace("-", "")):
        if candidate and candidate not in seen:
            seen.add(candidate)
            out.append(candidate)
    return out


def resolve(symbol: str) -> tuple[dict | None, str | None]:
    """Symbol -> sequence record, plus why it failed when it did.

    Returns `(record, None)` on success, `(None, "unknown")` when Ensembl has
    no such symbol, and `(None, "transient")` when the requests failed. Only
    "unknown" may be remembered, since a transient failure says nothing about
    the gene and blacklisting it loses data permanently.
    """
    transient = False
    for candidate in _candidate_symbols(symbol):
        try:
            gene = _get(f"/lookup/symbol/homo_sapiens/{urllib.parse.quote(candidate)}?expand=1")
        except Unknown:
            continue
        except Transient:
            transient = True
            continue
        time.sleep(PAUSE)

        transcripts = gene.get("Transcript") or []
        if not transcripts:
            continue
        canonical = next((x for x in transcripts if x.get("is_canonical")), None)
        chosen = canonical or max(transcripts, key=lambda x: x.get("length") or 0)

        # Gapmers frequently target pre-mRNA, and a cDNA has no introns: for
        # DGAT2, 24% of its ASOs locate in the cDNA and 100% in the genomic
        # span. Genomic is the right target unless cDNA is asked for.
        if GENOMIC[0]:
            query = f"/sequence/id/{gene.get('id')}?type=genomic"
        else:
            query = f"/sequence/id/{chosen['id']}?type=cdna"
        try:
            payload = _get(query, timeout=90)
        except Unknown:
            continue
        except Transient:
            transient = True
            continue
        time.sleep(PAUSE)

        sequence = (payload.get("seq") or "").upper()
        if not sequence:
            continue
        return {
            "symbol": symbol, "resolved_as": candidate, "gene_id": gene.get("id"),
            "transcript_id": chosen["id"], "canonical": bool(canonical),
            "kind": "genomic" if GENOMIC[0] else "cdna",
            "length": len(sequence),
            "sequence": sequence[: (MAX_GENOMIC if GENOMIC[0] else MAX_CDNA)],
            "truncated": len(sequence) > (MAX_GENOMIC if GENOMIC[0] else MAX_CDNA),
        }, None
    return None, ("transient" if transient else "unknown")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--retry-misses", action="store_true")
    parser.add_argument("--genomic", action="store_true",
                        help="fetch the genomic span including introns, which is "
                             "what gapmers target")
    args = parser.parse_args(argv)
    GENOMIC[0] = args.genomic

    if not ATLAS.exists():
        raise SystemExit(f"missing {ATLAS}; run scripts/ingest_aso_atlas.py first")

    import pandas as pd

    genes = (
        pd.read_parquet(ATLAS, columns=["target_gene"])["target_gene"]
        .dropna().astype(str).value_counts()
    )
    symbols = list(genes.index)[: args.limit] if args.limit else list(genes.index)

    # cDNA and genomic are different caches and must not share a manifest,
    # or a gene fetched as cDNA is skipped when genomic is asked for.
    cache_dir = GENOMIC_CACHE if args.genomic else CACHE
    globals()["CACHE"] = cache_dir
    manifest_path = cache_dir / "_manifest.json"
    globals()["MANIFEST"] = manifest_path
    cache_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"hits": {}, "misses": []}
    # A manifest written before transient and unknown were told apart holds
    # both kinds under "misses". Those entries are not trustworthy as
    # permanent, so they are retried once and re-sorted by what happens.
    misses = set(manifest.get("misses", []))
    unknown = set(manifest.get("unknown", []))
    legacy = misses - unknown
    if legacy:
        print(f"{len(legacy)} earlier miss(es) recorded before transient failures "
              "were separated; retrying them")

    fetched = skipped = unknown_now = transient_now = 0
    for i, symbol in enumerate(symbols, 1):
        path = cache_dir / f"{symbol.replace('/', '_')}.json"
        # Only a confirmed unknown is worth skipping. Anything else gets
        # another chance, which is the whole point of the change.
        if path.exists() or (symbol in unknown and not args.retry_misses):
            skipped += 1
            continue

        record, reason = resolve(symbol)
        if record is None:
            if reason == "unknown":
                unknown.add(symbol)
                unknown_now += 1
            else:
                transient_now += 1
            misses.add(symbol)
        else:
            path.write_text(json.dumps(record))
            manifest["hits"][symbol] = {"gene_id": record["gene_id"],
                                        "length": record["length"]}
            misses.discard(symbol)
            unknown.discard(symbol)
            fetched += 1
        if i % 25 == 0:
            print(f"  {i}/{len(symbols)}  fetched={fetched} cached={skipped} "
                  f"unknown={unknown_now} transient={transient_now}")
            manifest["misses"] = sorted(misses)
            manifest["unknown"] = sorted(unknown)
            manifest_path.write_text(json.dumps(manifest, indent=2))

    manifest["misses"] = sorted(misses)
    manifest["unknown"] = sorted(unknown)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nfetched {fetched}, already cached {skipped}, "
          f"not in Ensembl {unknown_now}, failed transiently {transient_now}")
    if transient_now:
        print("Transient failures are not blacklisted. Run again to pick them up.")
    print(f"{len(manifest['hits'])} gene(s) have a sequence in {cache_dir.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
