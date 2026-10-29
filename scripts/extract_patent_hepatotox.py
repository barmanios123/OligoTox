"""Pull per-sequence liver toxicity out of antisense patents.

Eight of thirteen endpoints sit at `insufficient_data` on 74 unique
sequences, and the internal corpus is not growing. Patents are where
per-sequence data with real readouts lives: ASO Atlas built 188,521 efficacy
rows from this same corpus. This goes after the toxicity half.

The chain, all three links of which have to resolve before a row is emitted:

    tox table     ISIS 233693, 50 mg/kg | ALT 52 | AST 79
    prose         "ISIS 233693 (SEQ ID NO: 131)"
    sequence table  SEQ ID NO 131 -> the actual 20-mer

Two things about Google Patents' markup that this depends on. It splits every
table into a header table followed immediately by a data table, so a header
row mentioning ALT means the *next* table holds the measurements. And a long
treatment label wraps onto its own single-cell row, so "ISIS 233693," and
"50 mg/kg" arrive separately and have to be rejoined.

**What this will not do.** It emits nothing unless the compound resolves to a
real sequence and the same table carries a control group to compare against.
A patent that says "levels associated with liver toxicity were not observed"
yields no rows, because that is not a measurement. Every rejection is counted
and reported, so the yield is visible rather than implied.

    python scripts/extract_patent_hepatotox.py --limit 20
    python scripts/extract_patent_hepatotox.py --publication US11225664B2

Writes data/external/patent_hepatotox.csv and a report beside it.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "external" / "patent_cache"
OUT_CSV = ROOT / "data" / "external" / "patent_hepatotox.csv"
OUT_JSON = ROOT / "reports" / "results" / "patent_hepatotox.json"

QUERY = 'antisense ALT AST hepatotoxicity gapmer'
SEARCH = "https://patents.google.com/xhr/query"
PATENT = "https://patents.google.com/patent/{pub}/en"
UA = "Mozilla/5.0 (compatible; oligotox research extractor)"
PAUSE = 1.0

OLIGO = re.compile(r"^[ACGT]{15,25}$")
ISIS_SEQID = re.compile(r"ISIS\s+(\d{5,7})\s*\(\s*SEQ\s+ID\s+NO:?\s*(\d+)\s*\)", re.I)
ISIS_LABEL = re.compile(r"ISIS\s+(\d{5,7})", re.I)
DOSE = re.compile(r"(\d+(?:\.\d+)?)\s*mg/kg", re.I)
CONTROL = re.compile(r"\b(saline|pbs|untreated|control|vehicle)\b", re.I)


def _get(url: str, timeout: int = 60) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", "replace")


def search(query: str, limit: int) -> list[str]:
    """Candidate publication numbers, US only.

    US publications keep their tables as markup. Other offices are often
    OCR'd text where the data rows are simply absent, which is why the first
    patent looked at produced nothing extractable.
    """
    found, page = [], 0
    while len(found) < limit and page < 5:
        url = (f"{SEARCH}?" + urllib.parse.urlencode(
            {"url": f"q={query}&country=US&num=100&page={page}"}))
        try:
            payload = json.loads(_get(url))
        except (urllib.error.URLError, json.JSONDecodeError, TimeoutError):
            break
        clusters = payload.get("results", {}).get("cluster", [])
        results = [r for c in clusters for r in c.get("result", [])]
        if not results:
            break
        for record in results:
            pub = record.get("patent", {}).get("publication_number")
            if pub and pub not in found:
                found.append(pub)
        page += 1
        time.sleep(PAUSE)
    return found[:limit]


def fetch(pub: str) -> str | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{pub}.html"
    if path.exists():
        return path.read_text(encoding="utf-8", errors="replace")
    try:
        html = _get(PATENT.format(pub=pub))
    except (urllib.error.URLError, TimeoutError, OSError):
        return None
    path.write_text(html, encoding="utf-8")
    time.sleep(PAUSE)
    return html


def _tables(html: str) -> list[list[list[str]]]:
    """Every table as a list of rows of cell text."""
    out = []
    for block in re.findall(r"<table[^>]*>.*?</table>", html, re.S | re.I):
        rows = []
        for row in re.findall(r"<tr[^>]*>(.*?)</tr>", block, re.S | re.I):
            cells = [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", c)).strip()
                     for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S | re.I)]
            if any(cells):
                rows.append(cells)
        out.append(rows)
    return out


def _rejoin(rows: list[list[str]]) -> list[list[str]]:
    """Put wrapped treatment labels back on their own row.

    "ISIS 233693," and "50 mg/kg" arrive as two rows because the label wraps.
    Without this the dose is lost and two different doses of one compound look
    like duplicates.
    """
    merged: list[list[str]] = []
    for row in rows:
        if len(row) == 1 and merged and DOSE.search(row[0]):
            merged[-1] = [merged[-1][0] + " " + row[0]] + merged[-1][1:]
        else:
            merged.append(list(row))
    return merged


def sequence_map(tables: list[list[list[str]]]) -> dict[int, str]:
    """SEQ ID NO -> sequence, from the sequence listing tables."""
    out: dict[int, str] = {}
    for i, rows in enumerate(tables):
        header = " ".join(" ".join(r) for r in rows).lower()
        if "seq" not in header or "sequence" not in header:
            continue
        for data in tables[i:i + 2]:
            for row in data:
                seq = next((c.upper().replace("U", "T") for c in row
                            if OLIGO.match(c.upper().replace("U", "T"))), None)
                if not seq:
                    continue
                ids = [int(c) for c in row if c.isdigit() and 0 < int(c) < 100000]
                if ids:
                    # The SEQ ID NO is the last numeric cell on the row.
                    out[ids[-1]] = seq
    return out


def tox_rows(tables: list[list[list[str]]]) -> list[dict]:
    """Treatment rows from tables whose header mentions ALT or AST."""
    found = []
    for i, rows in enumerate(tables):
        header = " ".join(" ".join(r) for r in rows)
        if not re.search(r"\bALT\b", header) and not re.search(r"\bAST\b", header):
            continue
        columns = [c.strip().upper() for r in rows for c in r]
        if i + 1 >= len(tables):
            continue
        data = _rejoin(tables[i + 1])
        if not data:
            continue
        try:
            alt_at = columns.index("ALT")
            ast_at = columns.index("AST") if "AST" in columns else None
        except ValueError:
            continue
        # The label occupies column 0, so a header column k is data cell k.
        for row in data:
            label = row[0] if row else ""
            values = row[1:]
            def value(pos):
                idx = pos - 1
                if idx is None or idx < 0 or idx >= len(values):
                    return None
                try:
                    return float(values[idx])
                except ValueError:
                    return None
            found.append({
                "label": label,
                "alt": value(alt_at),
                "ast": value(ast_at) if ast_at is not None else None,
                "table": i + 1,
            })
    return found


def extract(pub: str, html: str) -> tuple[list[dict], dict[str, int]]:
    text = re.sub(r"[ \t]+", " ", re.sub(r"<[^>]+>", " ", html))
    tables = _tables(html)
    seqs = sequence_map(tables)
    isis_to_seq = {int(a): int(b) for a, b in ISIS_SEQID.findall(text)}

    rejected = {"no_compound": 0, "no_seq_id": 0, "no_sequence": 0,
                "no_alt_value": 0, "no_control": 0}
    rows, by_table = [], {}
    for entry in tox_rows(tables):
        by_table.setdefault(entry["table"], []).append(entry)

    for table, entries in by_table.items():
        control = next((e for e in entries
                        if CONTROL.search(e["label"]) and e["alt"] is not None), None)
        for entry in entries:
            if CONTROL.search(entry["label"]):
                continue
            match = ISIS_LABEL.search(entry["label"])
            if not match:
                rejected["no_compound"] += 1
                continue
            isis = int(match.group(1))
            if entry["alt"] is None:
                rejected["no_alt_value"] += 1
                continue
            if isis not in isis_to_seq:
                rejected["no_seq_id"] += 1
                continue
            seq_id = isis_to_seq[isis]
            if seq_id not in seqs:
                rejected["no_sequence"] += 1
                continue
            if control is None:
                # Without a control in the same table an ALT of 52 means
                # nothing: normal mouse ALT spans a wide range between studies.
                rejected["no_control"] += 1
                continue
            dose = DOSE.search(entry["label"])
            rows.append({
                "publication": pub,
                "isis_number": isis,
                "seq_id": seq_id,
                "sequence": seqs[seq_id],
                "dose_mg_kg": float(dose.group(1)) if dose else None,
                "alt": entry["alt"],
                "ast": entry["ast"],
                "control_alt": control["alt"],
                "alt_fold_vs_control": round(entry["alt"] / control["alt"], 3)
                if control["alt"] else None,
            })
    return rows, rejected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--limit", type=int, default=25)
    parser.add_argument("--publication", action="append", default=None)
    parser.add_argument("--query", default=QUERY)
    args = parser.parse_args(argv)

    pubs = args.publication or search(args.query, args.limit)
    print(f"{len(pubs)} publication(s) to read\n")

    all_rows, totals, per_patent = [], {}, []
    for i, pub in enumerate(pubs, 1):
        html = fetch(pub)
        if html is None:
            print(f"  {pub:16s} could not fetch")
            continue
        rows, rejected = extract(pub, html)
        all_rows.extend(rows)
        for key, count in rejected.items():
            totals[key] = totals.get(key, 0) + count
        per_patent.append({"publication": pub, "rows": len(rows), "rejected": rejected})
        print(f"  {pub:16s} {len(rows):3d} row(s)"
              + (f"   rejected: {rejected}" if any(rejected.values()) else ""))

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    if all_rows:
        import csv
        with OUT_CSV.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
            writer.writeheader()
            writer.writerows(all_rows)

    sequences = {r["sequence"] for r in all_rows}
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps({
        "query": args.query,
        "publications_read": len(pubs),
        "rows": len(all_rows),
        "unique_sequences": len(sequences),
        "rejected": totals,
        "rule": "A row is emitted only when the compound resolves to a real "
                "sequence and its own table carries a control group. A patent "
                "that reports toxicity qualitatively yields nothing.",
        "per_patent": per_patent,
    }, indent=2) + "\n")

    print(f"\n{len(all_rows)} row(s), {len(sequences)} unique sequence(s)")
    if totals:
        print(f"rejected: {totals}")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
