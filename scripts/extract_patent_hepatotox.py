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
# How far back to read for a dose statement. An Example's treatment
# paragraph sits a few thousand characters above its first table.
PROSE_WINDOW = 6000
PAUSE = 1.0

OLIGO = re.compile(r"^[ACGT]{15,25}$")
# A gapmer motif such as 5-10-5 or 3-14-3: MOE wing, deoxy gap, MOE wing.
# Ionis routinely files several compounds on one base sequence that differ
# only in this, so two rows sharing a sequence is correct data rather than a
# mapping error, and dropping the column is what makes them look identical.
MOTIF = re.compile(r"^\d{1,2}-\d{1,2}-\d{1,2}$")
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


def _tables_with_spans(html: str) -> list[tuple[list[list[str]], int]]:
    """Every table as (rows, start offset in the html).

    The offset is what lets a dose stated in the Example prose be tied to the
    table it introduces, which is where the dose lives whenever the table
    itself has no dose column.
    """
    out = []
    for match in re.finditer(r"<table[^>]*>.*?</table>", html, re.S | re.I):
        rows = []
        for row in re.findall(r"<tr[^>]*>(.*?)</tr>", match.group(0), re.S | re.I):
            cells = [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", c)).strip()
                     for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S | re.I)]
            if any(cells):
                rows.append(cells)
        out.append((rows, match.start()))
    return out


def _tables(html: str) -> list[list[list[str]]]:
    """Every table as a list of rows of cell text."""
    return [rows for rows, _ in _tables_with_spans(html)]


def _rejoin(rows: list[list[str]]) -> list[list[str]]:
    """Put a wrapped dose back on the label it belongs to.

    Where a table has no dose column the dose sits inside the treatment
    label, and a long label wraps so that "ISIS 233693," and "50 mg/kg"
    arrive as two rows. Without this the dose is lost and three doses of one
    compound look like three unexplained ALT values.
    """
    merged: list[list[str]] = []
    for row in rows:
        if len(row) == 1 and merged and DOSE.search(row[0]):
            merged[-1] = [merged[-1][0] + " " + row[0]] + merged[-1][1:]
        else:
            merged.append(list(row))
    return merged


def _column_names(header_rows: list[list[str]]) -> list[str]:
    """One name per column, from the stacked header rows.

    Google Patents splits a two-line header such as `Dose / (mg/kg)` across
    two rows, so a column is only named by reading down. Single-cell rows are
    the table's caption, not a header, and are left out.
    """
    stacked = [r for r in header_rows if len(r) > 1]
    width = max((len(r) for r in stacked), default=0)
    names = []
    for k in range(width):
        parts = [r[k] for r in stacked if k < len(r) and r[k]]
        names.append(" ".join(parts).strip())
    return names


def _find_column(names: list[str], pattern: str) -> int | None:
    for k, name in enumerate(names):
        if re.search(pattern, name, re.I):
            return k
    return None


FRAGMENT = re.compile(r"^\d{5,7}$")


def _label_dose(label: str) -> float | None:
    """The dose stated inside a treatment label, as in "ISIS 233698, 50 mg/kg"."""
    found = {float(v) for v in DOSE.findall(label)}
    return found.pop() if len(found) == 1 else None


def _groups(data: list[list[str]], label_at: int) -> list[tuple[str, list[list[str]]]]:
    """Split data rows into treatment groups, label carried down.

    A dose-response table prints the compound once and leaves the label blank
    on its lower doses, and a long label wraps so that `ISIS` and `304299`
    land on separate rows. Both mean the label on a row is not the label of
    its group, which is why reading row[0] alone lost every dose-response
    table in the corpus.
    """
    groups: list[tuple[list[str], list[list[str]]]] = []
    for row in data:
        cell = row[label_at] if label_at < len(row) else ""
        if re.search(r"ISIS", cell, re.I) or CONTROL.search(cell):
            groups.append(([cell], [row]))
        elif groups and FRAGMENT.match(cell):
            # The tail of a label that wrapped onto this row.
            groups[-1][0].append(cell)
            groups[-1][1].append(row)
        elif groups and not cell:
            groups[-1][1].append(row)
    return [(" ".join(parts), rows) for parts, rows in groups]


def compound_details(tables: list[list[list[str]]]) -> dict[int, dict]:
    """ISIS number -> its sequence and motif, where a table lists both.

    This is the only place the chemistry is recorded. Two compounds on the
    same base sequence are distinguished by motif alone, so a row without it
    cannot be told apart from its sibling.
    """
    out: dict[int, dict] = {}
    for rows in tables:
        for row in rows:
            seq = next((c.upper().replace("U", "T") for c in row
                        if OLIGO.match(c.upper().replace("U", "T"))), None)
            if not seq:
                continue
            isis = next((int(c) for c in row if c.isdigit() and 100000 <= int(c) < 10000000), None)
            if isis is None:
                continue
            motif = next((c for c in row if MOTIF.match(c.strip())), None)
            out[isis] = {"sequence": seq, "motif": motif}
    return out


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


def tox_rows(tables: list[tuple[list[list[str]], int]]) -> list[dict]:
    """Treatment rows from tables whose header mentions ALT or AST.

    Columns are matched to data cells from the right, because Google Patents
    often emits a leading spacer cell in the data rows that the header does
    not have. Aligning from the left read the wrong column.
    """
    found = []
    for i, (rows, start) in enumerate(tables):
        header = " ".join(" ".join(r) for r in rows)
        if not re.search(r"\bALT\b", header) and not re.search(r"\bAST\b", header):
            continue
        if i + 1 >= len(tables):
            continue
        data = _rejoin([r for r in tables[i + 1][0] if r])
        names = _column_names(rows)
        alt_at = _find_column(names, r"\bALT\b")
        if alt_at is None or not data:
            continue
        ast_at = _find_column(names, r"\bAST\b")
        dose_at = _find_column(names, r"\bdose\b")

        width = max(len(r) for r in data)
        offset = width - len(names)
        if offset < 0:
            continue

        label_at = next((k for k in range(width)
                         if any(re.search(r"ISIS", r[k], re.I) or CONTROL.search(r[k])
                                for r in data if k < len(r))), None)
        if label_at is None:
            continue

        def cell(row, col):
            idx = col + offset
            return row[idx] if 0 <= idx < len(row) else ""

        def number(row, col):
            if col is None:
                return None
            try:
                return float(cell(row, col).replace(",", ""))
            except ValueError:
                return None

        for label, group in _groups(data, label_at):
            for row in group:
                in_row = _label_dose(row[label_at] if label_at < len(row) else "")
                found.append({
                    "label": label,
                    "alt": number(row, alt_at),
                    "ast": number(row, ast_at),
                    "dose": (number(row, dose_at) if number(row, dose_at) is not None
                             else in_row if in_row is not None
                             else _label_dose(label)),
                    "dose_from": ("table" if number(row, dose_at) is not None
                                  else "label" if (in_row is not None
                                                   or _label_dose(label) is not None)
                                  else None),
                    "table": i + 1,
                    "start": start,
                })
    return found


def prose_dose(html: str, before: int, isis: set[int]) -> float | None:
    """The dose from the Example that introduces a table, when unambiguous.

    Most screening tables carry no dose column: every animal got the same
    dose and the prose says so once. That dose is real data and worth having,
    but only where the passage leaves no choice about which value applies, so
    this returns None unless the window states exactly one.
    """
    window = re.sub(r"\s+", " ", re.sub(
        r"<[^>]+>", " ", html[max(0, before - PROSE_WINDOW):before]))
    at = window.rfind("Example")
    if at > 0:
        window = window[at:]
    if not (re.search(r"ISIS\s+oligonucleotide", window, re.I)
            or any(str(n) in window for n in isis)):
        return None
    # Anaesthetic doses are stated in mg/kg too and are not the treatment.
    window = re.sub(r"\b(ketamine|xylazine|warfarin)\b[^.]*", " ", window, flags=re.I)
    window = re.sub(r"\d+(?:\.\d+)?\s*mg/kg\s+(ketamine|xylazine)",
                    " ", window, flags=re.I)
    values = {float(v) for v in DOSE.findall(window)}
    return values.pop() if len(values) == 1 else None


def extract(pub: str, html: str) -> tuple[list[dict], dict[str, int]]:
    text = re.sub(r"[ \t]+", " ", re.sub(r"<[^>]+>", " ", html))
    spans = _tables_with_spans(html)
    tables = [rows for rows, _ in spans]
    seqs = sequence_map(tables)
    details = compound_details(tables)
    isis_to_seq = {int(a): int(b) for a, b in ISIS_SEQID.findall(text)}

    rejected = {"no_compound": 0, "no_seq_id": 0, "no_sequence": 0,
                "no_alt_value": 0, "no_control": 0}
    rows, by_table = [], {}
    for entry in tox_rows(spans):
        by_table.setdefault(entry["table"], []).append(entry)

    for table, entries in by_table.items():
        control = next((e for e in entries
                        if CONTROL.search(e["label"]) and e["alt"] is not None), None)
        numbers = {int(m) for e in entries
                   for m in ISIS_LABEL.findall(e["label"])}
        # One dose for the whole table, used only where the table itself has
        # no dose column to read.
        shared = prose_dose(html, entries[0]["start"], numbers)
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
            # A table that lists the compound directly is better than the
            # prose route, because it carries the motif too.
            detail = details.get(isis)
            if detail:
                seq_id = isis_to_seq.get(isis)
                sequence, motif = detail["sequence"], detail["motif"]
            else:
                if isis not in isis_to_seq:
                    rejected["no_seq_id"] += 1
                    continue
                seq_id = isis_to_seq[isis]
                if seq_id not in seqs:
                    rejected["no_sequence"] += 1
                    continue
                sequence, motif = seqs[seq_id], None
            if control is None:
                # Without a control in the same table an ALT of 52 means
                # nothing: normal mouse ALT spans a wide range between studies.
                rejected["no_control"] += 1
                continue
            dose = entry["dose"] if entry["dose"] is not None else shared
            source = entry["dose_from"] or ("prose" if dose is not None else None)
            rows.append({
                "publication": pub,
                "isis_number": isis,
                "seq_id": seq_id,
                "sequence": sequence,
                "motif": motif,
                "dose_mg_kg": dose,
                "dose_source": source,
                "table": table,
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
    by_source: dict[str, int] = {}
    for row in all_rows:
        key = row["dose_source"] or "unknown"
        by_source[key] = by_source.get(key, 0) + 1
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps({
        "query": args.query,
        "publications_read": len(pubs),
        "rows": len(all_rows),
        "unique_sequences": len(sequences),
        "dose_source": by_source,
        "dose_rule": "A dose is taken from the table's own dose column, else "
                     "from the treatment label, else from the Example prose "
                     "and only where that passage states exactly one dose. "
                     "Otherwise it is left empty, because a liver readout "
                     "without a dose cannot be compared with another study.",
        "rejected": totals,
        "rule": "A row is emitted only when the compound resolves to a real "
                "sequence and its own table carries a control group. A patent "
                "that reports toxicity qualitatively yields nothing.",
        "per_patent": per_patent,
    }, indent=2) + "\n")

    print(f"\n{len(all_rows)} row(s), {len(sequences)} unique sequence(s)")
    print(f"dose: {by_source}")
    if totals:
        print(f"rejected: {totals}")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
