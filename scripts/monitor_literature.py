"""Watch the literature for oligonucleotide toxicity datasets worth ingesting.

The continuous pipeline recomputes the same four sources forever unless new
data arrives. This is the inbound stream: a weekly query against Europe PMC
for papers that plausibly carry a usable dataset, scored and deduplicated so
only genuinely new candidates are reported.

It does not ingest anything. A paper is a lead, and turning one into a source
means reading it, establishing the readout direction and threshold, and
declaring it in `data/registry/sources.json`. That judgement is the thing the
registry exists to keep out of automation.

Scoring favours what makes a paper ingestible rather than what makes it
interesting:

    supplementary data     a dataset usually lives in a supplement
    open access            needed to read the supplement at all
    panel language         "screen", "library", "panel", "N gapmers" imply
                           many oligos rather than one drug
    a named tox readout    ALT, AST, caspase, FOB, calcium oscillation
    recency                anything already ingested is filtered out

State lives in reports/results/literature_seen.json, so each run reports only
what it has not reported before.

    python scripts/monitor_literature.py
    python scripts/monitor_literature.py --since 2025-01-01 --min-score 4
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "reports" / "results"
SEEN = RESULTS / "literature_seen.json"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

QUERIES = [
    '(antisense OR gapmer OR "ASO") AND (hepatotoxicity OR "liver injury" OR ALT)',
    '(siRNA OR "small interfering RNA") AND (toxicity OR immunostimulation OR "cytokine release")',
    '(oligonucleotide) AND (neurotoxicity OR "calcium oscillation" OR "functional observational")',
    '(gapmer OR "antisense oligonucleotide") AND (screen OR library OR panel) AND sequence',
]

# Signals that a paper carries a dataset rather than a single compound.
PANEL = re.compile(
    r"\b(screen(ing|ed)?|librar(y|ies)|panel|set of \d+|\d{2,}\s+(gapmers?|oligo\w*|siRNAs?|ASOs?))\b",
    re.I,
)
READOUT = re.compile(
    r"\b(ALT|AST|aminotransferase|caspase|FOB|functional observational|"
    r"calcium oscillation|TLR\d|cytokine release|thrombocytopeni\w+)\b",
    re.I,
)
ALREADY_HAVE = re.compile(
    r"\b(hagedorn|papargyri|mcquisten|ichihara|oligogym)\b", re.I
)


def _search(query: str, since: str, page_size: int) -> list[dict]:
    full = f'{query} AND (FIRST_PDATE:[{since} TO {date.today().isoformat()}])'
    url = f"{API}?{urllib.parse.urlencode({'query': full, 'format': 'json', 'pageSize': page_size, 'resultType': 'core'})}"
    try:
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = json.loads(response.read())
    except (urllib.error.URLError, json.JSONDecodeError) as exc:
        print(f"  query failed ({type(exc).__name__}), skipping: {query[:50]}")
        return []
    return payload.get("resultList", {}).get("result", [])


def _score(record: dict) -> tuple[int, list[str]]:
    text = " ".join(filter(None, [record.get("title"), record.get("abstractText")]))
    score, why = 0, []

    if record.get("hasSuppl") == "Y":
        score += 3
        why.append("has supplementary data")
    if record.get("isOpenAccess") == "Y":
        score += 2
        why.append("open access")
    if PANEL.search(text):
        score += 3
        why.append("describes a panel or screen, not one compound")
    readouts = sorted({m.group(0).upper() for m in READOUT.finditer(text)})
    if readouts:
        score += 2
        why.append(f"names a readout ({', '.join(readouts[:3])})")
    if ALREADY_HAVE.search(text):
        score -= 4
        why.append("looks like a source already ingested")
    return score, why


def _load_seen() -> set[str]:
    if not SEEN.exists():
        return set()
    return set(json.loads(SEEN.read_text()).get("ids", []))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--since", default="2024-01-01")
    parser.add_argument("--min-score", type=int, default=5)
    parser.add_argument("--page-size", type=int, default=25)
    parser.add_argument("--all", action="store_true",
                        help="include candidates reported on a previous run")
    args = parser.parse_args(argv)

    seen = set() if args.all else _load_seen()
    candidates: dict[str, dict] = {}

    for query in QUERIES:
        print(f"querying: {query[:62]}...")
        for record in _search(query, args.since, args.page_size):
            ident = record.get("pmid") or record.get("id")
            if not ident or ident in seen or ident in candidates:
                continue
            score, why = _score(record)
            if score < args.min_score:
                continue
            candidates[ident] = {
                "id": ident,
                "score": score,
                "year": record.get("pubYear"),
                "title": (record.get("title") or "").strip(),
                "journal": (record.get("journalTitle") or "").strip(),
                "doi": record.get("doi"),
                "open_access": record.get("isOpenAccess") == "Y",
                "has_supplementary": record.get("hasSuppl") == "Y",
                "why": why,
                "url": f"https://europepmc.org/article/MED/{ident}",
            }

    # Standing list: a quiet week must not erase leads found earlier.
    existing = {}
    candidates_path = RESULTS / "literature_candidates.json"
    if candidates_path.exists():
        existing = {c["id"]: c for c in
                    json.loads(candidates_path.read_text()).get("candidates", [])}
    new_ids = sorted(set(candidates) - set(existing))
    merged = {**existing, **candidates}

    ranked = sorted(merged.values(), key=lambda c: -c["score"])
    RESULTS.mkdir(parents=True, exist_ok=True)
    candidates_path.write_text(json.dumps({
        "checked": datetime.now(timezone.utc).isoformat(),
        "since": args.since,
        "min_score": args.min_score,
        "note": "Leads, not sources. A candidate becomes ingestible only once its "
                "readout direction and label threshold are declared in "
                "data/registry/sources.json.",
        "new_this_run": new_ids,
        "candidates": ranked,
    }, indent=2) + "\n")

    if not args.all:
        SEEN.write_text(json.dumps({
            "updated": datetime.now(timezone.utc).isoformat(),
            "ids": sorted(seen | set(candidates)),
        }, indent=2) + "\n")

    fresh = [c for c in ranked if c["id"] in new_ids]
    print(f"\n{len(fresh)} new candidate(s) this run, {len(ranked)} standing, "
          f"scoring >= {args.min_score}\n")
    for candidate in (fresh or ranked)[:12]:
        print(f"  [{candidate['score']}] {candidate['year']} {candidate['title'][:76]}")
        print(f"      {'; '.join(candidate['why'])}")
        print(f"      {candidate['url']}")
    if not fresh:
        print("  nothing new this run; the standing list above is unchanged.")
    print(f"\nwrote {(RESULTS / 'literature_candidates.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
