"""Decide what each literature candidate is, and record the decision.

The monitor finds papers; without this they accumulate unread. 34 of the
first 36 were still unassessed, which makes the monitor a backlog generator
rather than an inbound stream.

Each candidate is fetched once and classified from its abstract:

    dataset        a screen or panel with a measured readout, the only class
                   worth a human reading it
    single_study   one compound or construct characterised in depth; real
                   science, but it cannot supply training rows
    method         a platform, chemistry or delivery paper
    reference      a database or atlas that records events per drug rather
                   than outcomes per sequence
    review         a survey of other people's work; it has no rows of its own
    decided        already in data/registry/sources.json, adopted or rejected

`dataset` stays open, and so does `unknown`: where only a title came back
there is not enough to dismiss a paper on. Everything else is marked assessed
with its reason, so the monitor reports it once and never again.

The classifier is deliberately blunt and errs toward `dataset`: a false
positive costs one person five minutes, a false negative loses a dataset.
Its decision is recorded alongside the abstract it read, so a wrong call is
visible rather than buried.

    python scripts/triage_literature.py
    python scripts/triage_literature.py --recheck   # re-assess everything

Writes reports/results/literature_triage.json.
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
RESULTS = ROOT / "reports" / "results"
CANDIDATES = RESULTS / "literature_candidates.json"
TRIAGE = RESULTS / "literature_triage.json"
# Hand-written calls that outrank the classifier, so a paper someone actually
# read is never re-flagged by a regex.
DECISIONS = ROOT / "data" / "registry" / "literature_decisions.json"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
PAUSE = 0.2

# A screen reports many oligos with a measured outcome.
PANEL = re.compile(
    r"\b(screen(ed|ing)?|librar(y|ies)|panel of|a set of \d+|\d{2,}\s+"
    r"(gapmers?|oligonucleotides?|ASOs?|siRNAs?|sequences?)|systematic(ally)? (tested|evaluated))\b",
    re.I,
)
READOUT = re.compile(
    r"\b(ALT|AST|aminotransferase|caspase|FOB score|calcium oscillation|"
    r"knockdown|IC50|EC50|percent inhibition|cytotoxicity|viability)\b", re.I,
)
# One compound, characterised in depth.
SINGLE = re.compile(
    r"\b(we (report|describe|present) (a|an|the) (novel|new)?\s*\w+\s*(ASO|oligonucleotide|"
    r"compound)|first[- ]in[- ]class|case report|phase [I1V]+ (trial|study)|"
    r"targeting \w+ (ameliorates|improves|rescues|reduces))\b", re.I,
)
METHOD = re.compile(
    r"\b(platform|delivery system|conjugat|nanoparticle|formulation|"
    r"chemistry|modification strateg|synthes[ie]s|scaffold|bottlebrush)\b", re.I,
)
REFERENCE = re.compile(r"\b(atlas|database|registry|adverse[- ]event|pharmacovigilance)\b", re.I)
# A review uses the same panel-and-readout vocabulary as a screen while
# reporting no measurements of its own, which is how a CNS safety review was
# first misread as a dataset.
REVIEW = re.compile(
    r"\b(this review (summarizes|summarises|describes|covers)|we review|"
    r"review (summarizes|summarises)|current (approaches|strategies) for|"
    r"recent (advances|progress) in|here,? we discuss)\b", re.I,
)


def decided_pmids() -> dict[str, str]:
    """PMIDs already recorded in the source registry, with their status.

    The registry is the record of what has been adopted or rejected. Without
    this the triage re-flags a paper whose dataset is already shipping, which
    is how the ASO Atlas paper came back up as unassessed.
    """
    registry = ROOT / "data" / "registry" / "sources.json"
    if not registry.exists():
        return {}
    document = json.loads(registry.read_text())
    entries = document["sources"] if isinstance(document, dict) else document
    return {e["pmid"]: f"{e['name']} ({e['status']})"
            for e in entries if e.get("pmid")}


def fetch_abstracts(ids: list[str]) -> tuple[dict[str, str], set[str]]:
    """Text per id, and the ids that had a real abstract rather than a title."""
    out: dict[str, str] = {}
    with_abstract: set[str] = set()
    for start in range(0, len(ids), 20):
        batch = ids[start:start + 20]
        if not batch:
            continue
        query = " OR ".join(
            (f"EXT_ID:{i}" if i.isdigit() else f"PMCID:{i}") for i in batch
        )
        url = f"{API}?{urllib.parse.urlencode({'query': query, 'format': 'json', 'pageSize': 25, 'resultType': 'core'})}"
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = json.loads(response.read())
        except (urllib.error.URLError, json.JSONDecodeError, TimeoutError):
            continue
        for record in payload.get("resultList", {}).get("result", []):
            abstract = record.get("abstractText")
            text = re.sub(r"<[^>]+>", " ", " ".join(filter(
                None, [record.get("title"), abstract])))
            for ident in {record.get("pmid"), record.get("pmcid"), record.get("id")}:
                if ident:
                    out[ident] = text
                    if abstract:
                        with_abstract.add(ident)
        time.sleep(PAUSE)
    return out, with_abstract


def classify(text: str, has_abstract: bool = True) -> tuple[str, str]:
    if not text:
        return "unknown", "nothing retrieved; left open"
    if REVIEW.search(text):
        return "review", "a review of other people's work; it has no rows of its own"
    if REFERENCE.search(text) and not PANEL.search(text):
        return "reference", "a database or atlas; records events per drug, not outcomes per sequence"
    if PANEL.search(text) and READOUT.search(text):
        return "dataset", "describes a screen or panel with a measured readout"
    if SINGLE.search(text):
        return "single_study", "characterises one compound or construct; cannot supply training rows"
    if METHOD.search(text):
        return "method", "a platform, chemistry or delivery paper rather than a dataset"
    if PANEL.search(text):
        return "dataset", "panel language without a clear readout; kept open to be safe"
    if not has_abstract:
        return "unknown", ("only a title was retrieved and it shows no panel or "
                           "screen language; too little to dismiss it on")
    return "single_study", "no panel or screen language; most likely one compound"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--recheck", action="store_true")
    args = parser.parse_args(argv)

    if not CANDIDATES.exists():
        raise SystemExit("no candidates; run scripts/monitor_literature.py first")
    candidates = json.loads(CANDIDATES.read_text())["candidates"]
    previous = json.loads(TRIAGE.read_text()).get("assessed", {}) if TRIAGE.exists() else {}

    todo = [c for c in candidates if args.recheck or c["id"] not in previous]
    print(f"{len(candidates)} standing, {len(todo)} to assess")
    abstracts, with_abstract = (fetch_abstracts([c["id"] for c in todo])
                                if todo else ({}, set()))

    registry = decided_pmids()
    overrides = json.loads(DECISIONS.read_text()) if DECISIONS.exists() else {}
    assessed = dict(previous)
    for candidate in todo:
        if candidate["id"] in overrides:
            call = overrides[candidate["id"]]
            verdict, reason = call["verdict"], f"read by hand: {call['reason']}"
        elif candidate["id"] in registry:
            verdict = "decided"
            reason = f"already in the source registry as {registry[candidate['id']]}"
        else:
            verdict, reason = classify(abstracts.get(candidate["id"], ""),
                                       candidate["id"] in with_abstract)
        assessed[candidate["id"]] = {
            "title": candidate["title"], "year": candidate["year"],
            "url": candidate["url"], "verdict": verdict, "reason": reason,
            "read_from": ("abstract" if candidate["id"] in with_abstract
                          else "title only" if candidate["id"] in abstracts
                          else "nothing retrieved"),
        }

    TRIAGE.write_text(json.dumps({
        "note": "`dataset` and `unknown` stay open for a human to read. The classifier errs "
                "toward dataset, because a false positive costs five minutes and a "
                "false negative loses a dataset.",
        "counts": {v: sum(1 for a in assessed.values() if a["verdict"] == v)
                   for v in sorted({a["verdict"] for a in assessed.values()})},
        "assessed": assessed,
    }, indent=2) + "\n")

    counts: dict[str, int] = {}
    for record in assessed.values():
        counts[record["verdict"]] = counts.get(record["verdict"], 0) + 1
    print()
    for verdict, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {verdict:14s} {n}")

    worth_reading = [a for a in assessed.values()
                     if a["verdict"] in {"dataset", "unknown"}]
    print(f"\n{len(worth_reading)} worth a human reading:")
    for record in worth_reading[:10]:
        print(f"  {record['year']} {record['title'][:72]}")
        print(f"      {record['url']}")
    print(f"\nwrote {TRIAGE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
