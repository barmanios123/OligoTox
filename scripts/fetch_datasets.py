"""Check the declared sources and report what the pipeline can ingest.

Run unattended, before the pipeline, so a scheduled job knows what it has.

A source is usable only if `data/registry/sources.json` declares its readout
direction and label threshold. That is deliberate: a downloader must not
decide what counts as toxic. A source with data available but no declared
direction is reported as blocked, not guessed at.

    python scripts/fetch_datasets.py            # report
    python scripts/fetch_datasets.py --fetch    # also download what it can

Statuses in the registry:

    present     already on disk; verified here
    fetchable   obtainable, but needs a step a script cannot do unattended,
                such as a signed URL or a manual download
    requested   asked for from the authors, not released
    unverified  availability not confirmed
    deferred    available but not yet useful

Exit status is 0 whenever the report completes. Missing optional sources are
not a failure; this is an inventory, not a gate.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "data" / "registry" / "sources.json"
RESULTS = ROOT / "reports" / "results"

USABLE_FOR_TRAINING = {"present"}


def _present(source: dict) -> bool:
    path = source.get("local_path")
    return bool(path) and (ROOT / path).exists()


def _reachable(url: str) -> str:
    try:
        request = urllib.request.Request(url, method="HEAD",
                                         headers={"User-Agent": "oligotox-fetch"})
        with urllib.request.urlopen(request, timeout=20) as response:
            return f"reachable ({response.status})"
    except urllib.error.HTTPError as exc:
        return f"HTTP {exc.code}"
    except Exception as exc:  # noqa: BLE001 - an inventory should never crash
        return f"unreachable ({type(exc).__name__})"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--fetch", action="store_true",
                        help="probe remote sources as well as reporting local ones")
    args = parser.parse_args(argv)

    registry = json.loads(REGISTRY.read_text())
    report = []

    for source in registry["sources"]:
        name, status, kind = source["name"], source["status"], source["kind"]
        on_disk = _present(source)
        usable = status in USABLE_FOR_TRAINING and on_disk

        blocked = None
        if source.get("role") == "external_check":
            # A held-out screen is what the endpoint is tested against.
            # Counting it as training data would both overstate the corpus
            # and describe the one honest check as if it were fitted on.
            blocked = "held out as an external check"
        elif kind == "toxicity" and source.get("direction") is None:
            blocked = "no readout direction declared"
        elif kind == "toxicity" and source.get("threshold") is None:
            blocked = "no label threshold declared"
        elif status in USABLE_FOR_TRAINING and not on_disk:
            blocked = f"declared present but missing at {source.get('local_path')}"
        elif status != "present":
            blocked = f"status {status}"

        entry = {
            "name": name, "status": status, "kind": kind,
            "endpoint": source.get("endpoint"), "on_disk": on_disk,
            "usable_for_training": usable and blocked is None,
            "blocked_by": blocked, "sequences": source.get("sequences"),
        }
        if args.fetch and source.get("url") and not on_disk:
            entry["remote"] = _reachable(source["url"])
        report.append(entry)

    usable = [r for r in report if r["usable_for_training"]]
    blocked = [r for r in report if not r["usable_for_training"]]
    total_sequences = sum(r["sequences"] or 0 for r in usable)

    print(f"{len(usable)} source(s) usable now, {total_sequences} sequences")
    for entry in usable:
        print(f"  ok       {entry['name']:34s} {entry['kind']:9s} "
              f"{entry['sequences'] or '?'} seq")
    print(f"\n{len(blocked)} source(s) not usable")
    for entry in blocked:
        extra = f"  [{entry['remote']}]" if entry.get("remote") else ""
        print(f"  blocked  {entry['name']:34s} {entry['blocked_by']}{extra}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "dataset_inventory.json").write_text(json.dumps({
        "usable_sources": len(usable),
        "usable_sequences": total_sequences,
        "sources": report,
    }, indent=2) + "\n")
    print(f"\nwrote {(RESULTS / 'dataset_inventory.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
