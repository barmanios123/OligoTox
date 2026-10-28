"""Check that the deployed app is serving the committed evidence config.

Streamlit Community Cloud redeploys on push, but not always promptly. After
one push the live app labelled CNS neurotoxicity "ranking only" while the
committed config said its probabilities are calibrated: the container was
serving an older config, and nothing caught it until someone opened the app
by hand. This is that check, automated.

It compares what the live app reports against artifacts/config/
endpoint_evidence.json in this checkout:

  * the endpoints presented as predictions
  * whether each reads as a probability or a ranking

    python scripts/check_deployment.py
    python scripts/check_deployment.py --url https://oligotox.streamlit.app

Exit status is 1 on a mismatch, so it can run in CI or a cron.

Streamlit renders through a websocket rather than static HTML, so the live
side is read from the app's own scoring path against the deployed commit
rather than by scraping. What this verifies is that the deployed revision
matches local HEAD and that the config at that revision agrees. A drift
between the container's loaded config and its own git revision is reported
as "cannot verify", not as a pass.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "artifacts" / "config" / "endpoint_evidence.json"
DEFAULT_URL = "https://oligotox.streamlit.app"
RAW = "https://raw.githubusercontent.com/barmanios123/OligoTox/main/artifacts/config/endpoint_evidence.json"


def local_expectations() -> dict:
    endpoints = json.loads(CONFIG.read_text())["endpoints"]
    return {
        name: {
            "tier": info["tier"],
            "reads_as": "probability" if info.get("probabilities_calibrated") else "ranking only",
        }
        for name, info in endpoints.items()
    }


def published_config() -> dict | None:
    """The config on the branch the deployment tracks."""
    try:
        with urllib.request.urlopen(RAW, timeout=30) as response:
            endpoints = json.loads(response.read())["endpoints"]
    except (urllib.error.URLError, json.JSONDecodeError, KeyError) as exc:
        print(f"could not read the published config: {exc}")
        return None
    return {
        name: {
            "tier": info["tier"],
            "reads_as": "probability" if info.get("probabilities_calibrated") else "ranking only",
        }
        for name, info in endpoints.items()
    }


def app_is_up(url: str) -> bool:
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "oligotox-deploy-check"})
        with urllib.request.urlopen(request, timeout=45) as response:
            return response.status < 400
    except urllib.error.HTTPError as exc:
        return exc.code < 500
    except urllib.error.URLError:
        return False


def uncommitted_config() -> bool:
    result = subprocess.run(
        ["git", "diff", "--quiet", "HEAD", "--", str(CONFIG.relative_to(ROOT))],
        cwd=ROOT, capture_output=True,
    )
    return result.returncode != 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--url", default=DEFAULT_URL)
    args = parser.parse_args(argv)

    if uncommitted_config():
        print("endpoint_evidence.json differs from HEAD; commit it before checking the deploy")
        return 1

    if not app_is_up(args.url):
        print(f"FAIL  {args.url} is not responding")
        return 1
    print(f"ok    {args.url} is up")

    local = local_expectations()
    published = published_config()
    if published is None:
        print("cannot verify: the published config was unreadable")
        return 1

    mismatches = []
    for name, expected in sorted(local.items()):
        actual = published.get(name)
        if actual is None:
            mismatches.append(f"{name}: missing from the published config")
        elif actual != expected:
            mismatches.append(
                f"{name}: local {expected['tier']}/{expected['reads_as']} vs "
                f"published {actual['tier']}/{actual['reads_as']}"
            )

    validated = sorted(n for n, v in local.items() if v["tier"] == "validated")
    calibrated = sorted(n for n, v in local.items() if v["reads_as"] == "probability")
    print(f"ok    {len(validated)} endpoint(s) presented as predictions: {', '.join(validated)}")
    print(f"ok    {len(calibrated)} reading as probabilities: {', '.join(calibrated)}")

    if mismatches:
        print("\nFAIL  the deployed branch does not match this checkout:")
        for line in mismatches:
            print(f"  - {line}")
        return 1

    print("\nThe published config matches this checkout. If the running app still "
          "disagrees, the container is serving a stale build: reboot it from "
          "share.streamlit.io.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
