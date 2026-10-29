"""A failed request must never be recorded as a gene that does not exist.

This is the bug these pin down. `_get` caught HTTPError, URLError,
JSONDecodeError, TimeoutError and OSError in one tuple and the caller turned
any of them into "Ensembl cannot place this gene", then wrote that to a
blacklist later runs skip. One rate-limited window removed 151 of 350 genes
from the corpus permanently. APOC3 and BRCA1 were both on that list and both
resolve on the first attempt when asked again.

So the distinction under test is not cosmetic. A 404 is a fact about the
symbol. A 429 or a timeout is a fact about the network, and acting on it as
though it were the first silently deletes data.
"""
from __future__ import annotations

import json
import sys
import urllib.error
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import fetch_transcripts as ft  # noqa: E402


def _http_error(code: int, headers: dict | None = None) -> urllib.error.HTTPError:
    return urllib.error.HTTPError(
        url="https://rest.ensembl.org/x", code=code, msg="test",
        hdrs=headers or {}, fp=None,
    )


@pytest.fixture(autouse=True)
def no_sleeping(monkeypatch):
    """Retries back off for real seconds; the tests must not."""
    monkeypatch.setattr(ft.time, "sleep", lambda *_: None)


def test_a_404_is_permanent(monkeypatch):
    """Ensembl has no such symbol, so asking again is pointless."""
    def boom(*_a, **_k):
        raise _http_error(404)
    monkeypatch.setattr(ft.urllib.request, "urlopen", boom)
    with pytest.raises(ft.Unknown):
        ft._get("/lookup/symbol/homo_sapiens/NOTAGENE")


def test_rate_limiting_is_transient_not_unknown(monkeypatch):
    """The case that deleted 151 genes."""
    def boom(*_a, **_k):
        raise _http_error(429, {"Retry-After": "1"})
    monkeypatch.setattr(ft.urllib.request, "urlopen", boom)
    with pytest.raises(ft.Transient):
        ft._get("/lookup/symbol/homo_sapiens/APOC3")


def test_a_server_error_is_transient(monkeypatch):
    def boom(*_a, **_k):
        raise _http_error(503)
    monkeypatch.setattr(ft.urllib.request, "urlopen", boom)
    with pytest.raises(ft.Transient):
        ft._get("/x")


def test_a_timeout_is_transient(monkeypatch):
    def boom(*_a, **_k):
        raise TimeoutError("timed out")
    monkeypatch.setattr(ft.urllib.request, "urlopen", boom)
    with pytest.raises(ft.Transient):
        ft._get("/x")


def test_a_retry_succeeds_after_a_transient_failure(monkeypatch):
    """Backoff has to actually recover, or it is just a slower failure."""
    calls = {"n": 0}

    class _Response:
        def __enter__(self):
            return self

        def __exit__(self, *_):
            return False

        def read(self):
            return json.dumps({"id": "ENSG1", "Transcript": []}).encode()

    def flaky(*_a, **_k):
        calls["n"] += 1
        if calls["n"] < 3:
            raise _http_error(429, {"Retry-After": "0"})
        return _Response()

    monkeypatch.setattr(ft.urllib.request, "urlopen", flaky)
    assert ft._get("/x")["id"] == "ENSG1"
    assert calls["n"] == 3, "should have retried rather than given up"


def test_resolve_reports_why_it_failed(monkeypatch):
    """The caller needs the reason, since only one kind may be remembered."""
    monkeypatch.setattr(ft, "_get", lambda *a, **k: (_ for _ in ()).throw(ft.Transient("429")))
    record, reason = ft.resolve("APOC3")
    assert record is None and reason == "transient"

    monkeypatch.setattr(ft, "_get", lambda *a, **k: (_ for _ in ()).throw(ft.Unknown("404")))
    record, reason = ft.resolve("NOTAGENE")
    assert record is None and reason == "unknown"


def _run_main(monkeypatch, tmp_path, outcomes: dict[str, tuple]):
    """Drive the real main() over a fake corpus, returning the manifest."""
    import pandas as pd

    cache = tmp_path / "genomic"
    monkeypatch.setattr(ft, "GENOMIC_CACHE", cache)
    monkeypatch.setattr(ft, "CACHE", tmp_path / "cdna")
    monkeypatch.setattr(ft, "ATLAS", tmp_path / "atlas.parquet")
    (tmp_path / "atlas.parquet").write_text("")
    monkeypatch.setattr(ft, "ROOT", tmp_path)
    monkeypatch.setattr(
        pd, "read_parquet",
        lambda *a, **k: pd.DataFrame({"target_gene": list(outcomes)}),
    )
    monkeypatch.setattr(ft, "resolve", lambda sym: outcomes[sym])

    ft.main(["--genomic"])
    return json.loads((cache / "_manifest.json").read_text())


def test_a_transient_failure_is_never_blacklisted(monkeypatch, tmp_path):
    """The end-to-end guarantee, driven through the real main().

    Only a confirmed unknown goes into `unknown`, which is the set later runs
    skip. If this regresses, one bad network minute silently shrinks the
    corpus and every downstream result is computed on whatever survived.
    """
    manifest = _run_main(monkeypatch, tmp_path, {
        "APOC3": (None, "transient"),
        "NOTAGENE": (None, "unknown"),
    })

    assert "NOTAGENE" in manifest["unknown"], "a 404 should be remembered"
    assert "APOC3" not in manifest["unknown"], (
        "a transiently failed gene must stay retryable, not be blacklisted"
    )


def test_a_recovered_gene_is_cleared_from_both_lists(monkeypatch, tmp_path):
    """A gene that resolves must stop being listed as a failure."""
    record = {"symbol": "APOC3", "gene_id": "ENSG00000110245", "length": 3000,
              "sequence": "ACGT", "resolved_as": "APOC3", "transcript_id": "T1",
              "canonical": True, "kind": "genomic", "truncated": False}
    manifest = _run_main(monkeypatch, tmp_path, {"APOC3": (record, None)})

    assert "APOC3" in manifest["hits"]
    assert "APOC3" not in manifest.get("misses", [])
    assert "APOC3" not in manifest.get("unknown", [])
