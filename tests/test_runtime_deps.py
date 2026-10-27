"""The deployed app may only import what requirements.txt pins.

A deploy failed because oligotox.efficacy.__init__ eagerly imported an
untrained stub that imports torch, an optional extra. Importing
`oligotox.efficacy.potency` was enough to pull it in.
"""
import importlib
import sys

import pytest

OPTIONAL = ("torch", "RNA", "transformers", "plotly", "fastapi")

# Everything oligotox_app.py imports at module level.
APP_IMPORTS = (
    "oligotox.features.candidate",
    "oligotox.pipelines.score_candidate",
    "oligotox.pkpd.bridge",
    "oligotox.reporting.oswg_report",
    "oligotox.efficacy.potency",
    "oligotox.efficacy",
)


class _Blocker:
    """Make the optional extras look uninstalled."""

    def find_module(self, name, path=None):
        return self.find_spec(name, path)

    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in OPTIONAL:
            raise ImportError(f"{name} is blocked: not pinned in requirements.txt")
        return None


@pytest.fixture
def without_optional_extras():
    blocker = _Blocker()
    saved = {k: v for k, v in sys.modules.items()
             if k.split(".")[0] in OPTIONAL or k.startswith("oligotox")}
    for name in saved:
        del sys.modules[name]
    sys.meta_path.insert(0, blocker)
    try:
        yield
    finally:
        sys.meta_path.remove(blocker)
        for name in [k for k in sys.modules if k.startswith("oligotox")]:
            del sys.modules[name]
        sys.modules.update(saved)


@pytest.mark.parametrize("module", APP_IMPORTS)
def test_app_imports_without_optional_extras(without_optional_extras, module):
    importlib.import_module(module)


def test_the_untrained_stubs_say_what_they_need(without_optional_extras):
    efficacy = importlib.import_module("oligotox.efficacy")
    # Reachable by name, but only with the extra installed.
    with pytest.raises(ImportError, match="efficacy"):
        efficacy.OligoFormer
    # The shipped path is unaffected.
    assert callable(efficacy.score_efficacy)
