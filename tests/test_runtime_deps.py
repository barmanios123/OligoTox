"""The deployed app may only import what requirements.txt pins.

A deploy failed because oligotox.efficacy.__init__ eagerly imported an
untrained stub that imports torch, an optional extra. Importing
`oligotox.efficacy.potency` was enough to pull it in.
"""
import importlib
import sys

import pytest

# ViennaRNA (`RNA`) used to be here. It is required now: structure_features
# raises without it rather than falling back to a GC placeholder, so blocking
# it no longer tests anything real, and leaving it listed as optional is how
# its absence from pyproject went unnoticed.
OPTIONAL = ("torch", "transformers", "plotly", "fastapi")

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


def test_pyproject_declares_everything_requirements_pins():
    """`pip install .` must give a working package, not one missing a hard dep.

    requirements.txt pins what the deployed app and CI install; pyproject's
    core dependencies are what anyone installing the package gets. When a
    requirement is added to only the first, both environments keep working and
    nothing fails until someone installs the package on its own. That is how
    ViennaRNA came to be pinned in requirements.txt and absent from pyproject
    while the feature path hard-required it.
    """
    import re
    import tomllib
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent
    pinned = {
        re.split(r"[=<>!~]", line, maxsplit=1)[0].strip().lower()
        for line in (root / "requirements.txt").read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]
    declared = {
        re.split(r"[=<>!~\[;]", dep, maxsplit=1)[0].strip().lower()
        for group in [project.get("dependencies", [])]
        + list(project.get("optional-dependencies", {}).values())
        for dep in group
    }
    assert not (pinned - declared), (
        f"pinned in requirements.txt but undeclared in pyproject: "
        f"{sorted(pinned - declared)}"
    )
