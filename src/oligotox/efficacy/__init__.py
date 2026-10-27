"""Efficacy prediction.

`score_efficacy` is the shipped path and depends only on the runtime pins in
requirements.txt.

The OligoFormer, ENsiRNA and siDPT classes are untrained interface stubs that
import torch, an optional extra. Importing them here eagerly made
`from oligotox.efficacy.potency import score_efficacy` pull in torch, which
broke the deployed app, so they resolve lazily: touching one still raises a
clear ImportError when torch is absent, but importing this package does not.
"""
from .potency import load_efficacy_config, score_efficacy

__all__ = ["score_efficacy", "load_efficacy_config", "OligoFormer", "ENsiRNA", "siDPT"]

_LAZY = {
    "OligoFormer": ("oligoformer", "OligoFormer"),
    "ENsiRNA": ("ensirna", "ENsiRNA"),
    "siDPT": ("sidpt", "siDPT"),
}


def __getattr__(name: str):
    if name not in _LAZY:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr = _LAZY[name]
    try:
        module = __import__(f"{__name__}.{module_name}", fromlist=[attr])
    except ImportError as exc:  # pragma: no cover - depends on the optional extra
        raise ImportError(
            f"{name} is an untrained interface stub needing the optional deep learning "
            'extra. Install it with: pip install -e ".[efficacy]"'
        ) from exc
    return getattr(module, attr)


def __dir__():
    return sorted(__all__)
