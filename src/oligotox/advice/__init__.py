"""Design advice: what to change about a sequence, and what it costs.

Two things this deliberately does not do.

It does not advise from unvalidated endpoints. Ten of the thirteen toxicity
models are not presented as predictions, and an edit recommended to lower a
number that means nothing is worse than no advice at all. Only endpoints
tiered `validated` contribute, and the advice says which.

It does not read coefficients. The models are linear, so coefficients are
available, but with correlated features they describe the fit rather than
the effect of an edit. Every suggestion here is measured: the sequence is
actually mutated, actually re-scored through the same pipeline the app uses,
and the reported change is the difference between two real predictions.

    from oligotox.advice import explain, suggest
"""
from .advisor import Edit, SequenceAdvice, explain, suggest

__all__ = ["explain", "suggest", "Edit", "SequenceAdvice"]
