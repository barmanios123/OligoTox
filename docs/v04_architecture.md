# OligoTox v0.4 Architecture

OligoTox v0.4 adds a decision-support layer on top of the v0.3 mechanistic and ML stack.

## Modules

- `qAOP/`: Quantitative adverse outcome pathways.
- `redesign/`: Counterfactual chemistry and delivery redesign suggestions.
- `impurity_tox/`: Sequence variant and impurity liability assessment.
- `n_of_1/`: Human off-target screening and minimum safety package generation.
- `translational_failure/`: Preclinical-to-clinical translational gap scoring.

## Purpose

The v0.4 layer is designed to answer not only "what is the risk?" but also:

- Why is the risk elevated?
- What change may reduce risk?
- What impurity variants could be more dangerous?
- What minimum human-focused safety package is needed?
- Where is translation from preclinical species to humans most fragile?

## Status

This layer is currently heuristic and architecture-first.
It is intended to be strengthened with curated datasets in the next phase.
