# OligoTox Architecture

## Source of truth

All importable Python application code lives in:

- `src/oligotox/`

Use package imports like:

```python
from oligotox.pipelines import build_training_dataset
from oligotox.features.aso_optimizer_v03 import batch_aso_features
from oligotox.omics import OMICS_REGISTRY
```

Do not import from old root-level package names such as:

```python
from pipelines import ...
from features import ...
from omics import ...
from config import ...
```

## Repository layout

- `src/oligotox/` — reusable application and library code
- `scripts/` — thin command-line wrappers
- `tests/` — automated tests
- `app/`, `app_mvp/` — UI/application surfaces
- `data/`, `docs/`, `notebooks/`, `reports/`, `artifacts/` — project assets and outputs
- `archive/` — archived legacy code and historical material

## Rule of thumb

If code is meant to be imported by other code, it belongs in `src/oligotox/`.

If code is just a runnable helper, it belongs in `scripts/` and should import from `oligotox.*`.
