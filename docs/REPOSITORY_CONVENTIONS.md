# Repository Conventions

This document clarifies the intended meaning of the main repository areas and reduces ambiguity around source code, generated outputs, and local-only material.

## Intended public surface

The most important public entrypoints are:

- `oligotox_app.py` — app entrypoint
- `tests/` — pytest suite (run with `pytest`)
- `src/oligotox/` — package implementation
- `scripts/` — operational helpers
- `docs/` — explanatory documentation

## Source code vs local/generated material

### Source code

These areas should be treated as the main code surface:

- `src/`
- `scripts/`
- `oligotox_app.py`
- `tests/`
- `configs/`
- selected docs in `docs/`

### Local or generated material

These areas may contain local runtime outputs, generated files, or workflow-specific artifacts:

- `artifacts/`
- `data/`
- `reports/`
- `output/`
- `external/`
- `_local_only/`
- `app_backups/`

Not every file in those areas should be treated as curated source-controlled product code.

## App convention

The canonical app entrypoint is:

```bash
bash scripts/run_app.sh
```

This launches:

```bash
streamlit run oligotox_app.py --server.port 8502
```

## Validation convention

The canonical validation command is:

```bash
pytest -q
```

The broader operational health command is:

```bash
bash scripts/project_health.sh
```
