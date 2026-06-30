# OligoTox

OligoTox is a mechanism-aware safety, translational interpretation, and candidate triage platform for oligonucleotide therapeutics.

## What OligoTox does

The platform helps evaluate oligonucleotide candidates across key toxicity and translation risk domains relevant to ASOs, gapmers, siRNAs, lncRNA-targeting oligos, and related modalities. Its scientific framing is organized around three mechanistic classes:

- **D1:** exaggerated on-target pharmacology
- **D2:** hybridization-dependent off-target effects
- **D3:** chemistry, formulation, conjugate, and delivery-class liabilities

## What is working now

The current Streamlit app supports:

- Single-sequence toxicity analysis
- Modality, chemistry, and conjugate-aware input
- Multi-endpoint risk prediction through the OligoTox inference backend
- Confidence and uncertainty proxy display
- Intermediate safety calculations such as GC content, hybridization risk, innate immune proxy, platelet/PS proxy, hepatic gapmer proxy, and renal accumulation proxy
- Optional PK/PD context using dose and schedule inputs
- Assay recommendation generation
- Candidate comparison scaffolding
- Regulatory-style safety summary output
- JSON result export

## App entrypoint

Main app:

- `oligotox_app.py`

Core package code:

- `src/oligotox/`

Supporting areas include:

- `docs/` — architecture, status, and strategy documents
- `data/` — processed and source data
- `tests/` — test suite
- `scripts/` — utilities

## Run locally

### Launch the app

```bash
streamlit run oligotox_app.py --server.port 8502
```

Then open:

```text
http://localhost:8502
```

### Run training

```bash
# training entrypoint moved or removed from minimal runtime repo
```

If you prefer the package entrypoint and it is wired in your environment:

```bash
python src/oligotox/training/train.py
```

### Run tests

```bash
pytest
```

## How to use the app

1. Open the **Sequence Analysis** tab.
2. Paste an oligonucleotide sequence.
3. Choose modality, chemistry, and conjugate.
4. Optionally enter dose and schedule for PK/PD context.
5. Run prediction to view endpoint risks, mechanistic support metrics, confidence display, assay suggestions, and the regulatory-style summary.

## What was updated today

Today’s work focused on stabilizing and extending the app:

- Restored and verified a working `oligotox_app.py`
- Reconnected the Streamlit app to the `run_oligotox` inference backend
- Added modality-aware result context
- Added intermediate safety interpretation outputs
- Added assay recommendation logic
- Added comparator/report scaffolding
- Verified the app compiles and launches successfully
- Pushed the updated app state to GitHub

## Current status

The app is running and the latest Streamlit enhancement commit has been pushed to GitHub. The platform is functional for demo and iterative development, but some components are still evolving and should be described according to what is currently verified in code and local execution.

## Development note

Keep root-level clutter minimal. Longer status notes, architecture writeups, and provenance documents should live under `docs/` where possible.

---

## 📄 License

**Proprietary and Confidential**

This software is proprietary and confidential. It is provided for **evaluation purposes only** to authorized clients and collaborators.

- ❌ No commercial use without license
- ❌ No modification or derivative works
- ❌ No distribution or sublicensing
- ✅ Evaluation and collaboration discussions only

For commercial licensing inquiries, please contact the repository owner.

See [LICENSE](LICENSE) file for complete terms.

