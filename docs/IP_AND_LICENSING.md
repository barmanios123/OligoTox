# IP & Licensing Strategy — OligoTox

> **Status:** Proprietary platform · All rights reserved · Version 0.1.0

---

## 1. What OligoTox Is (IP Perspective)

OligoTox is a **software + curated-data + trained-model platform** for
sequence-aware oligonucleotide toxicity prediction. Unlike a text report or
generic ML notebook, it combines:

| Layer | IP Form |
|---|---|
| Source code (pipeline, API, app) | Copyright + Trade secret |
| Curated sequence->candidate mappings | Trade secret database |
| Trained model weights & configs | Trade secret |
| Feature-engineering recipes | Trade secret / potential patent |
| Labeled endpoint dataset | Proprietary compiled dataset |

Each layer can be protected independently; no single patent is required to
derive commercial value.

---

## 2. Licensing — What Is Possible Today

You can **license OligoTox for money or in-kind value right now**, without
any patents filed.

### Why
- You own the code, data curation, and model configuration (absent conflicting
  employer/university IP agreements).
- A granted patent is **not** a prerequisite for charging for software or
  services — it only adds enforcement leverage later.
- For early-stage AI drug-discovery tools the dominant deal structure is
  **access to a hosted tool or analytical services**, with models and training
  data kept as trade secrets.

### Realistic Near-Term License Structures

#### A. Services / Pilot Access (lowest friction)
> "I run OligoTox analyses on your oligo series and deliver structured reports;
> you do not receive raw code or model weights."

- Wrap in a simple **services agreement** (NDA + scope of work + IP ownership clause).
- You retain all platform rights; client owns the reports/results for internal use only.
- Pricing anchor: **$5k–$20k** for a 3–6 month pilot (small biotech or academic lab scale).

#### B. Collaboration / Sponsored Research
> Company or PI provides proprietary oligo data -> you extend OligoTox for
> their endpoints -> they compensate you (cash, equity, co-authorship, or co-IP).

- Agreement specifies who owns jointly developed improvements.
- Useful when a partner has scarce labeled data you need to grow the model.

#### C. On-Premise / Source-Code License (later stage)
- Higher ticket ($50k–$500k range for small biotech), but requires:
  - Reproducible deployment (Docker/CLI — already partially satisfied).
  - Clear documentation of inputs, outputs, and known limitations.
  - Validation story (benchmarks, case studies).

---

## 3. What Makes OligoTox Potentially Patentable

Software/AI patents require a **specific, concrete technical method** that
improves computer technology or another technical field — not just
"using ML to predict toxicity."

Candidate technical elements in OligoTox that may support a patent claim
(pending prior-art search):

1. **Sequence-recovery & mapping workflow** — the multi-source join logic,
   validation pipeline, and audit steps for heterogeneous oligo toxicity
   datasets with partial or missing sequence information.

2. **Chemistry-aware feature augmentation** — combining sequence-derived
   features with modality-specific chemical descriptors, then applying
   conservative feature/label dropping to preserve signal in extremely sparse
   labeled settings.

3. **Oligo-specific model architecture / training regime** — if the
   multi-task endpoint structure, calibration behavior, or uncertainty
   handling is sufficiently distinct from generic tox-ML literature
   (Tox21 models, MolToxPred, commercial AI-tox services).

A patent attorney would frame this as:
> *"A method and system for constructing and training sequence-aware
> oligonucleotide toxicity models from heterogeneous labeled sources with
> partial or missing sequence information."*

**Action trigger:** If a company wants to adopt or co-develop OligoTox,
engage an IP attorney for a prior-art search before disclosure.

---

## 4. Trade-Secret Strategy (Default Until Patents Filed)

- Keep this repository **private** on GitHub.
- Share only high-level descriptions and non-sensitive demos in interviews
  and public materials.
- Any demo should use synthetic/anonymized oligo sequences, not your curated
  mapping files or training data.
- Watermark or timestamp critical methodology documents (commit history
  serves as timestamped evidence of inventorship).

---

## 5. Gaps Before a Commercial Product License

To command premium pricing as a stand-alone product, partners will expect:

- [ ] **Input/output spec** — documented API schema (what sequences/formats
  go in, what endpoint scores come out). *(In progress: `api.py`)*
- [ ] **Validation story** — internal benchmarks reproducing known tox
  patterns across labeled endpoints.
- [ ] **Limitation disclosure** — current dataset size, sequence coverage
  gaps, modality restrictions.
- [ ] **Reproducible deployment** — Docker image or conda env + CLI entry
  point. *(Partially done via `pyproject.toml`)*

---

## 6. Practical Next Steps

| Timeline | Action |
|---|---|
| **This month** | Keep repo private; prepare a 2-3 page OligoTox one-pager (problem -> approach -> what works -> limitations -> pilot pricing) |
| **First serious inquiry** | Engage IP attorney for prior-art search; prepare NDA template |
| **First paid pilot** | Services agreement -> NDA + scope of work + IP ownership clause |
| **Co-development deal** | Sponsored-research-style agreement; negotiate co-IP terms upfront |
| **If novel enough** | Provisional patent filing (~$1,500-$3,000 attorney-assisted) to establish priority date |

---

## 7. Current License

This repository is released under a **proprietary license**.
See `LICENSE` in the repo root.

No part of this codebase, trained models, curated datasets, or derived works
may be reproduced, distributed, or used commercially without explicit written
permission from the author.

---

*Document authored by: Beshoy Armanios*  
*Last updated: June 2026*