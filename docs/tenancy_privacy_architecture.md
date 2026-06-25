# OligoTox Tenancy & Data Privacy Architecture

## Current Deployment: Single-Tenant (Local)
All data, models, and artifacts are stored in `~/oligotox/` and **never transmitted externally**.
Sequence data never appears in logs, URLs, or network calls.

## Multi-Tenant SaaS Roadmap
| Layer | Design |
|---|---|
| Data | Per-tenant S3/Azure Blob bucket with tenant-specific KMS key |
| Compute | Per-tenant isolated container (ECS/ACI) for training jobs |
| Models | Per-tenant registry; public baselines are read-only shared layer |
| API | JWT-authenticated endpoints; sequences never in logs |
| Audit | Immutable run log with dataset fingerprints only — no raw sequences |

## Data Privacy Guarantees
- Sequences are **never logged, transmitted, or included in URLs**
- Partner data trains **partner-specific models only** — never merged into public baseline without consent
- Audit trail stores **SHA-256 fingerprints** of datasets, not datasets themselves
- Exported `.joblib` model files contain **learned parameters only** — sequence reconstruction is infeasible

## Regulatory Alignment
| Framework | Status |
|---|---|
| 21 CFR Part 11 | Audit trail satisfies electronic records requirements with user authentication |
| GDPR / CCPA | No personal data processed; oligonucleotide sequences are not personal data |
| ICH M10 | Not applicable (computational predictions, not bioanalytical validation) |
| FDA AI/ML SaMD guidance | Conforms via model versioning, change control registry, CV AUROC monitoring |

## Partner Privacy Workflow (NDA/CDA-Safe)
1. Partner downloads `src/oligotox/partner_import.py` (runs entirely inside their firewall)
2. Partner runs import + validation + training locally — **no data leaves their environment**
3. Partner shares only: the audit report JSON + trained model binary
4. OligoTox integrates model artifact into benchmark comparison
5. Raw sequences stay inside the partner's firewall at all times
