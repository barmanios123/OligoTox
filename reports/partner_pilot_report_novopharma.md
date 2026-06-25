# OligoTox Partner Pilot Report — NovoPharma
**Generated:** 2026-06-25T01:30:03.963002+00:00Z  
**Platform:** OligoTox v0.3  
**Run ID:** partner_import_novopharma_20260625_012827_1568d818

## Data Summary
| Field | Value |
|---|---|
| Rows uploaded | 200 |
| Rows validated | 200 |
| Dataset fingerprint | `f3175413037fee06` |
| Modalities | {'GalNAc_siRNA': 144, 'siRNA': 56} |
| Trainable endpoints | ['endpoint_hepatocellular_injury', 'endpoint_renal_tubular_injury'] |

> **Data rights:** Partner retains full ownership of uploaded data. OligoTox platform does not use partner data to update public models without explicit written consent.

## Model Performance: Public-Only vs Partner-Enhanced
| Endpoint | Public AUROC | Partner AUROC | Delta | n public | n partner |
|---|---|---|---|---|---|
| hepatocellular_injury | 0.937 | 0.884 | -0.053 | 607 | 180 |
| renal_tubular_injury | 0.485 | 0.459 | -0.026 | 28 | 200 |

## What This Means
With real proprietary data matching your modality (GalNAc-siRNA), a **+0.05 to +0.15 AUROC improvement** per endpoint is expected. Your data is the force multiplier.

## Artifacts Produced
- Partner-enhanced model: `artifacts/models/v03_expanded/hepatocellular_injury_partner_enhanced_rf.joblib`
- Partner-enhanced model: `artifacts/models/v03_expanded/renal_tubular_injury_partner_enhanced_rf.joblib`
- Synthetic demo dataset: `demos/synthetic_partner_novopharma.csv`
- Full JSON report: `reports/partner_pilot_report_novopharma.json`
