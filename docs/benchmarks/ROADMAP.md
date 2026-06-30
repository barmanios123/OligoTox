# Benchmarking and Model Upgrade Roadmap

Phase 1
- Normalize OligoGym-style datasets into a common schema.
- Implement leakage-aware split generation.
- Benchmark current OligoTox baseline.

Phase 2
- Add RNA-FM and RiNALMo embedding feature extraction.
- Compare embedding-augmented baselines against current models.

Phase 3
- Integrate external model adapters for OligoFormer and ENsiRNA.
- Publish standardized benchmark report and comparison tables.

Phase 1.5
- Generate leakage-aware train/val/test splits grouped by sequence family or target gene.
- Verify zero group overlap across splits before any benchmark metrics are reported.

Phase 1.75
- Run a simple reproducible baseline on normalized, leakage-aware benchmark splits.
- Save metrics, predictions, and error analysis artifacts before adding advanced models.

Phase 1.9
- Add post-run error analysis for false positives, false negatives, and hardest benchmark cases.
- Inspect sequence-length, GC, and motif-linked failure patterns before upgrading models.
