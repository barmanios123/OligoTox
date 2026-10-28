# evidence-audit

Six checks on whether a model's score on a dataset means what it appears to
mean.

A PR-AUC is a number you can always compute and rarely interpret. Each check
here exists because it caught a result that looked good and was not. They were
written while building [OligoTox](https://github.com/barmanios123/OligoTox), an
oligonucleotide toxicity model, and every one of them corresponds to an entry
in that project's `docs/NEGATIVE_RESULTS.md`. None of them is specific to
oligonucleotides, or to biology.

## What it checks

| Check | The question a PR-AUC does not answer |
|---|---|
| `leakage` | Is a feature a relabelled copy of the target? |
| `grouping` | How many independent subjects are behind the row count? |
| `source_baseline` | Does a model seeing only *which dataset a row came from* do just as well? |
| `within_source` | Does anything survive inside a single source? |
| `simple_baselines` | Does one crude number match the whole feature matrix? |
| `calibration` | Are the probabilities probabilities, or only ranks? |

The verdicts are deliberately blunt: *passes*, *insufficient data*, *source
confounded*, *no better than a simple baseline*, *ranks but uncalibrated*.

## Install

```bash
pip install -e packages/evidence-audit
```

## Use

```python
from evidence_audit import audit_dataset

report = audit_dataset(X, y, groups=subject_ids, sources=study_ids)
print(report.summary())
```

Or from the command line, where the exit status is 1 if any endpoint fails a
check, so it can gate a pipeline:

```bash
evidence-audit data.csv --label outcome --group subject_id --source study
```

Everything not named by `--label`, `--group`, `--source` or `--exclude` is
treated as a feature.

## One thing to know about `--group`

The grouping key is what keeps replicates of one subject out of both sides of
a split, and it is the check most often got wrong. Group keys are compared as
upper-cased strings and otherwise left alone.

`--nucleotide` additionally folds `U` onto `T`, so RNA and DNA spellings of one
sequence count as one group. It is off by default, because the same
substitution would turn a subject id of `USER-1` into `TSER-1` and stop it
matching itself.

## Licence

MIT.
