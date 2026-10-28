# Negative results in oligonucleotide toxicity modelling

Failure modes found in this project, with the numbers. They are recorded
because none of them are published elsewhere: papers and platforms report the
endpoints that worked, so the same mistakes get made repeatedly on the same
kind of data.

Every result below is reproducible from a clone. The script that produced it
is named, and `README.md` has the commands.

---

## 1. Curated annotation columns leak the outcome

**What happened.** The first models scored PR-AUC 0.93 on hepatocellular
injury. 23 `label_*` / `endpoint_*` columns in the feature matrix were
curated annotations of the outcome; `label_kupffer_activation` correlated
0.98 with its own target. Two assay readouts, `cas_avg` and `kd_avg`, were
also features despite being measured after dosing, so they do not exist for a
design being screened.

**After removing them:** hepatocellular injury 0.93 → 0.29.

**Why it generalises.** Curated toxicity datasets are assembled by reading
papers and recording what they report. The annotation and the label come from
the same sentence. Any column whose value was decided by a human who already
knew the outcome is a label in disguise.

**Check for it:** correlate every feature with every target before training.
Anything above ~0.9 is an annotation, not a feature.
Script: `src/oligotox/training/train.py` (`_drop_leaky_features`).

---

## 2. Replicate rows make random splits meaningless

**What happened.** 1,614 rows covered only 74 unique sequences. A random
row-level train/test split put replicate measurements of the same oligo on
both sides.

**Row-level vs sequence-grouped CV**, same data, same model:

| Endpoint | Row-level | Sequence-grouped |
|---|---|---|
| kupffer cell activation | 0.73 | 0.48 |
| hepatocellular injury | 0.55 | 0.29 |
| hepatic vacuolation | 0.44 | 0.16 |

**Why it generalises.** Oligo datasets are usually dose series, replicates or
timepoints. Row count flatters the sample size: 1,614 rows is 74 independent
observations.

**Check for it:** count unique sequences, not rows, and group folds on the
normalised sequence. Treat RNA and DNA spellings as the same sequence.

---

## 3. The data source is almost the label in curated sets

**What happened.** Kupffer cell activation looked like the one strong
endpoint at 0.48 PR-AUC. A baseline seeing *only which dataset a row came
from* scored **0.69**, beating the real model.

The mechanism, measured: 1,572 of 1,614 rows come from two systematic
sources. The other 42 were curated one at a time from papers reporting a
toxicity finding, so **36 of those 42 carry a positive label, against 214 of
1,572** in the systematic sources. All 34 Kupffer positives sit in the
42-row bucket. "Which source is this row from" is close to the label by
construction.

**Why it generalises.** Any dataset mixing systematic screens with
hand-curated literature examples has this. The curated rows exist *because*
something was reported about them.

**Check for it:** add a baseline whose only feature is a one-hot of the
source. If it matches your model, your model is reading the source. Then
score within a single source and see if anything survives.
Script: `scripts/evaluate_baselines.py`.

---

## 4. A hash "embedding" is not a representation, and can hurt

**What happened.** 16 of 55 features were `embed_stub_hash_rna_00..15`, a
sha256 of the sequence, described in the project roadmap as integrated RNA
foundation model embeddings. Ablating the block **improved 9 of 12
endpoints**, mean +0.009 PR-AUC, off-target gene silencing by +0.073 — enough
to move that endpoint from `weak` to `validated`.

As a control on an independent split, the hash block scores ROC-AUC 0.506 and
0.507: exactly chance, as a hash should.

**Why it generalises.** Placeholder features get written early and survive
into production because nothing fails when they are wrong. With few unique
sequences they can act as an identity key.

**Check for it:** ablate every feature block. Anything you cannot explain
mechanistically should earn its place.
Script: `scripts/ablate_hash_embedding.py`.

---

## 5. Richer sequence representations did not help

**What happened.** With 3,338 CNS sequences and a cross-assay split, k-mer
representations were tested against 17 hand-designed design features.

| Representation | Features | MOE → Hagedorn | Hagedorn → MOE |
|---|---|---|---|
| Hand features | 17 | **0.828** | **0.783** |
| 3-mer counts | 64 | 0.791 | 0.772 |
| 3-mer + 4-mer | 320 | 0.792 | 0.767 |
| k-mer + hand | 337 | 0.827 | 0.779 |

**Seventeen hand features beat 320 k-mer features in both directions**, and
adding k-mers to them changes nothing.

**Why it matters.** The reflex when a model underperforms is a bigger
representation. Here sequence content beyond the design features carries
nothing, which is an argument for better labels rather than better encoders.
Script: `scripts/evaluate_sequence_representations.py`.

---

## 6. Ranking and calibration come apart at low base rates

**What happened.** Every model is a `CalibratedClassifierCV`, and calibration
was never reported. Measured:

| Endpoint | PR-AUC vs baselines | Brier skill |
|---|---|---|
| cns toxicity | clears | **+0.30** |
| offtarget gene silencing | clears | +0.03 |
| hepatocellular injury | clears | **−0.01** |
| injection site reaction | fails | −2.37 |

**Hepatocellular injury beats every baseline on ranking while its
probabilities are no better than quoting the base rate.** Unsupported
endpoints over-predict badly: injection site reaction averages 0.029 against
a 0.007 observed rate.

**Why it generalises.** At base rates of 0.2–13%, a model can order
candidates usefully and still be useless as a probability. A safety report
quoting "30% chance" relies on calibration, not discrimination.

**Check for it:** report Brier skill against the base rate alongside PR-AUC,
and label any uncalibrated score as a ranking.
Script: `scripts/evaluate_calibration.py`.

---

## 7. Internal validation did not predict external transfer

**What happened.** Hepatocellular injury clears its source-only baseline,
beats GC and length baselines, and survives a within-source check — the full
internal battery. On an independent screen of 149 LNA gapmers
([NAR 2022;50(13):7224](https://doi.org/10.1093/nar/gkac562)) it scored
**ROC-AUC 0.40**, ranking the 5 hepatotoxic sequences at 23, 76, 87, 120 and
139 of 149.

With 5 positives the interval spans most of the range, so this does not prove
the model wrong. It does mean passing every internal check bought no
demonstrated external validity.

Contrast CNS neurotoxicity, which was validated by training on one assay and
testing on another: 0.83 and 0.78, both directions, 3,338 sequences.

**Why it matters.** The difference between the two endpoints is not the
modelling. It is that one had a second independent dataset and the other did
not.
Script: `scripts/evaluate_hepatotox_external.py`.

---

## 8. A validated endpoint can still be someone else's published rule

CNS neurotoxicity transfers across assays, which is the strongest result
here. But the two features Hagedorn et al. name — G count and the distance
from the 3' end to the nearest G — score **the same or better** than the full
feature set in both directions (0.839 vs 0.828, 0.797 vs 0.783).

The endpoint is a reimplementation of a published rule, validated on an
independent assay. That is worth having and worth saying. The evidence string
the app displays says it.

---

## 9. Provenance drift between published and ingested data

The project trains on 772 rows attributed to Papargyri 2020, and the
published supplementary workbook for that study holds 772 rows. **Only 12 of
74 unique training sequences match the published sequences.** The workbook
also carries per-position LNA/DNA architecture that the pipeline never
ingested, because the file is a legacy `.xls` despite its `.xlsx` extension
and nothing had opened it.

**Why it generalises.** Attribution in a curated dataset is not the same as
derivation from the cited source. Check that what you ingested matches what
the citation contains.
Script: `scripts/ingest_papargyri_design.py`.

---

## 10. 49x more data made the model worse

**What happened.** ASO Atlas is 190,927 patent-extracted measurements over
165,555 unique sequences and 343 genes, with chemistry annotation. The
shipped ASO activity model is trained on McQuisten 2007: 3,867 sequences,
Spearman 0.26 under gene-grouped CV. Swapping one for the other should have
been the easiest win in the project.

| Setup | Spearman |
|---|---|
| McQuisten 2007, 3,867 sequences (shipped) | **0.26** |
| ASO Atlas pooled, 161,832 pairs | 0.18 |
| ASO Atlas + chemistry features | 0.16 |
| ASO Atlas + cell line, dose, transfection as features | 0.15 |
| ASO Atlas, transferred to McQuisten | **-0.06** |

More data scored lower, adding chemistry made it worse, and the model
transferred at *below chance* to the study the shipped model handles.

**Why.** The readout is not commensurable. ASO Atlas pools percent
inhibition across 81 cell lines, many dosages and several transfection
methods, extracted from 417 patents spanning 2001-2025. Restricting to one
experimental context recovers it:

| Context | Sequences | Spearman |
|---|---|---|
| Pooled, 81 cell lines | 161,832 | 0.18 |
| A431 only | 29,926 | 0.26 |
| HepG2 only | 19,400 | 0.15 |
| **A431 at a single dose** | 12,435 | **0.29** |

A single cell line at a single dose, on 8% of the data, beats the whole
corpus and beats the shipped model.

**The part worth noting.** Adding cell line, dose and transfection as
features made things *worse*, not better, and context alone predicts nothing
(0.01). The heterogeneity is not an offset that a covariate can absorb: the
relationship between sequence and potency differs between contexts. Fixing
it needs per-context models or explicit interactions, not more features.

**Why it generalises.** Aggregated public corpora are built by pooling
whatever was measured, and the pooling is invisible by the time it reaches a
single column. Before treating a large dataset as one population, check
whether a single experimental slice does better than the whole.

**How it was resolved.** Not by discarding the data. One model per
experimental context, each gated against its own gene-mean and GC baselines,
reaches 0.33 in the best context against 0.26 for the single-study model.
Five of seven contexts cleared the gate and two did not. A corpus that is
useless pooled can be useful partitioned, so the finding is about the
aggregation rather than the data.

Scripts: `scripts/evaluate_aso_atlas.py`, `scripts/train_context_potency.py`.
The data is Apache-2.0 at https://github.com/barneyhill/aso_atlas.

---

## 11. The mechanistic feature with the best story carried no signal

**What happened.** Gapmer hepatotoxicity has a published mechanism: Burel et
al. and Kasuya et al. showed it is RNase H1 cleaving unintended transcripts
the oligo is complementary to, especially long pre-mRNAs. That makes
off-target complementarity burden computable at design time, and it was the
most promising feature this project could add.

It was computed properly. Three reverse-complement probes per oligo (full
length, the central DNA gap RNase H1 requires, and the seed), matched in one
Aho-Corasick pass over the whole Ensembl human cDNA set: 465,769
transcripts, 1,276 Mb. The counts have real spread, gap hits ranging 0 to
10,851 with a median of 914.

| Feature set | Features | PR-AUC |
|---|---|---|
| Design only | 39 | **0.274** |
| Design + off-target | 45 | 0.222 |
| Design + 6 random columns | 45 | 0.244 |
| Off-target only | 6 | **0.125** |

Base rate is 0.126. **Off-target burden alone scores the base rate**, and
adding it costs more than adding six columns of noise.

**The control is the point.** An earlier version of this test used the 26
cached target transcripts as a background, 0.12 Mb against 1,280 Mb, and
produced a similar null. That null meant nothing: at that size the counts
are noise. Equally, without the random-feature control the drop would have
looked like evidence the feature is actively harmful, when most of it is the
dimensionality penalty of adding any six columns to 19 positive sequences.

**Why it probably failed, and what would test it.** The mechanism implicates
long **pre-mRNA**, and a cDNA set contains no introns. So the burden
measured here is against mature mRNA, which is the wrong compartment for the
published mechanism. This is the same limitation that caps target-site
accessibility coverage at 34%: both need genomic sequence, not cDNA.

So this is not "the mechanism is wrong". It is "the mechanism was measured
in the wrong compartment, and the right measurement needs intronic
sequence". Worth stating precisely, because the vaguer version would
discourage the experiment that could still work.

Scripts: `scripts/offtarget_transcriptome.py`,
`scripts/evaluate_offtarget_hepatotox.py`.

---

## 12. Accessibility helps, but only where it can be computed

**What happened.** Target-site accessibility, the probability that a binding
site is unpaired, computed with ViennaRNA's local partition function
(`pfl_fold_up`).

| Feature set | Spearman |
|---|---|
| Design only | 0.100 |
| Design + accessibility | **0.157** |
| Accessibility only | 0.020 |

A gain of +0.057, and notably accessibility predicts almost nothing on its
own. It is complementary rather than independent, which is what the
mechanism implies: a useful oligo needs both good design and an open site.

**The catch.** Only 34% of ASO Atlas oligos locate in the canonical cDNA.
Orientation is not the problem, since the reverse complement matches and the
forward sequence does not. Gapmers frequently target pre-mRNA, and cDNA has
no introns.

Within a single assay the result is unsettled: A431 gains +0.051, HepG2
loses 0.012, on 4 to 5 genes each. Too few genes to call, which is itself
worth recording rather than reporting the pooled number alone.

Script: `scripts/evaluate_accessibility.py`.

---

## What these add up to

Of 13 toxicity endpoints, 3 are presented as predictions and only 1 has
demonstrated external validity. A 190,927-row corpus that looked like the
obvious upgrade to the potency side made it worse, and the feature with the
best mechanistic story carried no signal at all when measured in the wrong
compartment. Every failure above was found by a cheap
check that is not standard practice: a source-only baseline, a grouped split,
an ablation, a calibration curve, an external set.

The platform's honest position is not that its models are good. It is that it
knows which of them are not.
