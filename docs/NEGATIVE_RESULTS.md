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

## 9. A provenance discrepancy that was my own measurement error

**What I reported.** That only 12 of 74 unique training sequences matched the
published Papargyri supplementary sequences, though both held 772 rows, and
that attribution in a curated dataset is not the same as derivation from the
cited source.

**What was actually true.** The workbook's 772 gapmers collapse to **12
unique base sequences**, because they are 12 sequences each carrying about 32
distinct LNA/DNA architectures. The training data holds the same 772 rows and
the same 12 sequences, and **all 12 match**. Joining on the gapmer identifier
matches **772 of 772**.

The 12-of-74 figure compared the workbook's sequences against every sequence
in the training set from all sources, not against the Papargyri subset. It
was a denominator error, and it described a provenance problem that does not
exist.

**What the investigation did find**, which matters more than the thing it
disproved. Chemistry is the variable Papargyri screened: 290 distinct
per-position LNA/DNA patterns over 12 sequences. Uppercasing a sequence
destroys exactly that information, so the platform sees 12 sequences where
the experiment ran 772 conditions. `lna_fraction` is recomputed from the
sequence, which cannot recover an architecture the sequence does not encode.

So the dataset that most deliberately varied chemistry contributes no
chemistry signal, and the fix is a join that was available all along.

**Why it generalises.** Check the denominator before reporting a mismatch,
and check whether the identifier columns join before concluding two datasets
disagree. A normalisation applied for good reasons elsewhere, uppercasing to
make RNA and DNA spellings comparable, silently discarded the experimental
variable here.

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

**The compartment was then tested, at the scale available.** Counting the
same probes against pre-mRNA (32 genomic spans, 2.12 Mb, introns included)
versus the same 32 genes' cDNA (0.14 Mb, exons only):

| Background | Median gap hits | Design | + off-target | Off-target only |
|---|---|---|---|---|
| pre-mRNA | 2 | 0.274 | 0.224 | 0.116 |
| cDNA, same genes | 0 | 0.274 | 0.236 | 0.111 |

Introns do produce more hits, as expected, and the prediction does not
improve. Switching compartment did not rescue it.

**This still does not settle it.** 2.12 Mb of pre-mRNA is not a
transcriptome; the cDNA test that produced the headline null used 1,276 Mb.
A genome-scale pre-mRNA background remains the untested version, and the
honest position is that the compartment hypothesis has been tested at small
scale and failed there, not that it has been refuted.

Scripts: `scripts/offtarget_transcriptome.py`,
`scripts/evaluate_offtarget_hepatotox.py`.

---

## 12. Fixing a coverage bias shrank the effect it was hiding

**What happened.** Target-site accessibility, P(binding site unpaired), from
ViennaRNA's local partition function.

Measured first against canonical cDNA, where only 34% of oligos locate, it
looked worthwhile: design features alone 0.100, with accessibility 0.157, a
gain of **+0.057**. That number was reported as evidence accessibility helps.

Then the coverage problem was fixed. Gapmers target pre-mRNA and a cDNA has
no introns, so the same oligos were relocated against genomic spans:

| Target sequence | Genes | Located | Design | + accessibility | Gain | Alone |
|---|---|---|---|---|---|---|
| canonical cDNA | 17 | 34% | 0.100 | 0.157 | **+0.057** | +0.020 |
| genomic span | 59 | 93% | 0.155 | 0.169 | **+0.014** | -0.023 |
| genomic span | 122 | 83% | 0.193 | 0.211 | **+0.018** | -0.004 |
| genomic span | 178 | 83% | 0.206 | 0.225 | **+0.019** | +0.032 |
| **genomic span** | **282** | **82%** | **0.229** | **0.248** | **+0.019** | **+0.062** |

**The gain nearly vanished, and then stayed put.** The correction from +0.057
to +0.014 held: against the cDNA sample the feature looked four times as
useful as it is. But the gain has not moved since, across a near fivefold
increase in genes, which is the behaviour of a small real effect rather than
noise.

**Why the first number was wrong.** The 34% that located in a cDNA were
exonic sites, which is a selected subset rather than a random one: they
differ in local composition and structure from the intronic majority. The
apparent benefit was partly that selection, and partly the instability of 17
genes under gene-grouped folds. An earlier run on 4 genes gave +0.155, which
should have been read as noise rather than a strong result.

**A claim in this entry was itself undersampled.** The line above used to
read "accessibility alone scores -0.023, below zero", offered as evidence the
feature carries nothing by itself. That number came from 59 genes, and those
59 were not a sample anyone chose: entry 15 describes the error handler that
had silently discarded 44% of the gene set. On 282 genes accessibility alone
scores **+0.062**, above the 0.05 the evaluation now requires before calling
a feature standalone signal.

So the verdict moved from "adds a little in combination, but carries no
signal on its own" to "adds signal", and it moved because the sample grew,
not because a threshold did. Both thresholds were set while the number was
+0.032 and failing. The correction here is narrower than the one it sits
inside: the coverage-bias finding stands, and a secondary claim made from the
truncated remainder did not.

**Why it generalises.** A missing-data problem that correlates with the
feature being tested inflates that feature's apparent value. The instinct is
to treat low coverage as a power problem to be apologised for; here it was a
bias problem, and fixing it reversed the conclusion rather than sharpening
it. If a feature can only be computed for a subset, establish why that
subset is a subset before trusting anything measured on it.

**What this does not settle: predicted is not measured.** The literature
monitor later surfaced Gallo et al. 2026 (iScience, PMC13651819), who
designed 24 gapmers against SARS-CoV-2 regions chosen for single-strandedness
and report a significant Pearson correlation between target SHAPE-MaP
reactivity and gapmer efficacy.

That is a claim in the opposite direction, and it is not in conflict with the
null above, because the two measure different things:

| | here | Gallo et al. |
|---|---|---|
| accessibility from | ViennaRNA prediction | SHAPE-MaP chemical probing |
| target | human pre-mRNA | SARS-CoV-2 genome |
| oligos | ~5,000 over 59 genes | ~10 assayed, of 24 designed |
| accessibility range | whatever the sites happened to be | pre-selected to be accessible |

Their result is small-n and range-restricted, since every candidate was
chosen for accessibility before efficacy was measured, so the correlation is
within an already-accessible set. It cannot carry much weight on its own. But
it does isolate the one variable this repository has never tested: whether
the null here is about accessibility, or only about *predicted* accessibility.
ViennaRNA's local partition function and a chemical probing experiment are
not interchangeable, and conflating them is the kind of substitution this
document exists to catch.

What would settle it is per-nucleotide probing data for human transcripts
these oligos target, scored the same way. That is tracked as
`shape_human_transcripts` in `data/registry/sources.json`. Until then the
honest claim is narrow: **predicted** accessibility adds +0.014 here, and
measured accessibility is untested.

Scripts: `scripts/evaluate_accessibility.py`,
`scripts/fetch_transcripts.py --genomic`.

---

## 13. Zero-filling a feature turns it into a source indicator

**What happened.** Entry 9 established that the Papargyri workbook carries
real per-position LNA/DNA architecture for 772 training rows, joining 772 of
772 on the gapmer identifier. That is the only genuine chemistry annotation
in the project, on the dataset that deliberately varied chemistry.

Adding it looked like the best single result the platform had produced:

| Feature set | Hepatocellular injury PR-AUC |
|---|---|
| Sequence features (39) | 0.274 |
| **+ chemistry architecture** | **0.564** |

A doubling. It was entirely an artifact.

**Why.** The chemistry columns exist for the 772 Papargyri rows and were
filled with zeros for the other 842. A column that is zero exactly when a row
is not Papargyri **is a source indicator**, and 193 of the 203 hepatocellular
positives are Papargyri rows: P(positive | Papargyri) = 0.250 against
P(positive | other) = 0.012.

Dropping the explicit `design_available` flag did not help, because the
zero-filled columns encode the same thing implicitly. The flag alone scores
0.201, which is how the problem was spotted.

**The clean test.** Evaluate inside Papargyri only, where chemistry varies
across 16 LNA fractions and 10 gap lengths and the source is constant:

| Feature set, within Papargyri | PR-AUC |
|---|---|
| Sequence features only | **0.541** |
| + real chemistry | 0.435 |
| Chemistry only | 0.449 |

Chemistry makes it **worse** where it is the actual experimental variable.
The entire apparent gain was provenance, not biology.

**Why it generalises.** A feature available for one subset and imputed for
the rest carries subset membership, whatever the imputation value. The usual
defences do not catch it: the folds were grouped by sequence, no label was
leaked, and the feature is mechanistically plausible. This is entry 3's
failure mode arriving through a new door, found only by holding the source
constant.

Any future feature that covers part of the dataset should be evaluated
within its covered subset before anyone looks at the pooled number.

Script: `scripts/ingest_papargyri_design.py`.

---

## 14. A threshold at exactly zero makes a verdict a coin flip

**What happened.** The audit harness decided calibration by the sign of the
Brier skill: at or below zero the probabilities are worth no more than
quoting the base rate, so the endpoint "ranks but is uncalibrated".

Hepatocellular injury sits almost exactly on that line. The committed run
recorded **-0.009** and the verdict *ranks but uncalibrated*. Re-running the
identical code on the pinned versions on a different machine gave **+0.002**
and the verdict *passes*. Nothing about the model, the data or the features
had changed. The endpoint had simply crossed zero by one hundredth, and the
published claim about it inverted.

**Why the sign test was wrong.** Zero is where the noise is. A metric that is
genuinely indistinguishable from the base rate will scatter either side of
zero across BLAS builds, library versions and platforms, so a test against
zero resolves that scatter into a confident verdict in whichever direction
the run happened to land. The question "are these probabilities better than
the base rate" needs a margin wide enough to survive the variation, not a
sign.

**The fix.** `MIN_BRIER_SKILL = 0.02`, which is the threshold
`scripts/derive_evidence_tiers.py` already used for the same decision. The
two had drifted apart: the gate that decides what the app shows demanded a
margin, while the harness that audits anyone's data accepted a sign. Both now
require 0.02, hepatocellular injury reads *ranks but uncalibrated* again, and
it does so for a reason that a rerun cannot flip.

**Why it generalises.** A verdict is only as stable as the distance between
its threshold and the noise floor. When a project defines the same threshold
in two places, the looser one is the one that will eventually be quoted.

Script: `packages/evidence-audit/src/evidence_audit/checks.py`.

---

## 15. An error handler deleted 44% of the dataset

**What happened.** The genomic fetch reported 151 of 350 target genes as
unresolvable. That number was accepted and the accessibility evaluation ran
on the 199 that remained.

The 151 were not unresolvable. APOC3 and BRCA1 were both on the list, and
both resolve on the first attempt when asked again.

**The cause.** `_get` caught five exception types in one tuple:

```python
NETWORK_ERRORS = (urllib.error.HTTPError, urllib.error.URLError,
                  json.JSONDecodeError, TimeoutError, OSError)
```

and the caller turned any of them into "Ensembl cannot place this gene",
which was written to a blacklist that later runs skip. A 404 means the symbol
does not exist. A 429, a 5xx or a timeout means the request lost and says
nothing whatsoever about the gene. One rate-limited window therefore removed
44% of the corpus permanently, and `--retry-misses` was the only way back.

**Why the first hypothesis was wrong, and worth recording.** The obvious
guess was that large genes were timing out, since the genomic request
downloads an entire span and APP is ~290kb. That is a tidy story and it is
false. The resolved genes run up to 582kb, and the "failures" include APOC3
at about 3kb. Size had no part in it. The temptation was to write the
plausible mechanism into the fix; the data refused it.

**What it cost.** Nothing was wrong with the recorded numbers, but they were
computed on an arbitrary 57% of the intended genes, and nothing in the output
said so. The coverage field read 83% located, which describes the oligos that
mapped within the genes that happened to survive, not the genes themselves. A
sample can be silently truncated upstream of the number that reports
coverage.

**The fix.** Unknown and transient are separate exceptions. Only unknown is
remembered. Transient retries with backoff that honours `Retry-After`, and
the manifest keeps the two in different keys so the distinction survives a
restart. Of the original 151, **17** are genuinely absent from Ensembl.

**Why it generalises.** An exception handler is a claim about what a failure
means. Catching several kinds together asserts they mean the same thing, and
when one of them is "try again later" and another is "this does not exist",
that assertion quietly destroys data. The blast radius is set by what the
handler writes down, and a handler that writes to a skip-list is making
permanent decisions from transient evidence.

Script: `scripts/fetch_transcripts.py`.

---

## What these add up to

Of 13 toxicity endpoints, 3 are presented as predictions and only 1 has
demonstrated external validity. A 190,927-row corpus that looked like the
obvious upgrade to the potency side made it worse, and the feature with the
best mechanistic story carried no signal at all when measured in the wrong
compartment. A second feature looked useful until a coverage bias was fixed,
which cut its apparent benefit to a quarter. Every failure above was found by a cheap
check that is not standard practice: a source-only baseline, a grouped split,
an ablation, a calibration curve, an external set. The last of them was
found by re-running a committed result on a second machine and watching the
verdict change.

The platform's honest position is not that its models are good. It is that it
knows which of them are not.
