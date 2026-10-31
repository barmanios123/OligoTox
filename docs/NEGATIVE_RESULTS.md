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

## 16. Six of seven contexts gained, and one seed was doing most of the work

**What happened.** Target-site accessibility helps pooled across assays
(entry 12, +0.019 on 282 genes). The shipped potency models are per-assay, so
the question is whether it improves them individually. Seven contexts, each
restricted to rows whose binding site actually locates, nothing zero-filled.

The first pass looked emphatic: six of seven gained, median about +0.05, two
of them above +0.08. Written up at that point it would have read as a clear
result.

**Then the same runs were repeated with different row samples.** Three seeds,
all seven contexts, same code:

| Context | Genes | seed 42 | seed 7 | seed 1234 | Alone |
|---|---|---|---|---|---|
| a431_2000_gymnosis | 9 | +0.104 | +0.069 | +0.075 | always + |
| a431_4000_gymnosis | 12 | +0.050 | +0.045 | +0.085 | always + |
| hepg2_4500_electroporat | 7 | +0.087 | +0.075 | +0.081 | always - |
| hepg2_1000_electroporat | 5 | +0.012 | +0.021 | +0.002 | ~0 |
| t24_anydose_other | 59 | +0.018 | +0.006 | +0.000 | ~0 |
| a549_anydose_other | 54 | +0.082 | +0.031 | +0.011 | unstable |
| shsy5y_4000_electroporat | 6 | -0.064 | -0.047 | -0.035 | stable - |

**Seed 42 produced the most favourable number in five of the seven.** That is
not a coincidence to shrug at. The first table was the luckiest draw of the
three, and it was the one written down. T24 collapses from +0.018 to exactly
zero. A549 swings sevenfold, and its own design baseline swings from -0.058
to +0.082, so the thing the gain is measured against is no steadier than the
gain.

**What survives.** Three contexts hold across every seed. The two A431
gymnosis contexts gain +0.045 to +0.104 *and* accessibility ranks positively
on its own every time. HepG2 at 4500 gains just as steadily, +0.075 to
+0.087, but accessibility is negative alone in all three runs, so there it is
correcting the design block rather than carrying signal. Those are different
strengths of evidence and should not be reported as one.

SH-SY5Y is negative at every seed, with accessibility strongly negative alone
(-0.233, -0.289, -0.342). That is reproducible harm, not instability.

**Two mechanisms proposed, both wrong.** First: accessibility should matter
under gymnosis, where an oligo finds its target unaided, and not under
electroporation, which forces entry. HepG2 at 4500 is electroporation and
gains as reliably as either gymnosis context. Second, after that: the
small-gene contexts would prove unstable and the well-powered ones would
hold. The opposite happened. HepG2 4500 on 7 genes is the steadiest row in
the table, while T24 and A549 on 54 and 59 genes are the emptiest. Gene count
did not predict stability, and neither did delivery method.

**Why it generalises.** A single cross-validated number is one draw from a
distribution nobody looked at. Repeating a cheap evaluation under a different
sample seed costs minutes and is the difference between a finding and a
coincidence. The failure mode is specific and seductive: run it once, get a
good number, write it down. Everything above was produced by code that was
already correct. The only thing that changed was running it three times.

A second rule follows. Four of seven contexts were stability-tested first,
chosen because they were the interesting ones. That is the same selection
problem one level up, and testing the remaining three is what revealed HepG2
4500 as the steadiest result rather than the throwaway it was assumed to be.
Scrutiny has to be applied to the whole table or it manufactures the answer
it was looking for.

Script: `scripts/evaluate_context_accessibility.py --only <context> --seed <n>`.

---

## 17. Every unmeasured row was being trained as a negative result

Adding a new source of liver measurements moved the score of all twelve
endpoints, including the ones with no liver data at all. Kupffer cell
activation went from 0.478 to 0.606 PR-AUC and injection site reaction from
0.209 to 0.317, and both were promoted by the evidence gate, on the strength
of 444 rows that are blank for both of them.

The cause was one line in the trainer:

```python
y = pd.to_numeric(Y[endpoint], errors="coerce").fillna(0).astype(int).values
```

A row that nobody assayed for thrombocytopenia became a row asserting that
thrombocytopenia did not occur. Thrombocytopenia was trained on 4 real
positives and 2,030 invented negatives, and reported `samples=2034` as
though that were its evidence. The same `fillna(0)` was in
`evaluate_baselines.py`, `evaluate_calibration.py` and
`endpoint_requirements.py`, so the baselines a model is gated against were
computed on the same fabricated labels.

What makes this one worth writing down is that it is invisible while the
corpus is static. Every endpoint carried the same padding, so the numbers
were wrong but stable, and nothing in the loop's health checks looks at
whether a sample count is larger than the number of labels. It only surfaced
because a new source changed the padding, and then it surfaced as *good*
news: two promotions.

After the fix each endpoint sees only the rows carrying its label. The
sparse endpoints shrink from 2,058 rows to between 24 and 28, which is what
they always were.

This is entry 13's rule applied to labels instead of features. A zero is a
measurement. If nobody measured it, there is no row.

Fixed in `src/oligotox/training/train.py` and the three evaluators; pinned by
`tests/test_unmeasured_is_not_negative.py`, which asserts that each
endpoint's reported sample count equals the number of rows carrying its
label.

The promotions of `kupffer_cell_activation` and `injection_site_reaction`
are recorded in `reports/results/endpoint_history.jsonl` under the trigger
`patent-ingest`. That log is append-only, so the entry stays and the
following run records the reversal. It is not evidence of anything.

---

## 18. The estimator was mis-specified, and no output said so

Hepatocellular injury had 296 unique sequences and still lost to GC content
alone: 0.371 PR-AUC against 0.428. A 39-feature model losing to one crude
number is a modelling result, not a data one, so
`scripts/evaluate_hepatotox_model.py` tried the alternatives on the same
sequence-grouped folds over three seeds.

| Configuration | PR-AUC | Spread over 3 seeds |
|---|---|---|
| shipped: `CalibratedClassifierCV(cv=2)`, unscaled | 0.394 | 0.041 |
| balanced logistic, scaled | 0.662 | 0.021 |
| gradient boosting, shallow | 0.674 | 0.011 |
| gradient boosting, shallow, no source-shaped features | **0.674** | **0.008** |

Decomposing the gap found two independent defects in three lines of code:

| Change | PR-AUC |
|---|---|
| shipped | 0.371 |
| + scale the features | 0.437 |
| + calibrate on 5 folds instead of 2 | **0.619** |

`cv=2` meant every calibrated sub-model was fitted on half the rows, and the
endpoint with the most data is the one that suffers most from being handed
half of it. Separately, `liblinear` applies one penalty across coefficients
whose natural scales run from 0.5 (`gc_content`) to 20 (`seq_length`), so
regularisation fell unevenly across the features.

A paired bootstrap over sequences, which is the right test rather than
comparing two marginal confidence intervals, puts the whole fix at **+0.195
[+0.024, +0.342], P(better) 0.994**. Brier skill goes from **-0.054**, worse
than quoting the base rate, to **+0.128**.

**What was not done.** Gradient boosting scores higher still, 0.679 against
0.619, but the paired test puts that difference at +0.042 [-0.032, +0.100],
P(better) 0.90. It is not distinguishable, and it would give up the
calibrated probabilities the rest of the platform is built on, so the
estimator was fixed rather than replaced. Dropping all 21 source-shaped
chemistry features left gradient boosting at 0.674 with a smaller seed
spread, so none of the gain was source confounding.

**Why it survived this long.** Nothing in the loop looks at the shape of the
estimator. The model trained, predicted, and reported a number, and the
number was simply much worse than it needed to be. It only came up because
the patent ingestion made the endpoint large enough that losing to one
feature was obviously wrong rather than plausibly just a small dataset.

The fix moved `offtarget_gene_silencing` from `weak` to `validated` and
`injection_site_reaction` from `weak` to `source_confounded`. It did not
rescue hepatocellular injury, which still fails to clear the simple
baseline (CI low 0.343 against 0.436) and still does not transfer to the
NAR screen (ROC-AUC 0.379 [0.236-0.547], where 0.5 is chance).

Pinned by `tests/test_estimator_shape.py`.

---

## 19. The models were briefly unloadable by anything but the trainer

Moving the estimator into a class broke every shipped artifact, silently, in
a way no test covered and no training run could notice.

The trainer runs as `python -m oligotox.training.train`, so its module is
`__main__`. A class defined beside the training loop pickles as
`__main__.AdaptiveCalibrated`, and joblib resolves that by import path when
loading. Training succeeded. Every `.joblib` it wrote raised
`AttributeError` in the app, in the scorer and in the external check.

It was caught by running the external hepatotoxicity check, which is the
only step that loads a model in a fresh process. The fix is to define the
estimator in `src/oligotox/models/estimator.py`, which gives the pickle a
stable import path.

The same move also lost `feature_names_in_`, which sklearn sets inside
`_validate_data` and a hand-written estimator never calls. The app's How it
works page lists each model's real inputs by reading that attribute off the
model, so every toxicity model would have appeared to consume nothing at
all. Both are now asserted by `tests/test_estimator_shape.py`, including a
test that opens every shipped model.

The lesson is narrow and worth having: a model that trains is not a model
that loads, and only loading it in another process tells you which.

---

## 20. A published motif does not replicate in the second source

Kasuya et al. (NAR 2022) report that TGC and TCC counts track
hepatotoxicity in LNA gapmers, and the Papargyri rows in this corpus agree
emphatically. It was the obvious feature to add once hepatocellular injury
had a real corpus behind it. Until the patent ingestion there was no second
in-vivo source to check it against.

Positive-minus-negative mean motif count, resampling whole sequences:

| Motif | Patents (369 seq) | Papargyri (12 seq) | oligotox_main (24 seq) | Pooled |
|---|---|---|---|---|
| TGC | **+0.002** | **-0.420** | -0.214 | -0.213 [-0.316, +0.015] |
| TCC | **-0.026** | **+0.425** | +0.229 | +0.145 [-0.051, +0.275] |

At 256 patent sequences TGC read +0.039 and TCC +0.083, both weakly
positive. At 369 both have collapsed to zero, and TCC has crossed into
disagreeing with Papargyri outright.

**TGC changes sign between sources.** It is higher in the toxic gapmers in
the patents and markedly lower in the toxic ones in Papargyri. The pooled
figure of -0.232 is Papargyri's 12 sequences outvoting 256, because the
pooling is over rows and Papargyri contributes 579 of them.

**TCC no longer agrees in direction either.** It was positive in all three
sources on the smaller corpus and is now slightly negative in the patents,
with every source's own interval spanning zero. The pooled figure no longer
excludes zero at all.

So the honest reading is that the motif association is a Papargyri
phenomenon that the patent corpus does not confirm, and TCC is at best a
weak effect the data here cannot resolve. Neither was added as a feature.
Adding TGC in particular would have taught the model to recognise which
source a sequence came from, by way of a mechanism with a citation attached
to it, which is the most persuasive possible form of entry 13.

The check is `scripts/evaluate_motif_replication.py`, and it is cheap: it
is the same number computed once per source instead of once.

---

## 21. One shared random generator made tiers depend on unrelated data

`evaluate_baselines.py` built a single `np.random.default_rng(42)` before
its loop and handed it to every endpoint in turn. Each endpoint's bootstrap
therefore drew from wherever the endpoints before it had left the stream.

Adding 1,008 hepatocellular rows from a second batch of patents moved
`offtarget_gene_silencing`'s lower bound from 0.209 to 0.188. Its model did
not change. Its data did not change. Its PR-AUC was identical to three
decimal places at 0.641. The endpoint was demoted from `validated` to
`weak` and stopped being presented, because that bound sits 0.004 from the
baseline it is compared against and the random draws had shifted underneath
it.

The seed was fixed, so every run was reproducible and the defect was
invisible in exactly the way a fixed seed is supposed to prevent. What was
not fixed was the *stream position*, which depends on how much work the
loop did before reaching each endpoint.

Each endpoint now seeds its own generator from its own name,
`default_rng([42, crc32(endpoint)])`, so its interval is reproducible and
independent of what else is in the corpus. The same fix was applied to the
two scripts written the same day, which had copied the pattern.

The wider point is about the gate rather than the generator. An endpoint
whose tier can be flipped by 0.004 is not really being measured; it is
being rounded. `offtarget_gene_silencing` has gone validated, weak,
validated, weak across four runs this session. That instability is the
finding, and it is more informative than whichever side of the line it
happens to land on.

Pinned by `tests/test_bootstrap_independence.py`.

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
verdict change. The one after it was found because a new source made a
long-standing bug move, and it announced itself as two promotions rather
than as a failure.

Two of the 13 are presented as predictions. The last five entries above all
came out of one afternoon's work on a single endpoint: adding data to it
exposed a labelling defect that had been inflating every endpoint, fixing
that exposed an estimator that had been throwing away half its training
data, and fixing that exposed models that could not be loaded outside the
process that wrote them. The one feature worth adding, a motif with a
published mechanism behind it, turned out to change sign between sources.
None of them was visible in any number the platform reported.

The platform's honest position is not that its models are good. It is that it
knows which of them are not.
