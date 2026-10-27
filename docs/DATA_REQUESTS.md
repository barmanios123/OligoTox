# Outstanding data requests

What the platform needs, why, and the state of each ask. The evidence for
each is reproducible from this repo, so these are requests with numbers
attached rather than general interest.

---

## 1. ASOptimizer rodent ALT/AST curation (highest value)

**What:** the ~6,000 curated rodent ALT/AST experiments the ASOptimizer
authors describe in [NAR 2025;53(W1):W39](https://doi.org/10.1093/nar/gkaf379),
listed as a future update rather than a current release.

**Why it matters here:** hepatocellular injury is this platform's most-used
endpoint and rests on 19 positive sequences from a single curation. It clears
every internal check and then scores ROC-AUC 0.40 on an independent screen of
149 gapmers. It is not a modelling problem; there is no second dataset. 6,000
labelled experiments is roughly a thousandfold increase over the external set
currently available.

**Status:** not released. Needs an email.

### Draft

> Subject: Request for the curated rodent ALT/AST dataset described in the ASOptimizer paper
>
> Dear Dr [Name],
>
> I'm a pharmacology/toxicology researcher working on sequence-based
> prediction of oligonucleotide toxicity. Your ASOptimizer paper mentions a
> curated set of approximately 6,000 publicly available rodent ALT and AST
> experiments, listed as a planned update. I'm writing to ask whether that
> set is available, or whether you'd be open to sharing it.
>
> The reason I'm asking is specific. I maintain an open evaluation of
> oligonucleotide toxicity models (github.com/barmanios123/OligoTox, MIT
> licensed). Its hepatocellular injury endpoint rests on 19 positive
> sequences from a single curated source. It passes every internal check I
> can construct — it beats a data-source-only baseline, beats GC and length
> baselines, and survives a within-source split — and then scores ROC-AUC
> 0.40 on the independent 149-gapmer screen from Nucleic Acids Research
> 2022;50(13):7224, ranking the five hepatotoxic sequences at 23, 76, 87,
> 120 and 139 of 149. Its probabilities are also uncalibrated, with a Brier
> skill of −0.01 against the base rate.
>
> I report all of that openly, including the negative results, because I
> think the field under-reports them. But the conclusion is that the endpoint
> is limited by labelled data rather than by method, and your curation is the
> largest I'm aware of.
>
> I'd be glad to share back anything useful: the evaluation harness is
> modality-agnostic and already reports source-confounding baselines,
> sequence-grouped splits, calibration and external transfer per endpoint. If
> it would be helpful to run it against your models, I'd be happy to.
>
> Thank you for considering it.
>
> Best regards,
> Beshoy Armanios, PhD, MPH

---

## 2. ASO Atlas gapmer corpus

**What:** 188,521 RNase H-mediated gapmer sequences extracted from 417 USPTO
patents, described in a [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.10.29.685292v1).

**Why:** efficacy-focused rather than toxicity, but it would take the potency
models from 2,838 sequences to a corpus comparable with anything proprietary.
The siRNA potency model already transfers between independent studies
(Spearman 0.46 and 0.49), so more data has a clear path to improving it.

**Status:** preprint. Check whether the corpus is released with it before
asking.

---

## 3. Open TG-GATEs / DrugMatrix liver data

**What:** liver histopathology and clinical chemistry for small molecules.

**Why:** not oligonucleotides, so not directly usable as labels. The use is a
hepatotoxicity prior — a model of what liver injury looks like in rodent
clinical chemistry, onto which oligo-specific features could be layered. Only
worth doing once there are enough oligo labels to fine-tune against, so it
ranks behind request 1.

**Status:** publicly downloadable. Not blocked, just not yet useful.

---

## 4. Measured exposure data for the PK/PD simulator

**What:** liver concentration-time data for GalNAc-siRNA or PS-ASO.

**Why:** the GalNAc PK/PD simulator in `src/oligotox/pkpd/` is mechanistic
and has never been fitted or validated against measured exposure. It is
currently the only component making quantitative claims with no evaluation
behind it at all.

**Status:** no source identified. Until there is one, the simulator should
keep saying it is unfitted rather than acquiring a second unvalidated layer.
