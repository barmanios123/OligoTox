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

**Status:** not released. Needs an email to the authors.

A drafted request is kept locally at `drafts/asoptimizer_request_email.md`,
which is gitignored. The points it has to make are the ones above: the
endpoint passes every internal check, fails externally, and is limited by
labelled data rather than method.

---

## 2. ASO Atlas gapmer corpus

**What:** 188,521 RNase H-mediated gapmer sequences extracted from 417 USPTO
patents, described in a [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.10.29.685292v1).

**Why:** efficacy-focused rather than toxicity, but it would take the potency
models from 2,838 sequences to a corpus comparable with anything proprietary.
The siRNA potency model already transfers between independent studies
(Spearman 0.46 and 0.49), so more data has a clear path to improving it.

**Status:** preprint, availability unconfirmed. An automated read of the
bioRxiv landing page returned HTTP 429, so the data availability statement
has not been checked. Before writing to the authors, open the preprint and
look for a Zenodo DOI, GitHub repository or Hugging Face dataset; the
OligoAI platform they describe may host it directly.

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

---

## 5. Redistribution terms for the OligoGym-derived datasets (open question)

**What:** the six processed datasets in `data/oligogym/`, 1.4 MB, which the
CNS neurotoxicity and potency validation depend on.

**Checked, and resolved for now.** The Hugging Face mirror of OligoGym
states CC BY 4.0, which permits redistribution with attribution. The
repository was redistributing the files without attributing them, which is
now fixed: see `data/oligogym/ATTRIBUTION.md`.

**Residual uncertainty:** the CC BY 4.0 label covers the OligoGym packaging.
Upstream rights in each underlying study were not independently verified and
no OligoGym LICENSE file was located. If an upstream owner's terms are
narrower, replace the files with a fetch script and checksums.

**If the terms do not permit redistribution:** replace the files with a
fetch script and checksums, so the numbers stay reproducible for anyone who
obtains the data themselves.

---

## A note on repository history

`data/raw` was tracked and public until commit 437556e: 114 files, 174 MB,
including 51 PDFs that are mostly publisher copyright, and two third-party
UniTox CSVs. It is removed from the current branch and gitignored, but those
objects remain reachable in history.

Purging them was considered and deliberately not done. The repository had
been public for a matter of hours, a force push breaks every clone, and
GitHub retains unreferenced objects for a period regardless, so a rewrite
would have been disruptive and incomplete. Recorded here so the decision is
visible rather than looking like an oversight.
