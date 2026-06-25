# OligoTox biology explainer

OligoTox is a research tool that helps users think about safety risk for RNA and oligonucleotide medicines. It combines curated biology knowledge, known chemistry patterns, and model-based predictions to estimate which toxicity endpoints may be most important for a candidate.

## What the app is trying to do

The app does not try to say whether a drug is "safe" or "unsafe" in an absolute sense. Instead, it helps answer a more practical question: **given this sequence, chemistry, and delivery setup, what kinds of toxicity are most biologically plausible and what should be checked next?**

In simple terms, the platform looks at a candidate and asks:

- What kind of oligonucleotide is this?
- How is it delivered?
- What tissues does it likely reach?
- Which toxicity mechanisms are known to be associated with this chemistry and platform?
- Which endpoints should be watched most closely?

## What PMO means

PMO stands for **phosphorodiamidate morpholino oligomer**. PMOs are a special oligonucleotide chemistry that do **not** have the phosphorothioate (PS) backbone used in many antisense oligonucleotides.

That matters biologically because PS chemistry is associated with several toxicity pathways, including platelet effects and complement-related immune effects. PMOs often show a different pattern: renal accumulation and kidney tubular findings can still happen, but some of the classic PS-driven toxicities are much less expected.

## What goes into a score

The app takes a small set of biologically meaningful inputs, such as:

- Sequence
- Modality, for example siRNA, ASO, or PMO
- Delivery vehicle, for example GalNAc or LNP
- Route, such as subcutaneous or intravenous
- Approximate chemistry fractions, such as PS fraction and LNA fraction
- Repeat dosing status
- Target tissue and related context when available

These inputs are turned into features that represent chemistry, platform, delivery, and likely mechanism. The app then uses trained endpoint models plus rule-based biological logic to estimate risk across multiple toxicity endpoints.

## What the endpoint scores mean

Each endpoint score is a model-backed estimate of relative concern for one toxicity category. A higher value means the platform sees more evidence that the endpoint could be relevant for the candidate.

Examples of endpoints include:

- Hepatocellular injury, meaning liver cell injury
- Renal tubular injury, meaning injury in kidney tubules
- Renal accumulation, meaning the chemistry tends to collect in kidney tissue
- Hypersensitivity, meaning allergy-like or infusion-related immune reactions
- Thrombocytopenia, meaning low platelet risk
- Complement activation, meaning innate immune activation through complement pathways
- Off-target gene silencing, meaning unintended RNA knockdown effects

These values are best interpreted comparatively, not as absolute clinical probabilities. In practice, users should ask: **which endpoints rose to the top, and do those top risks match what is biologically known for this class?**

## What the D1, D2, and D3 scores mean

The app groups toxicity logic into broad mechanism buckets so the output is easier to interpret.

A simple way to think about them is:

- **D1** = sequence- and hybridization-driven effects, such as off-target silencing or unwanted engagement with RNA biology
- **D2** = delivery, exposure, and tissue-distribution effects, such as where the molecule goes and what organs receive the highest burden
- **D3** = chemistry and backbone effects, such as PS-driven platelet or complement biology and other scaffold-linked toxicity patterns

These buckets are not clinical diagnoses. They are a way of explaining *why* the app is prioritizing certain endpoints.

## How to read a result

A useful reading order is:

1. Look at the **highest endpoint risk**.
2. Look at the **top D1/D2/D3 bucket** to understand the likely mechanism.
3. Check whether the result matches known class biology.
4. Use that pattern to decide what assays or safety checks should come next.

For example:

- A PMO-like profile showing kidney and hypersensitivity signals but very low thrombocytopenia and complement is biologically coherent.
- A PS-ASO-like profile showing more platelet or coagulation concern than a PMO control is biologically coherent.
- An LNP-delivered payload showing more innate immune, cytokine, or complement signal than a matched GalNAc payload is biologically coherent.

## Why some predictions are stronger than others

Not all endpoints are equally well trained.

Some endpoints have more data and stronger signal in the training set. Others are still sparse, especially rare safety outcomes such as thrombocytopenia, complement activation, some neurotoxicity labels, and reproductive or developmental toxicity.

That means the app should be used as a **research-grade prioritization tool**, not as a final decision-maker. Its strongest use cases are:

- Candidate ranking
- Mechanism-aware interpretation
- Early assay planning
- Highlighting which toxicity questions matter most for a given chemistry

## How to think about confidence

A confidence label such as Low, Medium, or High is a guide to how stable or supported the result appears within the current system. It should not be treated as regulatory-grade certainty.

Higher confidence generally means the model and rules are seeing a more recognizable biology pattern. Lower confidence often means the candidate is unusual, the endpoint is sparse, or the available training data are still limited.

## What the app is best at right now

The platform is most helpful when it is used to compare plausible biological patterns:

- PMO versus PS-backbone oligonucleotides
- GalNAc delivery versus LNP delivery
- Renal-driven versus immune-driven versus off-target-driven risk profiles

This is where the combination of curated biology and model scoring is most informative.

## What users should not over-interpret

Users should avoid treating the numbers as exact clinical probabilities. A score does not mean a patient will definitely experience that event. It means the platform sees a stronger or weaker mechanistic reason to watch that endpoint.

The right interpretation is:

- Which endpoints stand out?
- Does that pattern make biological sense for the modality and chemistry?
- What experiment, monitoring plan, or comparator should be used next?

## One-sentence summary

OligoTox is a biology-guided ranking tool: it uses chemistry, delivery, and sequence context to estimate which toxicity mechanisms and endpoints are most worth paying attention to for an RNA or oligonucleotide candidate.
