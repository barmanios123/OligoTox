# Attribution for the datasets in this directory

These files are derived from **OligoGym**, and are redistributed here under
the terms below. The CNS neurotoxicity endpoint and both potency models are
trained and validated on them, so they ship with the repository to keep the
published numbers reproducible.

## Source

> Rotrattanadumrong & De Donno. *OligoGym: Curated Datasets and Benchmarks
> for Oligonucleotide Drug Discovery.* NeurIPS 2025, Datasets and Benchmarks
> Track.
> https://papers.nips.cc/paper_files/paper/2025/hash/4927dd556ff6e49671e3a77ac1f12487-Abstract-Datasets_and_Benchmarks_Track.html

Dataset mirror: https://huggingface.co/datasets/CollageBio/oligogym

## Licence

The Hugging Face mirror states **CC BY 4.0**, which permits redistribution
with attribution. This file is that attribution.

**Known limit of this check:** the CC BY 4.0 label covers the OligoGym
packaging. Upstream rights in each underlying study have not been
independently verified, and no OligoGym `LICENSE` file or data-use agreement
was located. If an upstream owner's terms turn out to be narrower, the
remedy is to replace these files with a fetch script and checksums, so the
numbers stay reproducible for anyone who obtains the data themselves.

## Underlying studies

Each file carries data from a published study, credited individually in
`docs/provenance/DATA_AND_MODEL_PROVENANCE.md`:

| File | Study |
|---|---|
| `MOE_Neurotox_1_processed.csv.gz` | MOE neurotoxicity in vivo panel |
| `Hagedorn_2022_1_processed.csv.gz` | Hagedorn et al., Nucleic Acid Ther. 2022, doi:10.1089/nat.2021.0071 |
| `McQuisten_2007_1_processed.csv.gz` | McQuisten & Peek, BMC Bioinformatics 2007;8:184 |
| `Ichihara_2007_{1,2}_processed.csv.gz` | Ichihara et al. 2007 |
| `Papargyri_2020_1_processed.csv.gz` | Papargyri et al. 2020 |

Cite the original studies, and OligoGym for the curation. Not this
repository.
