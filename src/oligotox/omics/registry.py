from __future__ import annotations

from oligotox.omics.schemas import OmicsModalitySpec


OMICS_REGISTRY: dict[str, OmicsModalitySpec] = {
    "genomics": OmicsModalitySpec(
        name="genomics",
        file_name="genomics_features.csv",
        prefix="g_",
        reduce_method="variance_topk",
        max_features=50,
    ),
    "epigenomics": OmicsModalitySpec(
        name="epigenomics",
        file_name="epigenomics_features.csv",
        prefix="epi_",
        reduce_method="variance_topk",
        max_features=50,
    ),
    "transcriptomics": OmicsModalitySpec(
        name="transcriptomics",
        file_name="transcriptomics_features.csv",
        prefix="tx_",
        reduce_method="variance_topk",
        max_features=50,
    ),
    "proteomics": OmicsModalitySpec(
        name="proteomics",
        file_name="proteomics_features.csv",
        prefix="prot_",
        reduce_method="variance_topk",
        max_features=50,
    ),
    "metabolomics": OmicsModalitySpec(
        name="metabolomics",
        file_name="metabolomics_features.csv",
        prefix="met_",
        reduce_method="variance_topk",
        max_features=50,
    ),
    "lipidomics": OmicsModalitySpec(
        name="lipidomics",
        file_name="lipidomics_features.csv",
        prefix="lip_",
        reduce_method="variance_topk",
        max_features=50,
    ),
    "cytokines": OmicsModalitySpec(
        name="cytokines",
        file_name="cytokines_features.csv",
        prefix="cyto_",
        reduce_method="passthrough",
        max_features=None,
    ),
    "complement": OmicsModalitySpec(
        name="complement",
        file_name="complement_features.csv",
        prefix="comp_",
        reduce_method="passthrough",
        max_features=None,
    ),
    "microbiome": OmicsModalitySpec(
        name="microbiome",
        file_name="microbiome_features.csv",
        prefix="mb_",
        reduce_method="variance_topk",
        max_features=50,
    ),
    "clinical": OmicsModalitySpec(
        name="clinical",
        file_name="clinical_features.csv",
        prefix="clin_",
        reduce_method="passthrough",
        max_features=None,
    ),
}
