from __future__ import annotations
from typing import Dict
from oligotox.integrations.oligograph.adapter import OligoGraphAdapter
from oligotox.integrations.geneformer.adapter import GeneformerAdapter

def get_integrations() -> Dict[str, object]:
    return {
        "oligograph": OligoGraphAdapter(),
        "geneformer": GeneformerAdapter(),
    }
