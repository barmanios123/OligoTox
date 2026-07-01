from __future__ import annotations
from pathlib import Path
from typing import Any, Dict
from oligotox.integrations.base import BaseIntegration, IntegrationResult

class GeneformerAdapter(BaseIntegration):
    name = "geneformer"

    def __init__(self, repo_path: str = "./external/Geneformer"):
        self.repo_path = Path(repo_path)

    def is_available(self) -> bool:
        return self.repo_path.exists()

    def predict(self, payload: Dict[str, Any]) -> IntegrationResult:
        if not self.is_available():
            return IntegrationResult(
                name=self.name,
                enabled=False,
                available=False,
                score=None,
                outputs={},
                message="Geneformer repo not present at ./external/Geneformer",
            )

        return IntegrationResult(
            name=self.name,
            enabled=False,
            available=True,
            score=None,
            outputs={
                "supported": False,
                "role": "context_layer",
                "notes": "Geneformer should be used later for cell-context consequence interpretation, not primary oligo sequence scoring."
            },
            message="Geneformer detected but not activated in scoring path.",
        )
