from __future__ import annotations
from pathlib import Path
from typing import Any, Dict
from oligotox.integrations.base import BaseIntegration, IntegrationResult

class OligoGraphAdapter(BaseIntegration):
    name = "oligograph"

    def __init__(self, repo_path: str = "./external/OligoGraph"):
        self.repo_path = Path(repo_path)

    def is_available(self) -> bool:
        return self.repo_path.exists()

    def predict(self, payload: Dict[str, Any]) -> IntegrationResult:
        if not self.is_available():
            return IntegrationResult(
                name=self.name,
                enabled=True,
                available=False,
                score=None,
                outputs={},
                message="OligoGraph repo not present at ./external/OligoGraph",
            )

        seq = str(payload.get("sequence", ""))
        modality = str(payload.get("modality", "")).lower()

        if modality != "sirna":
            return IntegrationResult(
                name=self.name,
                enabled=True,
                available=True,
                score=None,
                outputs={"supported": False, "reason": "Current OligoGraph integration is limited to siRNA flows."},
                message="Unsupported modality for OligoGraph.",
            )

        heuristic = min(0.95, max(0.05, 0.45 + (seq.count("G") + seq.count("C")) / max(len(seq), 1) * 0.2))
        return IntegrationResult(
            name=self.name,
            enabled=True,
            available=True,
            score=round(heuristic, 4),
            outputs={
                "supported": True,
                "engine": "oligograph_stub",
                "efficacy_score": round(heuristic, 4),
                "notes": "Stub adapter active. Replace with true RiNALMo/OligoGraph inference once repo + weights are installed."
            },
            message="Stub inference completed.",
        )
