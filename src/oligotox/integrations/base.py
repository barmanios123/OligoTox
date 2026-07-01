from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Dict, Optional

@dataclass
class IntegrationResult:
    name: str
    enabled: bool
    available: bool
    score: Optional[float]
    outputs: Dict[str, Any]
    message: str = ""

class BaseIntegration:
    name = "base"

    def is_available(self) -> bool:
        return False

    def predict(self, payload: Dict[str, Any]) -> IntegrationResult:
        return IntegrationResult(
            name=self.name,
            enabled=False,
            available=False,
            score=None,
            outputs={},
            message="Integration not implemented.",
        )
