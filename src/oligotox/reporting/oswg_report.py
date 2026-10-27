"""OSWG-aligned safety reporting for oligonucleotides"""
from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from datetime import datetime
import json

@dataclass
class OSWGSafetyReport:
    """FDA/OSWG-aligned safety assessment report.

    `evidence_tiers` maps each endpoint to its tier from
    artifacts/config/endpoint_evidence.json. Only endpoints tiered "validated"
    reach the risk score and the regulatory notes; the rest are listed as
    unsupported. Passing nothing means no evidence information is available,
    and the report says so rather than treating every endpoint as sound.
    """
    compound_id: str
    sequence: str
    assessment_date: str
    toxicity_predictions: Dict[str, float]
    pkpd_summary: Dict
    efficacy_score: Optional[float] = None
    evidence_tiers: Optional[Dict[str, str]] = None

    def _supported_endpoints(self) -> List[str]:
        """Endpoints the evaluation backs, so safe to act on."""
        if self.evidence_tiers is None:
            return list(self.toxicity_predictions)
        return [e for e in self.toxicity_predictions if self.evidence_tiers.get(e) == "validated"]

    def _supported_predictions(self) -> Dict[str, float]:
        return {e: self.toxicity_predictions[e] for e in self._supported_endpoints()}

    def generate_report(self, format='json') -> str:
        """Generate formatted safety report.
        
        Args:
            format: 'json', 'pdf', or 'xml'
            
        Returns:
            Formatted report string
        """
        if format == 'json':
            return self._generate_json()
        elif format == 'pdf':
            return self._generate_pdf()
        else:
            return self._generate_xml()
    
    def _generate_json(self) -> str:
        """Generate JSON report."""
        supported = self._supported_endpoints()
        report = {
            "compound_id": self.compound_id,
            "sequence": self.sequence,
            "assessment_date": self.assessment_date,
            "safety_profile": {
                "toxicity_endpoints": self.toxicity_predictions,
                "overall_risk": self._calculate_risk_score(),
                "classification": self._classify_safety()
            },
            "evidence": {
                "endpoints_supporting_assessment": supported,
                "endpoints_not_supported": [
                    e for e in self.toxicity_predictions if e not in supported
                ],
                "tiers": self.evidence_tiers,
                "note": (
                    "Only endpoints validated against source, GC and length baselines "
                    "contribute to overall_risk, classification and regulatory_notes. "
                    "toxicity_endpoints lists every model output, including unsupported ones."
                    if self.evidence_tiers is not None
                    else "No evidence tiers supplied; every endpoint was treated as scoring, "
                         "which this build cannot justify. Supply evidence_tiers."
                ),
            },
            "conclusions": self._conclusions(),
            "pkpd_characteristics": self.pkpd_summary,
            "regulatory_notes": self._generate_regulatory_notes()
        }
        if self.efficacy_score is not None:
            report["efficacy"] = {"predicted_score": self.efficacy_score}
        return json.dumps(report, indent=2)
    
    # Relative severity, used only for the endpoints that reach the score. The
    # weights are normalised over whatever is actually present, so dropping an
    # unsupported endpoint doesn't cap the score below its own range.
    SEVERITY_WEIGHTS = {
        'hepatocellular_injury': 0.25,
        'complement_activation': 0.20,
        'thrombocytopenia': 0.20,
        'renal_accumulation': 0.15,
        'cytokine_release': 0.20,
    }


    # What a reader may and may not take from this report. A safety assessment
    # that states the limits of its own evidence is more use to a reviewer than
    # one that scores every endpoint with equal confidence.
    EVIDENCE_STANDARD = (
        "An endpoint supports a conclusion only if it has at least 10 positive "
        "sequences, beats a baseline using only the originating dataset, beats "
        "single-feature baselines, and retains signal within one source. "
        "Cross-validation folds are grouped by sequence."
    )

    def _conclusions(self) -> Dict[str, Any]:
        supported = self._supported_endpoints()
        unsupported = [e for e in self.toxicity_predictions if e not in supported]
        readable = [e.replace("_", " ") for e in supported]

        can = []
        cannot = []
        if supported:
            can.append(
                "Relative ranking of candidates on " + ", ".join(readable)
            )
            can.append(
                "Comparison of this candidate against others scored by the same models"
            )
        else:
            cannot.append("Any toxicity conclusion: no endpoint met the evidence standard")
        if unsupported:
            cannot.append(
                f"Any conclusion about {len(unsupported)} endpoint(s) whose models did "
                "not meet the evidence standard, however they scored"
            )
        cannot.append("An absolute probability of toxicity in any species")
        cannot.append("A prediction of clinical outcome, dose selection or margin of safety")
        cannot.append("Substitution for in vivo toxicology required by regulation")

        return {
            "evidence_standard": self.EVIDENCE_STANDARD,
            "supports": can,
            "does_not_support": cannot,
            "endpoints_meeting_standard": supported,
            "endpoints_below_standard": unsupported,
            "intended_use": (
                "Design-time triage to deprioritise candidates. Not a safety "
                "determination and not a substitute for study data."
            ),
        }

    def _calculate_risk_score(self) -> float:
        """Weighted mean risk over the endpoints the evaluation supports."""
        supported = self._supported_predictions()
        pairs = [
            (prob, self.SEVERITY_WEIGHTS[endpoint])
            for endpoint, prob in supported.items()
            if endpoint in self.SEVERITY_WEIGHTS
        ]
        total_weight = sum(w for _, w in pairs)
        if not pairs or total_weight <= 0:
            return 0.0
        return round(sum(p * w for p, w in pairs) / total_weight, 3)

    def _classify_safety(self) -> str:
        """Classify safety profile."""
        if not self._supported_endpoints():
            return "UNCLASSIFIED_NO_VALIDATED_ENDPOINT"
        risk = self._calculate_risk_score()
        if risk < 0.2:
            return "LOW_RISK"
        elif risk < 0.5:
            return "MODERATE_RISK"
        else:
            return "HIGH_RISK"

    def _generate_regulatory_notes(self) -> List[str]:
        """Recommendations, raised only from endpoints the evaluation supports.

        A monitoring recommendation driven by a model fit on a handful of
        positive sequences is worse than no recommendation, so unsupported
        endpoints raise nothing here however high they score.
        """
        supported = self._supported_predictions()
        notes = []

        if supported.get('hepatocellular_injury', 0) > 0.5:
            notes.append("Liver function monitoring recommended in clinical trials")

        if supported.get('complement_activation', 0) > 0.4:
            notes.append("Immunogenicity assessment required per OSWG guidelines")

        if supported.get('renal_accumulation', 0) > 0.5:
            notes.append("Renal function monitoring advised")

        unsupported = [e for e in self.toxicity_predictions if e not in supported]
        if unsupported:
            notes.append(
                f"{len(unsupported)} of {len(self.toxicity_predictions)} endpoint models are not "
                "supported by the held-out evaluation and raised no recommendation here, whatever "
                "they scored. See the evidence block."
            )
        if not supported:
            notes.append(
                "No endpoint passed the evaluation, so this report supports no safety conclusion."
            )

        return notes
    
    def _generate_pdf(self) -> str:
        """Generate PDF report (placeholder)."""
        return f"PDF report for {self.compound_id} would be generated here"
    
    def _generate_xml(self) -> str:
        """Generate XML report for regulatory submission."""
        return f'<oswg_report compound="{self.compound_id}"></oswg_report>'
