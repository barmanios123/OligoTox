"""OSWG-aligned safety reporting for oligonucleotides"""
from dataclasses import dataclass
from typing import Dict, List, Optional
from datetime import datetime
import json

@dataclass
class OSWGSafetyReport:
    """FDA/OSWG-aligned safety assessment report."""
    compound_id: str
    sequence: str
    assessment_date: str
    toxicity_predictions: Dict[str, float]
    pkpd_summary: Dict
    efficacy_score: Optional[float] = None
    
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
        report = {
            "compound_id": self.compound_id,
            "sequence": self.sequence,
            "assessment_date": self.assessment_date,
            "safety_profile": {
                "toxicity_endpoints": self.toxicity_predictions,
                "overall_risk": self._calculate_risk_score(),
                "classification": self._classify_safety()
            },
            "pkpd_characteristics": self.pkpd_summary,
            "regulatory_notes": self._generate_regulatory_notes()
        }
        if self.efficacy_score is not None:
            report["efficacy"] = {"predicted_score": self.efficacy_score}
        return json.dumps(report, indent=2)
    
    def _calculate_risk_score(self) -> float:
        """Calculate overall safety risk score."""
        weights = {
            'hepatocellular_injury': 0.25,
            'complement_activation': 0.20,
            'thrombocytopenia': 0.20,
            'renal_accumulation': 0.15,
            'systemic_inflammation': 0.20
        }
        score = sum(self.toxicity_predictions.get(k, 0) * v for k, v in weights.items())
        return round(score, 3)
    
    def _classify_safety(self) -> str:
        """Classify safety profile."""
        risk = self._calculate_risk_score()
        if risk < 0.2:
            return "LOW_RISK"
        elif risk < 0.5:
            return "MODERATE_RISK"
        else:
            return "HIGH_RISK"
    
    def _generate_regulatory_notes(self) -> List[str]:
        """Generate regulatory recommendations."""
        notes = []
        
        # Check hepatotoxicity
        if self.toxicity_predictions.get('hepatocellular_injury', 0) > 0.5:
            notes.append("Liver function monitoring recommended in clinical trials")
        
        # Check immune activation
        if self.toxicity_predictions.get('complement_activation', 0) > 0.4:
            notes.append("Immunogenicity assessment required per OSWG guidelines")
        
        # Check renal
        if self.toxicity_predictions.get('renal_accumulation', 0) > 0.5:
            notes.append("Renal function monitoring advised")
        
        return notes
    
    def _generate_pdf(self) -> str:
        """Generate PDF report (placeholder)."""
        return f"PDF report for {self.compound_id} would be generated here"
    
    def _generate_xml(self) -> str:
        """Generate XML report for regulatory submission."""
        return f'<oswg_report compound="{self.compound_id}"></oswg_report>'
