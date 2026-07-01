"""OSWG (Oligonucleotide Safety Working Group) regulatory mapping."""
from typing import Dict, List, Optional
import json

class OSWGSafetyReport:
    """Generate OSWG-compliant safety reports for regulatory submission."""
    
    def __init__(self,
                 compound_id: str,
                 sequence: str,
                 modality: str = "ASO",
                 conjugate: Optional[str] = None):
        """Initialize OSWG safety report.
        
        Args:
            compound_id: Unique compound identifier
            sequence: Oligonucleotide sequence
            modality: Oligonucleotide modality (ASO, siRNA, GalNAc, etc.)
            conjugate: Conjugate chemistry if applicable
        """
        self.compound_id = compound_id
        self.sequence = sequence
        self.modality = modality
        self.conjugate = conjugate
        self.toxicity_data = {}
        self.pkpd_data = {}
    
    def add_toxicity_data(self, endpoint: str, value: float, unit: str = ""):
        """Add toxicity endpoint data."""
        self.toxicity_data[endpoint] = {"value": value, "unit": unit}
    
    def add_pkpd_data(self, metric: str, value: float, unit: str = ""):
        """Add PK/PD metric data."""
        self.pkpd_data[metric] = {"value": value, "unit": unit}
    
    def generate_report(self) -> Dict:
        """Generate complete OSWG safety report."""
        report = {
            "compound_id": self.compound_id,
            "sequence": self.sequence,
            "modality": self.modality,
            "conjugate": self.conjugate,
            "toxicity_profile": self.toxicity_data,
            "pkpd_metrics": self.pkpd_data,
            "regulatory_classification": self._classify_risk(),
            "recommendations": self._generate_recommendations()
        }
        return report
    
    def _classify_risk(self) -> str:
        """Classify overall risk profile."""
        # Simplified risk classification
        if not self.toxicity_data:
            return "insufficient_data"
        
        hepato_risk = self.toxicity_data.get("hepatocellular_injury", {}).get("value", 0)
        if hepato_risk > 0.5:
            return "high_risk"
        elif hepato_risk > 0.2:
            return "moderate_risk"
        else:
            return "low_risk"
    
    def _generate_recommendations(self) -> List[str]:
        """Generate regulatory recommendations."""
        recommendations = []
        
        risk = self._classify_risk()
        if risk == "high_risk":
            recommendations.append("Enhanced hepatotoxicity monitoring recommended")
            recommendations.append("Consider dose reduction or alternative design")
        elif risk == "moderate_risk":
            recommendations.append("Regular liver function monitoring advised")
        
        if self.modality == "GalNAc":
            recommendations.append("GalNAc-specific biodistribution studies recommended")
        
        return recommendations
    
    def export_json(self, filepath: str):
        """Export report as JSON."""
        report = self.generate_report()
        with open(filepath, 'w') as f:
            json.dump(report, f, indent=2)
    
    def export_xml(self) -> str:
        """Export report as XML for regulatory submission."""
        report = self.generate_report()
        xml = f'<oswg_report compound="{self.compound_id}">'
        xml += f'<sequence>{self.sequence}</sequence>'
        xml += f'<modality>{self.modality}</modality>'
        xml += '</oswg_report>'
        return xml
