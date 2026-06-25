"""
Unit tests for Phase 2A PK/PD Bridge
Tests GalNAc-siRNA PK/PD simulator and toxicity risk adjustment
"""
import pytest
import numpy as np
from oligotox.pkpd import (
    GalNAcPKPDSimulator,
    compute_pkpd_enriched_features,
    adjust_toxicity_risk_with_pkpd
)


class TestGalNAcPKPDSimulator:
    """Test suite for GalNAc PK/PD 2-compartment model"""
    
    def test_simulator_initialization(self):
        """Test simulator initializes with default parameters"""
        sim = GalNAcPKPDSimulator()
        assert sim.species == "human"
        assert sim.CL_liver > 0
        assert sim.Vd_liver > 0
        assert sim.f_liver > 0
        
    def test_simulate_basic_Q4W(self):
        """Test basic simulation for Q4W dosing"""
        sim = GalNAcPKPDSimulator()
        result = sim.simulate(
            dose_mg_kg=3.0,
            schedule="Q4W"
        )
        
        # Check output structure
        assert "liver_auc_pmol_g_h" in result
        assert "liver_cmax_pmol_g" in result
        assert "liver_half_life_h" in result
        
        # Check reasonable values for 3 mg/kg Q4W
        assert 5000 < result["liver_auc_pmol_g_h"] < 50000
        assert 100 < result["liver_cmax_pmol_g"] < 2000
        assert 20 < result["liver_half_life_h"] < 100
        
    def test_dose_proportionality(self):
        """Test that PK metrics scale with dose"""
        sim = GalNAcPKPDSimulator()
        
        result_3mg = sim.simulate(dose_mg_kg=3.0, schedule="Q4W")
        result_6mg = sim.simulate(dose_mg_kg=6.0, schedule="Q4W")
        
        # AUC should approximately double with dose (allowing some nonlinearity)
        auc_ratio = result_6mg["liver_auc_pmol_g_h"] / result_3mg["liver_auc_pmol_g_h"]
        assert 1.7 < auc_ratio < 2.3
        
    def test_schedule_effects(self):
        """Test different dosing schedules"""
        sim = GalNAcPKPDSimulator()
        
        q4w = sim.simulate(dose_mg_kg=3.0, schedule="Q4W")
        q8w = sim.simulate(dose_mg_kg=3.0, schedule="Q8W")
        qw = sim.simulate(dose_mg_kg=3.0, schedule="QW")
        
        # Weekly dosing should have higher AUC than monthly
        assert qw["liver_auc_pmol_g_h"] > q4w["liver_auc_pmol_g_h"]
        
        # Q8W should have lower AUC than Q4W (less frequent dosing)
        assert q8w["liver_auc_pmol_g_h"] < q4w["liver_auc_pmol_g_h"]


class TestPKPDFeatureComputation:
    """Test PK/PD feature enrichment"""
    
    def test_feature_computation_basic(self):
        """Test basic PK/PD feature computation"""
        features = compute_pkpd_enriched_features(
            sequence="UGCUAGCUAGCUAGCU",
            modality="siRNA",
            chemistry="PS",
            dose_mg_kg=3.0,
            schedule="Q4W"
        )
        
        # Check all required features present
        assert "liver_auc_pmol_g_h" in features
        assert "liver_cmax_pmol_g" in features
        assert "exposure_score" in features
        assert "peak_knockdown_pct" in features
        
        # Check reasonable ranges
        assert 0 < features["exposure_score"] < 5.0
        assert 30 < features["peak_knockdown_pct"] < 100
        
    def test_chemistry_modifiers(self):
        """Test that chemistry affects PK/PD metrics"""
        features_ps = compute_pkpd_enriched_features(
            sequence="UGCUAGCUAGCUAGCU",
            modality="siRNA",
            chemistry="PS",
            dose_mg_kg=3.0,
            schedule="Q4W"
        )
        
        features_po = compute_pkpd_enriched_features(
            sequence="UGCUAGCUAGCUAGCU",
            modality="siRNA",
            chemistry="PO",
            dose_mg_kg=3.0,
            schedule="Q4W"
        )
        
        # PS should have higher exposure than PO (more stable)
        assert features_ps["exposure_score"] > features_po["exposure_score"]


class TestToxicityRiskAdjustment:
    """Test PK/PD bridge integration with OligoTox risk predictions"""
    
    def setup_method(self):
        """Create mock OligoTox result"""
        self.mock_result = {
            "risk_scores": {
                "hepatotoxicity": 0.65,
                "nephrotoxicity": 0.30,
                "thrombocytopenia": 0.45,
                "immune_activation": 0.20
            },
            "dominant_endpoint": "hepatotoxicity",
            "attribution": {
                "D1": 0.3,
                "D2": 0.5,
                "D3": 0.2
            }
        }
    
    def test_pkpd_adjustment_preserves_structure(self):
        """Test that PK/PD adjustment maintains result structure"""
        adjusted = adjust_toxicity_risk_with_pkpd(
            oligotox_result=self.mock_result.copy(),
            sequence="UGCUAGCUAGCUAGCU",
            modality="siRNA",
            chemistry="PS",
            dose_mg_kg=3.0,
            schedule="Q4W"
        )
        
        # Check structure preserved
        assert "risk_scores" in adjusted
        assert "pkpd_metrics" in adjusted
        assert "dominant_endpoint" in adjusted
        
        # Check PK/PD metrics added
        assert "liver_auc_pmol_g_h" in adjusted["pkpd_metrics"]
        assert "exposure_score" in adjusted["pkpd_metrics"]
        
    def test_hepatotox_adjustment_direction(self):
        """Test that hepatotox risk increases with high exposure"""
        result_low_dose = adjust_toxicity_risk_with_pkpd(
            oligotox_result=self.mock_result.copy(),
            sequence="UGCUAGCUAGCUAGCU",
            modality="siRNA",
            chemistry="PS",
            dose_mg_kg=1.0,  # Low dose
            schedule="Q12W"  # Infrequent
        )
        
        result_high_dose = adjust_toxicity_risk_with_pkpd(
            oligotox_result=self.mock_result.copy(),
            sequence="UGCUAGCUAGCUAGCU",
            modality="siRNA",
            chemistry="PS",
            dose_mg_kg=6.0,  # High dose
            schedule="QW"  # Frequent
        )
        
        # High dose/exposure should increase hepatotox risk
        original_hepatotox = self.mock_result["risk_scores"]["hepatotoxicity"]
        low_dose_hepatotox = result_low_dose["risk_scores"]["hepatotoxicity"]
        high_dose_hepatotox = result_high_dose["risk_scores"]["hepatotoxicity"]
        
        # Check adjustment directions
        assert high_dose_hepatotox >= original_hepatotox
        assert low_dose_hepatotox <= original_hepatotox
        
    def test_risk_scores_bounded(self):
        """Test that adjusted risk scores remain in [0, 1]"""
        adjusted = adjust_toxicity_risk_with_pkpd(
            oligotox_result=self.mock_result.copy(),
            sequence="UGCUAGCUAGCUAGCU",
            modality="siRNA",
            chemistry="PS",
            dose_mg_kg=10.0,  # Very high dose
            schedule="QW"
        )
        
        for endpoint, score in adjusted["risk_scores"].items():
            assert 0.0 <= score <= 1.0, f"{endpoint} score out of bounds: {score}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
