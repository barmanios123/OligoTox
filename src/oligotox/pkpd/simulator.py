"""
GalNAc-siRNA PK/PD Simulator
Phase 2A: 2-compartment model for liver-targeted oligonucleotides
"""
import numpy as np
from typing import Dict, Literal


class GalNAcPKPDSimulator:
    """
    2-compartment PK model for GalNAc-siRNA conjugates.
    
    Models:
    - Plasma compartment with first-order clearance
    - Liver compartment with GalNAc-mediated uptake
    - Tissue distribution and terminal elimination
    """
    
    def __init__(self, species: str = "human", body_weight_kg: float = 70.0):
        self.species = species
        self.body_weight_kg = body_weight_kg
        
        # Human PK parameters (typical GalNAc-siRNA values)
        self.CL_liver = 0.15  # L/h/kg - liver clearance
        self.Vd_liver = 0.5   # L/kg - liver volume of distribution
        self.ka_liver = 2.5   # 1/h - absorption rate to liver
        self.f_liver = 0.85   # fraction reaching liver via GalNAc
        
    def simulate(
        self,
        dose_mg_kg: float,
        schedule: Literal["QW", "Q2W", "Q4W", "Q8W", "Q12W"]
    ) -> Dict[str, float]:
        """
        Simulate PK/PD for given dose and schedule.
        
        Args:
            dose_mg_kg: Dose in mg/kg
            schedule: Dosing frequency (QW=weekly, Q4W=monthly, etc.)
        
        Returns:
            Dict with PK/PD metrics: AUC, Cmax, t1/2, knockdown
        """
        # Convert schedule to dosing interval (hours)
        interval_map = {
            "QW": 168,    # 7 days
            "Q2W": 336,   # 14 days
            "Q4W": 672,   # 28 days
            "Q8W": 1344,  # 56 days
            "Q12W": 2016  # 84 days
        }
        dosing_interval_h = interval_map[schedule]
        
        # Calculate liver exposure
        dose_total_mg = dose_mg_kg * self.body_weight_kg
        dose_liver_mg = dose_total_mg * self.f_liver
        
        # Convert to pmol (assuming MW ~13,000 Da for siRNA duplex)
        MW_DALTONS = 13000
        dose_liver_pmol = (dose_liver_mg / MW_DALTONS) * 1e9
        
        # Liver concentration (pmol/g, assuming 1.8 kg liver ~ 1800 g)
        liver_mass_g = 1800
        C0_liver = dose_liver_pmol / liver_mass_g
        
        # Calculate PK metrics
        k_elim = self.CL_liver / self.Vd_liver  # elimination rate constant
        t_half = 0.693 / k_elim
        
        # AUC over dosing interval (trapezoid)
        if dosing_interval_h < 10 * t_half:  # steady-state considerations
            auc_liver = C0_liver * (1 - np.exp(-k_elim * dosing_interval_h)) / k_elim
        else:
            auc_liver = C0_liver / k_elim  # complete elimination
        
        # Cmax (peak concentration)
        cmax_liver = C0_liver
        
        # PD: estimate target knockdown (simplified Emax model)
        EC50 = 50  # pmol/g (typical IC50 for potent siRNA)
        knockdown_pct = 100 * (cmax_liver / (cmax_liver + EC50))
        
        return {
            "liver_auc_pmol_g_h": auc_liver,
            "liver_cmax_pmol_g": cmax_liver,
            "liver_half_life_h": t_half,
            "peak_knockdown_pct": min(knockdown_pct, 95)  # cap at 95%
        }
