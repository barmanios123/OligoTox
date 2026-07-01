"""GalNAc-conjugated oligonucleotide PKPD simulator."""
import numpy as np
from typing import Dict, List, Optional

class GalNAcPKPDSimulator:
    """PBPK/PD model for GalNAc-conjugated oligonucleotides."""
    
    def __init__(self):
        # Physiological parameters
        self.body_weight = 70  # kg
        self.liver_volume = 1.5  # L
        self.plasma_volume = 3.0  # L
        
        # PK parameters for GalNAc-ASO
        self.ka = 0.5  # Absorption rate (h^-1)
        self.ke_plasma = 0.05  # Plasma elimination (h^-1)
        self.ke_liver = 0.001  # Liver elimination (h^-1)
        self.liver_uptake_fraction = 0.8
        
    def simulate(self,
                 dose_mg: float,
                 duration_hours: float = 168,
                 dosing_schedule: Optional[List[float]] = None) -> Dict:
        """Simulate PK/PD over time.
        
        Args:
            dose_mg: Dose in mg
            duration_hours: Simulation duration
            dosing_schedule: List of dosing times (hours)
            
        Returns:
            Dictionary with 'time', 'plasma', and 'liver' concentration arrays
        """
        if dosing_schedule is None:
            dosing_schedule = [0]  # Single dose at time 0
        
        # Time grid
        dt = 0.5  # hours
        time = np.arange(0, duration_hours + dt, dt)
        
        # Initialize compartments
        plasma = np.zeros_like(time)
        liver = np.zeros_like(time)
        
        # Convert dose to concentration
        dose_per_dose = dose_mg / self.body_weight  # mg/kg
        
        # Simulate each dose
        for dose_time in dosing_schedule:
            dose_idx = int(dose_time / dt)
            if dose_idx >= len(time):
                continue
                
            # Simple 2-compartment model
            for i in range(dose_idx + 1, len(time)):
                t_since_dose = time[i] - dose_time
                
                # Plasma concentration (absorption and elimination)
                plasma[i] += (dose_per_dose * self.ka * 
                             np.exp(-self.ke_plasma * t_since_dose) * 
                             (1 - np.exp(-self.ka * t_since_dose)))
                
                # Liver accumulation
                liver[i] += (dose_per_dose * self.liver_uptake_fraction * 
                            np.exp(-self.ke_liver * t_since_dose))
        
        return {
            'time': time,
            'plasma': plasma,
            'liver': liver
        }
