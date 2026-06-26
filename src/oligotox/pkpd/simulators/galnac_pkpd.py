"""GalNAc-conjugated oligonucleotide PKPD simulator"""
import numpy as np
from scipy.integrate import odeint
from typing import Dict, List

class GalNacPKPD:
    """Physiologically-based PK/PD model for GalNAc-oligonucleotides."""
    
    def __init__(self):
        # Default PK parameters for GalNAc-siRNA
        self.params = {
            'ka': 0.5,  # absorption rate
            'CL_plasma': 5.0,  # plasma clearance (L/h)
            'V_plasma': 3.0,  # plasma volume (L)
            'k_liver_uptake': 2.0,  # liver uptake rate
            'k_endosome': 1.5,  # endosomal release
            'k_risc': 0.8,  # RISC loading
            'k_deg': 0.05,  # degradation
            'IC50': 10.0,  # nM
        }
    
    def simulate(self, 
                 dose_mg: float,
                 duration_hours: float = 168,
                 dosing_schedule: List[float] = None) -> Dict:
        """Simulate PK/PD over time.
        
        Args:
            dose_mg: Dose in mg
            duration_hours: Simulation duration
            dosing_schedule: List of dosing times (hours)
            
        Returns:
            Dict with time, plasma_conc, liver_conc, target_knockdown
        """
        if dosing_schedule is None:
            dosing_schedule = [0]  # Single dose
        
        # Time points
        t = np.linspace(0, duration_hours, 1000)
        
        # Initial conditions [plasma, liver, endosome, RISC, target]
        y0 = [0, 0, 0, 0, 100]  # 100% baseline target
        
        # Solve ODE
        solution = odeint(self._pk_pd_model, y0, t, args=(self.params, dose_mg, dosing_schedule))
        
        return {
            'time_hours': t,
            'plasma_conc_nM': solution[:, 0],
            'liver_conc_nM': solution[:, 1],
            'risc_loaded_nM': solution[:, 3],
            'target_knockdown_pct': 100 - solution[:, 4],
            'auc_plasma': np.trapz(solution[:, 0], t),
            'max_knockdown': 100 - solution[:, 4].min()
        }
    
    def _pk_pd_model(self, y, t, params, dose, schedule):
        """ODE system for GalNAc PK/PD."""
        plasma, liver, endo, risc, target = y
        
        # Dosing function
        dose_rate = 0
        for dose_time in schedule:
            if abs(t - dose_time) < 0.1:
                dose_rate = dose * 1000 / params['V_plasma']  # mg to nM
        
        # ODEs
        dplasma = dose_rate - params['CL_plasma']*plasma/params['V_plasma'] - params['k_liver_uptake']*plasma
        dliver = params['k_liver_uptake']*plasma - params['k_endosome']*liver
        dendo = params['k_endosome']*liver - params['k_risc']*endo
        drisc = params['k_risc']*endo - params['k_deg']*risc
        
        # Target knockdown (mRNA)
        knockdown_rate = risc / (params['IC50'] + risc)
        dtarget = -knockdown_rate * target + 0.05 * (100 - target)  # synthesis
        
        return [dplasma, dliver, dendo, drisc, dtarget]
