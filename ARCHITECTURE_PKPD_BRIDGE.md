# GalNAc PK/PD Bridge: Technical Architecture and Implementation Roadmap

## Overview

This document provides the technical architecture and phased implementation plan for integrating GalNAc PK/PD simulation capabilities into the OligoTox platform as a **modular, optional enhancement** for liver-targeted, ASGPR-mediated siRNA therapeutics.

**Design Principles:**
1. **Modular**: GalNAc PK/PD bridge operates as an independent module, not a core rewrite
2. **Optional**: Users can run OligoTox without PK/PD for non-GalNAc candidates
3. **Backward-compatible**: Existing toxicity prediction workflow remains unchanged
4. **Validated incrementally**: Each phase includes validation against holdout data or prospective predictions

---

## System Architecture

### High-Level Component Diagram

```
┌────────────────────────────────────────────────────────────────────────────┐
│                         OligoTox Platform (v2.0)                           │
├────────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  ┌─────────────────────────────────────────────────────────────────────┐  │
│  │                   Core Toxicity Engine (v1.x)                       │  │
│  │  - src/oligotox/inference/api.py: run_oligotox()                    │  │
│  │  - Multi-endpoint models (12+ endpoints)                            │  │
│  │  - D1/D2/D3 mechanistic attribution                                 │  │
│  │  - Risk bands, thresholds, decision rules                           │  │
│  └─────────────────────────────────────────────────────────────────────┘  │
│                              ▲                                             │
│                              │                                             │
│                              │ PK/PD-enriched features                    │
│                              │ (liver_auc, cmax, t_half, exposure_score)  │
│                              │                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐  │
│  │          GalNAc PK/PD Bridge Module (NEW - v2.0)                    │  │
│  │  src/oligotox/pkpd/                                                 │  │
│  │    ├── simulator.py        # Core PK/PD simulation engine           │  │
│  │    ├── parameters.py       # Species-specific PK parameters         │  │
│  │    ├── models.py           # Compartmental model implementations    │  │
│  │    ├── sensitivity.py      # Parameter sensitivity analysis         │  │
│  │    └── bridge.py           # Integration with tox engine            │  │
│  └─────────────────────────────────────────────────────────────────────┘  │
│                              │                                             │
│                              │ Simulated PK/PD outputs                    │
│                              ▼                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐  │
│  │             Enhanced Risk Interpretation Layer (NEW)                │  │
│  │  - Dose-aware safety scoring                                        │  │
│  │  - Regimen comparison (e.g., weekly vs. daily)                      │  │
│  │  - Cross-species safety projection                                  │  │
│  │  - Joint efficacy-safety visualization                              │  │
│  └─────────────────────────────────────────────────────────────────────┘  │
│                                                                            │
└────────────────────────────────────────────────────────────────────────────┘

                                    ▲
                                    │
                                    │ User inputs via Streamlit UI
                                    │

┌────────────────────────────────────────────────────────────────────────────┐
│                    Streamlit UI (oligotox_app.py)                          │
├────────────────────────────────────────────────────────────────────────────┤
│  Existing Tabs:                          NEW Tab:                          │
│  - Sequence Analysis                     - Translational PK/PD (GalNAc)    │
│  - Overview                                                                │
│  - Curation                              Inputs:                           │
│  - Runs                                  - Sequence + modality/chem        │
│  - Models                                - Dose (mg/kg)                    │
│  - Predictions                           - Schedule (QW, Q2W, Q4W)         │
│  - Evaluation                            - Species (mouse, rat, NHP, human)│
│                                          - PK parameters (or use defaults) │
│                                                                            │
│                                          Outputs:                          │
│                                          - Plasma PK curve                 │
│                                          - Liver PK curve                  │
│                                          - mRNA knockdown (%)              │
│                                          - Protein suppression (%)         │
│                                          - Dose-adjusted tox risk profile  │
│                                          - Sensitivity analysis plots      │
└────────────────────────────────────────────────────────────────────────────┘
```

---

## Module Specifications

### 1. PK/PD Simulator Core (`src/oligotox/pkpd/simulator.py`)

**Purpose**: Simulate GalNAc-siRNA pharmacokinetics and pharmacodynamics across compartments and species.

**Key Classes**:

```python
class GalNAcPKPDSimulator:
    """Main simulation engine for GalNAc-siRNA PK/PD."""
    
    def __init__(self, species: str = "human", parameters: Optional[Dict] = None):
        """Initialize with species-specific default parameters."""
        pass
    
    def simulate_plasma_pk(self, dose_mg_kg: float, schedule: str) -> pd.DataFrame:
        """Simulate plasma concentration over time.
        
        Returns:
            DataFrame with columns: time_h, conc_nM, auc_h_nM, cmax_nM
        """
        pass
    
    def simulate_liver_pk(self, dose_mg_kg: float, schedule: str) -> pd.DataFrame:
        """Simulate liver concentration over time (total siRNA + RISC-loaded).
        
        Returns:
            DataFrame with columns: time_h, total_conc_pmol_g, risc_conc_pmol_g, 
                                     auc_h_pmol_g, cmax_pmol_g
        """
        pass
    
    def simulate_target_knockdown(self, dose_mg_kg: float, schedule: str) -> pd.DataFrame:
        """Simulate mRNA knockdown and protein suppression.
        
        Returns:
            DataFrame with columns: time_h, mrna_remaining_pct, protein_remaining_pct
        """
        pass
    
    def get_exposure_summary(self, dose_mg_kg: float, schedule: str) -> Dict:
        """Compute summary PK/PD metrics for toxicity bridge.
        
        Returns:
            {
                'plasma_auc_h_nM': float,
                'liver_auc_h_pmol_g': float,
                'liver_cmax_pmol_g': float,
                'liver_t_half_h': float,
                'peak_knockdown_pct': float,
                'duration_above_80pct_knockdown_days': float,
                'exposure_score': float  # Composite risk-relevant metric
            }
        """
        pass
```

**Compartmental Model** (simplified from MTNA paper):

1. **Plasma compartment**: First-order absorption → distribution to liver/tissues
2. **Liver compartment**: ASGPR-mediated uptake → endosomal processing → cytoplasmic release
3. **RISC compartment**: siRNA loading into RISC → mRNA cleavage → degradation
4. **Target mRNA/protein**: Knockdown kinetics with turnover rates

**Default Parameters** (tuned to inclisiran/givosiran data):
- `k_abs`: absorption rate constant (1/h)
- `CL_plasma`: plasma clearance (L/h/kg)
- `V_plasma`: plasma volume of distribution (L/kg)
- `k_uptake_liver`: ASGPR uptake rate (1/h)
- `f_escape`: endosomal escape fraction (0.01-0.05)
- `k_deg_liver`: liver siRNA degradation rate (1/h)
- `k_risc_loading`: RISC loading rate (1/h)
- `k_mrna_turnover`: Target mRNA turnover rate (1/h)
- `k_protein_turnover`: Target protein turnover rate (1/h)

---

### 2. Species-Specific Parameters (`src/oligotox/pkpd/parameters.py`)

**Purpose**: Provide validated PK/PD parameters for cross-species translation.

```python
SPECIES_PARAMETERS = {
    "mouse": {
        "body_weight_kg": 0.025,
        "liver_weight_g": 1.2,
        "asgpr_expression_rel": 1.0,
        "CL_plasma_L_h_kg": 0.5,
        "V_plasma_L_kg": 0.05,
        # ... other parameters
    },
    "rat": { ... },
    "cynomolgus": { ... },
    "human": {
        "body_weight_kg": 70,
        "liver_weight_g": 1800,
        "asgpr_expression_rel": 1.0,
        "CL_plasma_L_h_kg": 0.2,
        "V_plasma_L_kg": 0.07,
        # ... other parameters from MTNA validation
    },
}

def get_species_parameters(species: str, custom_overrides: Optional[Dict] = None) -> Dict:
    """Retrieve species parameters with optional custom overrides."""
    pass

def scale_dose_across_species(
    dose_mg_kg: float, 
    from_species: str, 
    to_species: str, 
    scaling_method: str = "auc"
) -> float:
    """Scale dose from one species to another (AUC-based or body surface area)."""
    pass
```

---

### 3. Toxicity Bridge (`src/oligotox/pkpd/bridge.py`)

**Purpose**: Connect PK/PD simulation outputs to toxicity risk adjustments.

```python
def compute_pkpd_enriched_features(
    sequence: str,
    modality: str,
    chemistry: str,
    dose_mg_kg: float,
    schedule: str,
    species: str
) -> Dict:
    """Compute PK/PD-derived features for toxicity model input.
    
    Returns:
        {
            'liver_auc_h_pmol_g': float,
            'liver_cmax_pmol_g': float,
            'liver_t_half_h': float,
            'exposure_score': float,
            'peak_knockdown_pct': float,
            'sustained_exposure_days': float
        }
    """
    simulator = GalNAcPKPDSimulator(species=species)
    exposure_summary = simulator.get_exposure_summary(dose_mg_kg, schedule)
    
    # Compute composite exposure score (log-scaled, normalized)
    exposure_score = compute_exposure_score(
        auc=exposure_summary['liver_auc_h_pmol_g'],
        cmax=exposure_summary['liver_cmax_pmol_g'],
        duration=exposure_summary['duration_above_80pct_knockdown_days']
    )
    
    return {
        'liver_auc': exposure_summary['liver_auc_h_pmol_g'],
        'liver_cmax': exposure_summary['liver_cmax_pmol_g'],
        'liver_t_half': exposure_summary['liver_t_half_h'],
        'exposure_score': exposure_score,
        'peak_knockdown_pct': exposure_summary['peak_knockdown_pct'],
        'sustained_exposure_days': exposure_summary['duration_above_80pct_knockdown_days']
    }

def adjust_toxicity_risk_with_pkpd(
    base_risk_profile: Dict,
    pkpd_features: Dict,
    adjustment_model: str = "multiplicative"
) -> Dict:
    """Adjust toxicity risk scores using PK/PD-derived exposure metrics.
    
    Args:
        base_risk_profile: Original OligoTox risk scores (per endpoint)
        pkpd_features: PK/PD-derived features from compute_pkpd_enriched_features()
        adjustment_model: 'multiplicative', 'additive', or 'learned'
    
    Returns:
        Adjusted risk profile with PK/PD context
    """
    adjusted_profile = base_risk_profile.copy()
    
    # Example: hepatotoxicity increases with sustained liver exposure
    if 'hepatic_vacuolation' in adjusted_profile:
        # Multiplicative adjustment based on exposure score
        exposure_factor = 1.0 + (pkpd_features['exposure_score'] - 1.0) * 0.3
        adjusted_profile['hepatic_vacuolation']['prob'] *= exposure_factor
        adjusted_profile['hepatic_vacuolation']['risk_band'] = compute_risk_band(
            adjusted_profile['hepatic_vacuolation']['prob']
        )
    
    # Similar adjustments for other endpoints (kupffer_cell, renal, immune)
    
    return adjusted_profile
```

---

### 4. Sensitivity Analysis (`src/oligotox/pkpd/sensitivity.py`)

**Purpose**: Identify which PK/PD parameters most influence toxicity risk.

```python
def run_sensitivity_analysis(
    sequence: str,
    dose_mg_kg: float,
    schedule: str,
    species: str,
    parameters_to_vary: List[str],
    variation_range: Tuple[float, float] = (0.5, 2.0),
    n_steps: int = 20
) -> pd.DataFrame:
    """Perform one-at-a-time sensitivity analysis.
    
    Args:
        parameters_to_vary: e.g., ['f_escape', 'k_deg_liver', 'k_uptake_liver']
        variation_range: (min_factor, max_factor) relative to default
    
    Returns:
        DataFrame with columns: parameter, value, endpoint, risk_score, knockdown_pct
    """
    pass

def plot_efficacy_safety_frontier(
    sequence: str,
    dose_range: Tuple[float, float],
    schedule: str,
    species: str,
    target_endpoints: List[str] = ['hepatic_vacuolation', 'kupffer_cell_activation']
) -> plt.Figure:
    """Plot joint efficacy (knockdown %) vs. safety (risk scores) across dose range."""
    pass
```

---

## Implementation Roadmap

### Phase 2A: Minimal PK/PD Bridge (Weeks 1-4)

**Goal**: Prove value of PK/PD enrichment with minimal complexity.

#### Week 1: Setup and Core Simulator
- [ ] Create `src/oligotox/pkpd/` module structure
- [ ] Implement simplified 2-compartment PK model (plasma + liver)
- [ ] Define default human parameters based on MTNA paper (inclisiran, givosiran)
- [ ] Unit tests for plasma/liver PK simulation

#### Week 2: Exposure Metrics and Bridge
- [ ] Implement `get_exposure_summary()` to compute AUC, Cmax, t_half
- [ ] Implement `compute_pkpd_enriched_features()` for toxicity bridge
- [ ] Simple multiplicative adjustment for hepatotoxicity endpoints
- [ ] Unit tests for bridge functions

#### Week 3: UI Integration
- [ ] Add "Advanced: PK/PD-Aware Scoring (GalNAc only)" checkbox to Sequence Analysis tab
- [ ] If enabled, collect dose (mg/kg) and schedule (dropdown: QW, Q2W, Q4W)
- [ ] Display PK/PD-adjusted risk profile alongside base profile
- [ ] Show exposure summary table (AUC, Cmax, t_half)

#### Week 4: Validation
- [ ] Retrospective validation: Do PK/PD adjustments improve prediction for known hepatotoxic vs. safe GalNAc-siRNAs?
- [ ] Compute AUROC improvement for hepatotoxicity endpoints
- [ ] Document results in validation report
- [ ] **Go/No-Go Decision**: Proceed to Phase 2B only if validation shows ≥5% AUROC improvement

---

### Phase 2B: Cross-Species Scaling (Weeks 5-10)

**Goal**: Enable dose projection and safety assessment across preclinical species.

#### Week 5-6: Species Parameters
- [ ] Implement `parameters.py` with mouse, rat, cynomolgus, human parameters
- [ ] Validate species scaling against published GalNAc-siRNA data
- [ ] Implement `scale_dose_across_species()` with AUC-matching

#### Week 7-8: Cross-Species Risk Projection
- [ ] Extend bridge to accept species parameter
- [ ] Add species dropdown to PK/PD UI tab
- [ ] Display dose-escalation table: "Mouse safe dose → Rat equivalent → NHP equivalent → Human projected"
- [ ] For each species, show adjusted risk profile

#### Week 9-10: Validation and Documentation
- [ ] Validate cross-species scaling with known IND-enabling studies
- [ ] Create case study: "Species-aware dose selection for [Candidate X]"
- [ ] Update PRODUCTION_READINESS.md with cross-species validation results

---

### Phase 2C: Full Sensitivity Analysis (Weeks 11-18)

**Goal**: Enable mechanistic parameter optimization for safety-guided design.

#### Week 11-13: Sensitivity Analysis Engine
- [ ] Implement `run_sensitivity_analysis()` with parameter sweeps
- [ ] Implement `plot_efficacy_safety_frontier()` for dose-response visualization
- [ ] Add Monte Carlo sampling option for parameter uncertainty

#### Week 14-16: New UI Tab - "Translational PK/PD"
- [ ] Create dedicated Streamlit tab (separate from Sequence Analysis)
- [ ] Inputs: sequence, dose range, schedule, species, custom parameters (advanced)
- [ ] Outputs:
  - Interactive PK/PD curve plots (plasma, liver, knockdown)
  - Sensitivity heatmaps (parameter vs. endpoint risk)
  - Efficacy-safety frontier plot
  - Recommended dose range with safety justification

#### Week 17-18: Partner Validation and Polish
- [ ] Beta test with internal GalNAc programs or external collaborators
- [ ] Incorporate feedback on parameter defaults and UI usability
- [ ] Finalize documentation: "GalNAc PK/PD Bridge User Guide"
- [ ] Prepare slide deck for partner presentations

---

## Data Requirements

### For Validation (Phase 2A)
- **GalNAc-siRNA training set** (if available):
  - Sequences with known hepatotoxicity outcomes
  - Dose/schedule information from IND or preclinical studies
  - Species-specific safety signals
- **External benchmark**: Published GalNAc-siRNA PK/PD data (inclisiran, givosiran, fitusiran, etc.)

### For Cross-Species Scaling (Phase 2B)
- Species-specific PK parameters from MTNA paper or internal databases
- Validation cohort: At least 3-5 GalNAc-siRNAs with cross-species tox data

---

## Testing Strategy

### Unit Tests (`tests/pkpd/`)
```python
def test_plasma_pk_simulation():
    """Test that plasma PK matches expected kinetics."""
    simulator = GalNAcPKPDSimulator(species="human")
    plasma_pk = simulator.simulate_plasma_pk(dose_mg_kg=3.0, schedule="Q4W")
    assert plasma_pk['cmax_nM'].max() > 0
    assert plasma_pk['auc_h_nM'].sum() > 0

def test_exposure_score_calculation():
    """Test exposure score is monotonic with AUC and Cmax."""
    # ...

def test_toxicity_adjustment_increases_risk_with_high_exposure():
    """Test that high exposure increases hepatotoxicity risk."""
    # ...
```

### Integration Tests
```python
def test_end_to_end_pkpd_enriched_prediction():
    """Test full workflow: sequence → PK/PD sim → adjusted tox risk."""
    result = run_oligotox_with_pkpd(
        sequence="UGCAUGCUAGCUAGCUAGCU",
        modality="siRNA",
        chemistry="PS",
        conjugate="GalNAc",
        dose_mg_kg=3.0,
        schedule="Q4W",
        species="human"
    )
    assert 'pkpd_adjusted_risk' in result
    assert result['pkpd_adjusted_risk']['hepatic_vacuolation']['prob'] != result['base_risk']['hepatic_vacuolation']['prob']
```

---

## Risk Mitigation

### Technical Risks
1. **Parameter identifiability**: PK/PD models may be underspecified without rich training data
   - Mitigation: Start with literature-validated defaults; allow user-specified overrides
   
2. **Overfitting to small GalNAc dataset**: Toxicity adjustments may not generalize
   - Mitigation: Phase 2A validation with strict go/no-go criteria; external prospective validation

3. **UI complexity**: Too many parameters may confuse users
   - Mitigation: Default "Simple Mode" with dose/schedule only; "Advanced Mode" for custom parameters

### Operational Risks
1. **Scope creep**: Expanding to non-GalNAc modalities prematurely
   - Mitigation: Strict modality gating (error if user selects LNP-siRNA or ASO with PK/PD enabled)

2. **Validation delays**: Insufficient data to validate cross-species scaling
   - Mitigation: Phase 2B is optional; can launch Phase 2A alone if validation succeeds

---

## Success Metrics

### Phase 2A
- ✅ AUROC improvement ≥5% for hepatotoxicity endpoints with PK/PD enrichment
- ✅ User feedback: "PK/PD-aware scoring adds value" (≥4/5 rating from beta testers)

### Phase 2B
- ✅ Cross-species dose projections within 2-fold of observed safe doses (≥80% cases)
- ✅ At least 1 prospective validation: predicted human dose aligns with IND-selected dose

### Phase 2C
- ✅ Sensitivity analysis identifies ≥2 actionable design parameters per candidate
- ✅ Partner adoption: ≥1 external collaborator integrates GalNAc PK/PD module into internal pipeline

---

## Maintenance and Governance

### Model Updates
- **Quarterly review** of default PK/PD parameters as new GalNAc-siRNA data becomes available
- **Version control**: PK/PD parameter sets tagged with validation dataset version

### Documentation
- **User Guide**: `docs/GALNAC_PKPD_USER_GUIDE.md` with worked examples
- **API Reference**: Auto-generated from docstrings (Sphinx)
- **Validation Reports**: Versioned markdown files in `docs/validation/`

---

## Appendix: Example Use Cases

### Use Case 1: Dose Optimization for Clinical Candidate
**Scenario**: Lead GalNAc-siRNA shows 90% target knockdown at 5 mg/kg Q4W in NHP, but hepatotoxicity signal at 10 mg/kg.

**Workflow**:
1. Run Sequence Analysis → Base toxicity risk: hepatic_vacuolation = MODERATE (prob=0.55)
2. Enable PK/PD mode → Input: 5 mg/kg Q4W, NHP
3. PK/PD adjustment → liver AUC = 15000 h·pmol/g → exposure_score = 1.2 → adjusted risk = MODERATE-HIGH (prob=0.66)
4. Test alternate schedule: 2.5 mg/kg Q2W → exposure_score = 0.9 → adjusted risk = LOW-MODERATE (prob=0.49)
5. **Decision**: Recommend Q2W schedule for better safety margin

### Use Case 2: Cross-Species Dose Projection for IND
**Scenario**: Mouse safe dose = 10 mg/kg QW, need to project human starting dose.

**Workflow**:
1. Input sequence + dose/schedule + species = mouse
2. Simulate mouse PK → liver AUC = 8000 h·pmol/g
3. Scale to human using AUC-matching → projected dose = 0.8 mg/kg QW
4. Run human PK/PD simulation at 0.8 mg/kg → exposure_score = 0.7 → hepatic risk = LOW
5. **Decision**: Propose 0.8 mg/kg QW as starting dose with adequate safety margin

### Use Case 3: Chemistry Optimization
**Scenario**: Compare two GalNAc linker designs (high escape vs. high stability).

**Workflow**:
1. Sequence A: f_escape = 0.05, k_deg = 0.02 → Run sensitivity analysis
2. Sequence B: f_escape = 0.02, k_deg = 0.01 → Run sensitivity analysis
3. Compare efficacy-safety frontiers:
   - Sequence A achieves 90% knockdown at 3 mg/kg but MODERATE hepatotoxicity risk
   - Sequence B requires 5 mg/kg for 90% knockdown but LOW hepatotoxicity risk
4. **Decision**: Select Sequence B for better safety profile despite higher dose requirement

---

**Version**: 1.0  
**Last Updated**: 2025-01-20  
**Status**: Planning (Phase 2A starts post-Phase 1 completion)
