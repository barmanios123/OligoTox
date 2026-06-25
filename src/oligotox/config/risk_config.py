# config/risk_config.py

ENDPOINT_SEVERITY_WEIGHTS = {
    "hepatocellular_injury": 1.00,
    "hepatic_vacuolation": 0.55,
    "renal_tubular_injury": 0.95,
    "renal_accumulation": 0.50,
    "thrombocytopenia": 0.95,
    "coagulation_effect": 0.90,
    "complement_activation": 0.85,
    "innate_immune_activation": 0.70,
    "cytokine_release": 0.85,
    "kupffer_cell_activation": 0.65,
    "cns_toxicity": 1.00,
    "offtarget_gene_silencing": 0.75,
    "injection_site_reaction": 0.30,
}

DEFAULT_THRESHOLD_POLICY = {
    "method": "f1",
    "min_threshold": 0.20,
    "max_threshold": 0.80,
    "fallback_threshold": 0.50,
}

RISK_BANDS = {
    "low": 0.25,
    "medium": 0.50,
    "high": 0.75,
}