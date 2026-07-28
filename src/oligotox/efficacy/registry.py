"""
Efficacy Model Registry - Unified interface for all efficacy prediction models
"""
from typing import Dict, List, Optional
import numpy as np

class EfficacyModelRegistry:
    """Registry for all integrated efficacy prediction models."""
    
    def __init__(self):
        """Initialize registry with all available models."""
        from oligotox.efficacy.oligoformer import OligoFormer
        from oligotox.efficacy.ensirna import ENsiRNA
        from oligotox.efficacy.sidpt import siDPT
        
        self.models = {
            "OligoFormer": OligoFormer(),
            "ENsiRNA": ENsiRNA(),
            "siDPT": siDPT(),
            "PFRED": self._pfred_model(),
            "ASOptimizer": self._asoptimizer_model(),
            "siRNADiscovery": self._sirna_discovery_model(),
        }
    
    def _pfred_model(self):
        """Placeholder for PFRED rules-based model."""
        class PFREDModel:
            def predict(self, sequence, modality="siRNA", chemistry="PS"):
                # Rule-based scoring
                gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
                return 0.7 + 0.2 * gc_content
        return PFREDModel()
    
    def _asoptimizer_model(self):
        """Placeholder for ASO optimizer."""
        class ASOptimizerModel:
            def predict(self, sequence, modality="ASO", chemistry="PS"):
                # Simple ASO scoring
                return 0.65 + 0.15 * np.random.random()
        return ASOptimizerModel()
    
    def _sirna_discovery_model(self):
        """Placeholder for siRNA Discovery."""
        class siRNADiscoveryModel:
            def predict(self, sequence, modality="siRNA", chemistry="PS"):
                return 0.72 + 0.1 * np.random.random()
        return siRNADiscoveryModel()
    
    def list_models(self) -> List[str]:
        """Return list of available model names."""
        return list(self.models.keys())
    
    def get_available_models(self):
        """Get list of available efficacy model names."""
        return list(self.models.keys())

    def get_model(self, name: str):
        """Get a specific model by name."""
        if name not in self.models:
            raise ValueError(f"Unknown model: {name}. Available: {self.list_models()}")
        return self.models[name]
    
    # Models with a sequences-list-only interface: predict(self, sequences)
    _SEQ_LIST_MODELS = {"OligoFormer", "ENsiRNA", "siDPT"}

    def _score_from_seq_list_model(self, name, model, sequence):
        """Call sequence-list style models and reduce output to a single float score."""
        output = model.predict([sequence])
        if isinstance(output, dict):
            scores = output.get("efficacy_scores")
            if scores is None:
                scores = next(iter(output.values()))
            value = float(np.asarray(scores).reshape(-1)[0])
            if value <= 1.0:
                value *= 100.0
            return value
        else:
            value = float(np.asarray(output).reshape(-1)[0])
            return value

    def predict_all(self, sequence: str, modality: str = "siRNA",
                     chemistry: str = "PS") -> Dict[str, float]:
        """Run all compatible models and return results dict."""
        results = {}
        for name, model in self.models.items():
            try:
                if hasattr(model, 'modalities') and modality not in model.modalities:
                    continue
                if name in self._SEQ_LIST_MODELS:
                    results[name] = self._score_from_seq_list_model(name, model, sequence)
                else:
                    results[name] = model.predict(sequence, modality, chemistry)
            except Exception as e:
                print(f"Warning: {name} failed: {e}")
                results[name] = 0.5  # Default fallback
        return results
