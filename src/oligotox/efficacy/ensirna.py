"""ENsiRNA: siRNA efficacy prediction model"""
import numpy as np
from sklearn.ensemble import GradientBoostingRegressor

class ENsiRNA:
    """Ensemble model for siRNA on-target efficacy."""
    
    def __init__(self):
        self.model = GradientBoostingRegressor(
            n_estimators=100,
            learning_rate=0.1,
            max_depth=5,
            random_state=42
        )
        self.is_trained = False
    
    def predict(self, sequences):
        """Predict siRNA efficacy scores (0-100)."""
        # Extract features (simplified for stub)
        features = self._extract_features(sequences)
        
        if not self.is_trained:
            # Return random scores if not trained
            return np.random.uniform(30, 95, len(sequences))
        
        return self.model.predict(features)
    
    def _extract_features(self, sequences):
        """Extract siRNA-specific features."""
        features = []
        for seq in sequences:
            feat = [
                seq.count('G') / len(seq),  # GC content
                seq.count('C') / len(seq),
                seq[:7].count('A'),  # Position-specific
                seq[:7].count('U'),
                len(seq)
            ]
            features.append(feat)
        return np.array(features)
