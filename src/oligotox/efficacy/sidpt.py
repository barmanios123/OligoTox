"""siDPT: Deep learning siRNA scoring"""
import torch
import torch.nn as nn
import numpy as np

class siDPT(nn.Module):
    """Deep neural network for siRNA design and potency prediction."""
    
    def __init__(self, input_dim=100):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )
    
    def forward(self, x):
        return torch.sigmoid(self.net(x)) * 100  # 0-100 score
    
    def predict(self, sequences):
        """Predict siRNA potency scores."""
        self.eval()
        # Simplified feature extraction
        features = [self._featurize(seq) for seq in sequences]
        x = torch.tensor(np.array(features), dtype=torch.float32)
        with torch.no_grad():
            scores = self(x)
        return scores.cpu().numpy()
    
    def _featurize(self, seq):
        """Convert sequence to feature vector."""
        # Simplified one-hot + positional encoding
        feat = np.zeros(100)
        for i, nuc in enumerate(seq[:25]):
            idx = {'A':0, 'C':1, 'G':2, 'T':3, 'U':3}.get(nuc, 4)
            if i*4+idx < 100:
                feat[i*4+idx] = 1
        return feat
