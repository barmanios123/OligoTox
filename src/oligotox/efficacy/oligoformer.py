
"""OligoFormer: Transformer-based oligonucleotide design model"""
import torch
import torch.nn as nn
import numpy as np
from typing import Dict, List, Optional

class OligoFormer(nn.Module):
    """Transformer model for oligo efficacy prediction.

    Based on transformer architecture optimized for oligonucleotide sequences.
    Predicts on-target efficacy and off-target potential.
    """

    def __init__(self, 
                 vocab_size: int = 5,  # A, C, G, T, N
                 d_model: int = 256,
                 nhead: int = 8,
                 num_layers: int = 6,
                 max_len: int = 50):
        super().__init__()

        self.embedding = nn.Embedding(vocab_size, d_model)
        self.pos_encoding = nn.Parameter(torch.randn(1, max_len, d_model))

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=d_model*4,
            dropout=0.1,
            batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers)

        # Prediction heads
        self.efficacy_head = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(128, 1),
            nn.Sigmoid()
        )

        self.offtarget_head = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(128, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        # x: [batch, seq_len] with nucleotide indices
        embedded = self.embedding(x) + self.pos_encoding[:, :x.size(1), :]
        transformed = self.transformer(embedded)

        # Use [CLS]-like token (first position) for classification
        pooled = transformed[:, 0, :]

        efficacy = self.efficacy_head(pooled)
        offtarget = self.offtarget_head(pooled)

        return {
            'efficacy_score': efficacy,
            'offtarget_score': offtarget,
            'embeddings': transformed
        }

    def predict(self, sequences: List[str]) -> Dict[str, np.ndarray]:
        """Predict efficacy for oligonucleotide sequences.

        Args:
            sequences: List of oligo sequences (strings)

        Returns:
            Dict with efficacy_scores and offtarget_scores
        """
        self.eval()
        # Convert sequences to indices
        seq_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

        indices = []
        for seq in sequences:
            idx = [seq_to_idx.get(n.upper(), 4) for n in seq]
            indices.append(idx)

        # Pad sequences
        max_len = max(len(idx) for idx in indices)
        padded = [idx + [4] * (max_len - len(idx)) for idx in indices]

        x = torch.tensor(padded, dtype=torch.long)

        with torch.no_grad():
            out = self(x)

        return {
            'efficacy_scores': out['efficacy_score'].cpu().numpy(),
            'offtarget_scores': out['offtarget_score'].cpu().numpy()
        }
