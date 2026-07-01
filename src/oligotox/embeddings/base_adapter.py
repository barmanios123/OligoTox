"""Base embedding adapter interface."""
from abc import ABC, abstractmethod
from typing import List, Optional
import numpy as np

class BaseEmbeddingAdapter(ABC):
    """Abstract base class for sequence embedding adapters."""
    
    def __init__(self, model_name: str, embedding_dim: int = 512):
        self.model_name = model_name
        self.embedding_dim = embedding_dim
    
    @abstractmethod
    def encode(self, sequence: str) -> np.ndarray:
        """Encode a single oligonucleotide sequence.
        
        Args:
            sequence: Oligo sequence string
            
        Returns:
            Embedding vector as numpy array
        """
        pass
    
    @abstractmethod
    def encode_batch(self, sequences: List[str]) -> np.ndarray:
        """Encode multiple sequences in batch.
        
        Args:
            sequences: List of oligo sequences
            
        Returns:
            Matrix of embeddings (n_sequences x embedding_dim)
        """
        pass
    
    def get_model_info(self):
        """Get model information."""
        return {
            'name': self.model_name,
            'embedding_dim': self.embedding_dim
        }
