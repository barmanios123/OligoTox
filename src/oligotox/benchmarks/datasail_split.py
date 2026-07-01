"""
DataSAIL-inspired Data Splitting
Leakage-resistant splitting strategies for oligonucleotide data
"""
import pandas as pd
import numpy as np
from typing import Tuple, List, Optional
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split

class DataSAILSplitter:
    """Leakage-resistant data splitting for oligonucleotide therapeutics.
    
    Implements strategies inspired by DataSAIL to prevent data leakage
    through sequence similarity or temporal relationships.
    """
    
    def __init__(self, random_state: int = 42):
        """Initialize DataSAIL splitter.
        
        Args:
            random_state: Random seed for reproducibility
        """
        self.random_state = random_state
        np.random.seed(random_state)
    
    def split(self,
             data: pd.DataFrame,
             strategy: str = "cluster",
             test_size: float = 0.2,
             val_size: float = 0.1,
             sequence_col: str = "sequence") -> Tuple[List[int], List[int], List[int]]:
        """Split data using specified strategy.
        
        Args:
            data: DataFrame with sequences
            strategy: Splitting strategy ('random', 'cluster', 'temporal')
            test_size: Fraction for test set
            val_size: Fraction for validation set
            sequence_col: Name of sequence column
        
        Returns:
            Tuple of (train_indices, val_indices, test_indices)
        """
        n_samples = len(data)
        indices = np.arange(n_samples)
        
        if strategy == "random":
            return self._random_split(indices, test_size, val_size)
        elif strategy == "cluster":
            return self._cluster_split(data, indices, test_size, val_size, sequence_col)
        elif strategy == "temporal":
            return self._temporal_split(data, indices, test_size, val_size)
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
    
    def _random_split(self,
                     indices: np.ndarray,
                     test_size: float,
                     val_size: float) -> Tuple[List[int], List[int], List[int]]:
        """Standard random split."""
        train_val_idx, test_idx = train_test_split(
            indices,
            test_size=test_size,
            random_state=self.random_state
        )
        
        if val_size > 0:
            val_size_adjusted = val_size / (1 - test_size)
            train_idx, val_idx = train_test_split(
                train_val_idx,
                test_size=val_size_adjusted,
                random_state=self.random_state
            )
        else:
            train_idx = train_val_idx
            val_idx = np.array([], dtype=int)
        
        return train_idx.tolist(), val_idx.tolist(), test_idx.tolist()
    
    def _cluster_split(self,
                      data: pd.DataFrame,
                      indices: np.ndarray,
                      test_size: float,
                      val_size: float,
                      sequence_col: str) -> Tuple[List[int], List[int], List[int]]:
        """Cluster-based split to prevent sequence similarity leakage."""
        # Create sequence features (simple k-mer encoding)
        features = self._encode_sequences(data[sequence_col])
        
        # Cluster sequences
        n_clusters = max(10, int(len(data) * 0.1))  # ~10% clusters
        kmeans = KMeans(n_clusters=n_clusters, random_state=self.random_state)
        clusters = kmeans.fit_predict(features)
        
        # Split by clusters
        unique_clusters = np.unique(clusters)
        n_test_clusters = int(len(unique_clusters) * test_size)
        n_val_clusters = int(len(unique_clusters) * val_size)
        
        np.random.shuffle(unique_clusters)
        
        test_clusters = unique_clusters[:n_test_clusters]
        val_clusters = unique_clusters[n_test_clusters:n_test_clusters + n_val_clusters]
        train_clusters = unique_clusters[n_test_clusters + n_val_clusters:]
        
        train_idx = indices[np.isin(clusters, train_clusters)]
        val_idx = indices[np.isin(clusters, val_clusters)]
        test_idx = indices[np.isin(clusters, test_clusters)]
        
        return train_idx.tolist(), val_idx.tolist(), test_idx.tolist()
    
    def _temporal_split(self,
                       data: pd.DataFrame,
                       indices: np.ndarray,
                       test_size: float,
                       val_size: float) -> Tuple[List[int], List[int], List[int]]:
        """Temporal split (assumes data is time-ordered)."""
        n = len(indices)
        n_test = int(n * test_size)
        n_val = int(n * val_size)
        
        test_idx = indices[-n_test:]
        val_idx = indices[-(n_test + n_val):-n_test] if n_val > 0 else np.array([], dtype=int)
        train_idx = indices[:-(n_test + n_val)] if n_val > 0 else indices[:-n_test]
        
        return train_idx.tolist(), val_idx.tolist(), test_idx.tolist()
    
    def _encode_sequences(self, sequences: pd.Series) -> np.ndarray:
        """Simple k-mer encoding for clustering."""
        # Create 3-mer features
        kmers = ['AAA', 'AAU', 'AAG', 'AAC', 'AUA', 'AUU', 'AUG', 'AUC',
                 'AGA', 'AGU', 'AGG', 'AGC', 'ACA', 'ACU', 'ACG', 'ACC',
                 'UAA', 'UAU', 'UAG', 'UAC', 'UUA', 'UUU', 'UUG', 'UUC',
                 'UGA', 'UGU', 'UGG', 'UGC', 'UCA', 'UCU', 'UCG', 'UCC',
                 'GAA', 'GAU', 'GAG', 'GAC', 'GUA', 'GUU', 'GUG', 'GUC',
                 'GGA', 'GGU', 'GGG', 'GGC', 'GCA', 'GCU', 'GCG', 'GCC',
                 'CAA', 'CAU', 'CAG', 'CAC', 'CUA', 'CUU', 'CUG', 'CUC',
                 'CGA', 'CGU', 'CGG', 'CGC', 'CCA', 'CCU', 'CCG', 'CCC']
        
        features = np.zeros((len(sequences), len(kmers)))
        
        for i, seq in enumerate(sequences):
            for j, kmer in enumerate(kmers):
                features[i, j] = seq.count(kmer) / max(1, len(seq) - 2)
        
        return features

