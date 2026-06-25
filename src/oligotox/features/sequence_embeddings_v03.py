"""
Integration 1: Nucleotide/RNA Foundation Model Embeddings
==========================================================
Supports: Nucleotide Transformer (500M / 2.5B), RNAErnie, RiNALMo, BiRNA-BERT

Ref: Dalla-Torre et al., Nat Methods 2024 (Nucleotide Transformer)
     Shuai et al., 2023 (RNAErnie)
     Penic et al., 2024 (RiNALMo)

Install: pip install transformers torch
"""
from __future__ import annotations
import numpy as np
import logging
from typing import List, Optional, Dict

logger = logging.getLogger(__name__)

HF_REPOS = {
    "nucleotide_transformer": "InstaDeepAI/nucleotide-transformer-500m-human-ref",
    "nucleotide_transformer_2500m": "InstaDeepAI/nucleotide-transformer-2.5b-multi-species",
    "rna_ernie":  "multimolecule/rnaernie",
    "rina_lmo":   "lbcb-sci/RiNALMo",
    "birna_bert": "multimolecule/rnamsm",
}


class SequenceEmbedder:
    """
    Generates per-sequence embeddings using HuggingFace RNA/DNA foundation models.
    Falls back to k-mer frequency features if transformers not available.
    """

    def __init__(
        self,
        model: str = "nucleotide_transformer",
        device: str = "cpu",
        n_pca: int = 32,
        batch_size: int = 8,
    ):
        self.model_name = model
        self.device = device
        self.n_pca = n_pca
        self.batch_size = batch_size
        self._model = None
        self._tokenizer = None
        self._pca = None
        self._loaded = False

    def _load(self):
        if self._loaded:
            return
        try:
            import torch
            from transformers import AutoTokenizer, AutoModel
            repo = HF_REPOS.get(self.model_name, self.model_name)
            logger.info(f"Loading {self.model_name} from {repo}...")
            self._tokenizer = AutoTokenizer.from_pretrained(repo, trust_remote_code=True)
            self._model = AutoModel.from_pretrained(repo, trust_remote_code=True).to(self.device)
            self._model.eval()
            self._loaded = True
            logger.info(f"Model loaded on {self.device}")
        except ImportError:
            logger.warning("transformers/torch not installed. Using k-mer fallback.")
            self._loaded = "fallback"
        except Exception as e:
            logger.warning(f"Model load failed ({e}). Using k-mer fallback.")
            self._loaded = "fallback"

    def _kmer_fallback(self, sequence: str, k: int = 4) -> np.ndarray:
        """4-mer frequency vector — fallback when transformers unavailable."""
        seq = sequence.upper().replace("U", "T")
        kmers = {}
        for i in range(len(seq) - k + 1):
            km = seq[i:i+k]
            kmers[km] = kmers.get(km, 0) + 1
        total = max(sum(kmers.values()), 1)
        # Use alphabetically sorted 4-mers as fixed vocabulary
        import itertools
        vocab = ["".join(x) for x in itertools.product("ACGT", repeat=k)]
        vec = np.array([kmers.get(km, 0) / total for km in vocab], dtype=np.float32)
        return vec

    def embed_sequences(self, sequences: List[str]) -> Dict[str, np.ndarray]:
        """
        Returns dict of feature arrays (one per embedding dimension).
        Keys: emb_0, emb_1, ..., emb_{n_pca-1}
        """
        self._load()

        if self._loaded == "fallback":
            raw = np.array([self._kmer_fallback(s) for s in sequences])
            # PCA on k-mer vectors if enough samples
            if raw.shape[0] > self.n_pca:
                from sklearn.decomposition import PCA
                pca = PCA(n_components=self.n_pca)
                raw = pca.fit_transform(raw)
            else:
                raw = raw[:, :self.n_pca]
            return {f"emb_{i}": raw[:, i] for i in range(raw.shape[1])}

        import torch
        embeddings = []
        for i in range(0, len(sequences), self.batch_size):
            batch = sequences[i:i + self.batch_size]
            tokens = self._tokenizer(
                batch, return_tensors="pt", padding=True, truncation=True, max_length=512
            ).to(self.device)
            with torch.no_grad():
                out = self._model(**tokens)
                # Mean pool over sequence dimension
                emb = out.last_hidden_state.mean(dim=1).cpu().numpy()
            embeddings.append(emb)
        raw = np.vstack(embeddings)

        # PCA to n_pca dimensions
        from sklearn.decomposition import PCA
        if self._pca is None and raw.shape[0] > 1:
            self._pca = PCA(n_components=min(self.n_pca, raw.shape[0] - 1, raw.shape[1]))
            raw = self._pca.fit_transform(raw)
        elif self._pca is not None:
            raw = self._pca.transform(raw)
        else:
            raw = raw[:, :self.n_pca]

        return {f"emb_{i}": raw[:, i] for i in range(raw.shape[1])}


if __name__ == "__main__":
    seqs = ["ACGUACGUACGUACGUACGU", "GCGCGCGCGCGCGCGCGCGC"]
    emb = SequenceEmbedder(model="nucleotide_transformer", device="cpu", n_pca=8)
    out = emb.embed_sequences(seqs)
    print(f"Embedding keys: {list(out.keys())}")
    print(f"First embedding: {out['emb_0']}")
