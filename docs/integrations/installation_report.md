# Integration Environment Setup Report
**Date:** June 30, 2026  
**Python Version:** 3.13  
**Conda Environment:** base

---

## ✅ Successfully Installed Packages

### Core ML/DL Stack
- **PyTorch:** 2.12.1 (with CUDA/MPS support)
- **TorchVision:** 0.27.1
- **TorchAudio:** 2.11.0
- **PyTorch Geometric:** 2.8.0 (for graph neural networks)

### NLP/Transformers
- **Transformers:** 5.12.1 (Hugging Face)
- **Tokenizers:** 0.22.2
- **SentencePiece:** 0.2.1
- **Datasets:** 5.0.0
- **Accelerate:** 1.14.0
- **SafeTensors:** 0.7.0

### Scientific Computing
- **NumPy:** 2.4.6
- **Pandas:** 2.3.3
- **SciPy:** 1.18.0
- **Scikit-learn:** 1.9.0
- **PyArrow:** 21.0.0

### Bioinformatics
- **Biopython:** 1.87

### Web/API Development
- **FastAPI:** 0.138.0
- **Streamlit:** 1.51.0
- **Pydantic:** 2.12.4
- **Starlette:** 1.3.1

### Utilities
- **PyYAML:** 6.0.3
- **Requests:** 2.34.2
- **tqdm:** 4.67.1
- **Matplotlib:** 3.10.6
- **Seaborn:** 0.13.2
- **Pytest:** 8.4.2

---

## ❌ Failed Installation

### DataSAIL
**Error:** `No matching distribution found for datasail`  
**Root Cause:** DataSAIL requires Python <3.13.0, incompatible with current Python 3.13  
**Affected Versions:**
- 1.1.3, 1.2.x, 1.3.0 all require Python >=3.9,<3.13.0

**Impact:**
- Cannot use DataSAIL for train/test splitting with molecular similarity clustering
- OligoGym benchmarking integration delayed

---

## 🔧 Workarounds for DataSAIL

### Option 1: Create Python 3.12 Conda Environment (RECOMMENDED)
```bash
conda create -n oligotox_py312 python=3.12 -y
conda activate oligotox_py312
pip install datasail
# Re-install all other packages in this env
```

### Option 2: Implement Custom Splitting Logic
```python
# src/data/custom_datasail_alternative.py
from sklearn.cluster import AgglomerativeClustering
from rdkit.Chem import AllChem, DataStructs
import numpy as np

def similarity_aware_split(
    sequences: list[str],
    test_size: float = 0.2,
    similarity_threshold: float = 0.7
) -> tuple[np.ndarray, np.ndarray]:
    """
    Split oligonucleotides ensuring test set has low similarity to train.
    Mimics DataSAIL behavior using RDKit molecular fingerprints.
    """
    # Compute pairwise Tanimoto similarity
    fps = [AllChem.GetMorganFingerprintAsBitVect(seq_to_mol(s), 2) for s in sequences]
    n = len(fps)
    similarity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            similarity_matrix[i, j] = similarity_matrix[j, i] = sim
    
    # Cluster sequences by similarity
    distance_matrix = 1 - similarity_matrix
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=1-similarity_threshold,
        metric='precomputed',
        linkage='average'
    ).fit(distance_matrix)
    
    # Assign clusters to train/test ensuring no overlap
    clusters = clustering.labels_
    unique_clusters = np.unique(clusters)
    np.random.shuffle(unique_clusters)
    n_test_clusters = int(len(unique_clusters) * test_size)
    test_clusters = set(unique_clusters[:n_test_clusters])
    
    test_idx = np.array([i for i, c in enumerate(clusters) if c in test_clusters])
    train_idx = np.array([i for i, c in enumerate(clusters) if c not in test_clusters])
    
    return train_idx, test_idx
```

### Option 3: Wait for DataSAIL Python 3.13 Support
Monitor: https://github.com/kalininalab/DataSAIL/issues

---

## 📋 Integration Readiness

| Module | Status | Blocker |
|--------|--------|--------|
| **OligoFormer** | ⚠️ Partial | Need model weights |
| **RNA-FM** | ✅ Ready | Transformers installed |
| **DataSAIL** | ❌ Blocked | Python version |
| **Hepatotoxicity** | ✅ Ready | Can proceed |
| **RNase H1** | ✅ Ready | Can proceed |
| **TLR7/8** | ✅ Ready | Can proceed |
| **GalNAc PK/PD** | ✅ Ready | SciPy/NumPy available |
| **PBPK** | ✅ Ready | Can proceed |
| **Geneformer** | ⚠️ Partial | Need model weights |
| **PyG Models** | ✅ Ready | torch-geometric installed |

---

## 🎯 Next Steps

1. **Immediate:** Document current platform model training pipeline (DONE)
2. **Priority 1:** Create Python 3.12 environment for DataSAIL
3. **Priority 2:** Download RNA-FM and OligoFormer pre-trained weights
4. **Priority 3:** Implement hepatotoxicity/KC activation/TLR modules
5. **Priority 4:** Integrate GalNAc PK/PD calculator
6. **Priority 5:** Build OligoGraph representation layer

---

## 📦 Full Package Inventory
