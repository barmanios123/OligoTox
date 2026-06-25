import os
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.omics.column_mapper import standardize_deg_columns


def test_column_mapper_handles_common_deg_names():
    df = pd.DataFrame(
        {
            "Gene Symbol": ["ISG15", "OAS1"],
            "log2FoldChange": [2.1, 1.8],
            "padj": [1e-6, 2e-5],
            "baseMean": [120, 100],
        }
    )

    out = standardize_deg_columns(df)
    assert "gene_symbol" in out.columns
    assert "log2fc" in out.columns
    assert "pval_adj" in out.columns
    assert "base_mean" in out.columns

    print("Mapped columns:", list(out.columns))


if __name__ == "__main__":
    test_column_mapper_handles_common_deg_names()
