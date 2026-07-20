"""Helpers for deriving protein identifiers from a gene count matrix."""

import pandas as pd


def extract_gene_symbols_from_count_matrix(count_matrix: pd.DataFrame) -> list[str]:
    """Extract the gene symbol row index from a gene count matrix.

    Args:
        count_matrix: A count matrix indexed by gene symbol (e.g. the output
            of `zscore_filter.read_count_matrix`).

    Returns:
        The row index as a plain list of gene symbol strings.
    """

    return count_matrix.index.tolist()
