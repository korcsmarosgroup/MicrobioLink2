"""Helpers for reading a gene count matrix and deriving protein identifiers from it."""

from pathlib import Path
from typing import Union

import pandas as pd

PathLike = Union[str, Path]


def read_count_matrix(
    filename: PathLike,
    index_col: int | str = 0,
) -> pd.DataFrame:
    """Read a gene/protein count matrix from disk.

    Args:
        filename: Path to a comma-separated count matrix.
        index_col: Column to use as the row index.

    Returns:
        A pandas data frame indexed by gene or protein identifier.
    """

    return pd.read_csv(filename, index_col=index_col)


def extract_gene_symbols_from_count_matrix(count_matrix: pd.DataFrame) -> list[str]:
    """Extract the gene symbol row index from a gene count matrix.

    Args:
        count_matrix: A count matrix indexed by gene symbol (e.g. the output
            of `read_count_matrix`).

    Returns:
        The row index as a plain list of gene symbol strings.
    """

    return count_matrix.index.tolist()


def extract_expressed_gene_symbols_from_count_matrix(count_matrix: pd.DataFrame) -> list[str]:
    """Extract gene symbols with at least one non-zero, non-NaN expression value.

    A safety-net filter for when the count matrix hasn't already been
    filtered by z-score (see `zscore_filter.filter_counts_by_zscore`).

    Args:
        count_matrix: A count matrix indexed by gene symbol (e.g. the output
            of `read_count_matrix`).

    Returns:
        The row index, restricted to genes with at least one expressed
        value, as a plain list of gene symbol strings.
    """

    expressed_mask = count_matrix.notna().any(axis=1) & (count_matrix != 0).any(axis=1)
    return count_matrix.index[expressed_mask].tolist()
