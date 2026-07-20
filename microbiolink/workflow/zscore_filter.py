"""Z-score based filtering of gene/protein count matrices."""

from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from ..utils import gene_matrix

PathLike = Union[str, Path]


def _filter_column_by_zscore(
    values: pd.Series,
    zscore_threshold: float,
) -> pd.Series:
    """Apply the MicrobioLink z-score filter to one sample column.

    Fits a Gaussian kernel density estimate to the non-NaN values, takes the
    mode of the KDE as mu, derives sigma from the mean of values above mu,
    and keeps only values whose z-score is strictly above the threshold.

    Args:
        values: One expression column.
        zscore_threshold: Minimum z-score to retain a value.

    Returns:
        A filtered numeric series with dropped values represented as NaN.
    """

    count = values.to_numpy(dtype=float)
    count_filtered = count[np.logical_not(np.isnan(count))]

    kernel = gaussian_kde(count_filtered)

    xi = np.linspace(count_filtered.min(), count_filtered.max(), 100)
    yi = kernel.evaluate(xi)
    mu = xi[np.argmax(yi)]

    upper_tail = count_filtered[count_filtered > mu]
    sigma = (upper_tail.mean() - mu) * np.sqrt(np.pi / 2)

    zcount = (count - mu) / sigma
    kept = np.where(zcount > zscore_threshold, count, np.nan)

    return pd.Series(kept, index=values.index)


def filter_counts_by_zscore(
    count_matrix: pd.DataFrame,
    zscore_threshold: float = -3,
) -> pd.DataFrame:
    """Filter a count matrix by the MicrobioLink z-score rule.

    Each column (sample/condition) is filtered independently: a value is
    kept as-is if its column-wise z-score is strictly above
    zscore_threshold, otherwise it becomes NaN.

    Args:
        count_matrix: Input count matrix with genes or proteins on rows and
            samples/conditions on columns.
        zscore_threshold: Minimum z-score to retain a value.

    Returns:
        A filtered copy of the count matrix, same shape/index/columns as the
        input.
    """

    filtered_matrix = pd.DataFrame(index=count_matrix.index)

    for column_name in count_matrix.columns:
        filtered_matrix[column_name] = _filter_column_by_zscore(
            count_matrix[column_name],
            zscore_threshold=zscore_threshold,
        )

    return filtered_matrix


def filter_count_matrix_file(
    input_file: PathLike,
    zscore_threshold: float = -3,
    output_file: PathLike | None = None,
) -> pd.DataFrame:
    """Read, filter, and optionally write a count matrix.

    Args:
        input_file: Path to the input count matrix.
        zscore_threshold: Minimum z-score to retain a value.
        output_file: Optional path to write the filtered matrix.

    Returns:
        The filtered count matrix.
    """

    count_matrix = gene_matrix.read_count_matrix(input_file)
    filtered_matrix = filter_counts_by_zscore(
        count_matrix,
        zscore_threshold=zscore_threshold,
    )

    if output_file is not None:
        filtered_matrix.to_csv(output_file)

    return filtered_matrix
