#!/usr/bin/env python

"""Expression filtering utilities for the library-style MicrobioLink API."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg
from scipy.stats import gaussian_kde


PathLike = str | Path


def read_count_matrix(
    filename: PathLike,
    index_col: int | str = 0,
) -> pd.DataFrame:
    """Read a host count matrix from disk.

    Args:
        filename: Path to a comma-separated count matrix.
        index_col: Column to use as the row index.

    Returns:
        A pandas data frame indexed by gene or protein identifier.
    """

    return pd.read_csv(filename, index_col = index_col)


def _filter_numeric_series(
    values: pd.Series,
    zscore_threshold: float,
) -> pd.Series:
    """Apply the legacy z-score expression filter to one sample column.

    The original CLI implementation fits a Gaussian kernel density estimate to
    the observed values, finds the modal expression value, and keeps only
    measurements with z-scores above the requested threshold. When there is not
    enough variance to fit a KDE reliably, the function keeps the observed
    values as-is.

    Args:
        values: One expression column.
        zscore_threshold: Minimum z-score to retain a value.

    Returns:
        A filtered numeric series with dropped values represented as NaN.
    """

    numeric_values = pd.to_numeric(values, errors = 'coerce')
    filtered_values = numeric_values.dropna()

    if filtered_values.size < 2:
        return numeric_values

    if np.isclose(filtered_values.max(), filtered_values.min()):
        return numeric_values

    try:
        kernel = gaussian_kde(filtered_values.to_numpy(dtype = float))
    except (linalg.LinAlgError, ValueError):
        return numeric_values

    xi = np.linspace(filtered_values.min(), filtered_values.max(), 100)
    yi = kernel.evaluate(xi)

    mu = float(xi[np.argmax(yi)])
    upper_tail = filtered_values[filtered_values > mu]

    if upper_tail.empty:
        return numeric_values

    sigma = float((upper_tail.mean() - mu) * np.sqrt(np.pi / 2))
    if np.isclose(sigma, 0.0):
        return numeric_values

    zscores = (numeric_values - mu) / sigma
    return numeric_values.where(zscores > zscore_threshold)


def filter_counts_by_zscore(
    count_matrix: pd.DataFrame,
    zscore_threshold: float = -3,
) -> pd.DataFrame:
    """Filter a count matrix by the MicrobioLink z-score rule.

    Args:
        count_matrix: Input count matrix with genes or proteins on rows and
            samples on columns.
        zscore_threshold: Minimum z-score to retain a value.

    Returns:
        A filtered copy of the count matrix.
    """

    filtered_matrix = pd.DataFrame(index = count_matrix.index)

    for column_name in count_matrix.columns:
        filtered_matrix[column_name] = _filter_numeric_series(
            count_matrix[column_name],
            zscore_threshold = zscore_threshold,
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

    count_matrix = read_count_matrix(input_file)
    filtered_matrix = filter_counts_by_zscore(
        count_matrix,
        zscore_threshold = zscore_threshold,
    )

    if output_file is not None:
        filtered_matrix.to_csv(output_file)

    return filtered_matrix
