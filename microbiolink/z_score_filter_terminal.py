#!/usr/bin/env python

"""Gene expression filtering based on z-score thresholding."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def main(args: argparse.Namespace) -> None:
    df = pd.read_csv(args.input_file, index_col=0)

    columns = df.columns
    rows = df.index.values

    df_cells_log2_filtered = pd.DataFrame(columns=columns, index=rows)

    for column in columns:
        count = np.array(df[column].tolist(), dtype=float)
        count_filtered = count[np.logical_not(np.isnan(count))]

        kernel = gaussian_kde(count_filtered)
        xi = np.linspace(count_filtered.min(), count_filtered.max(), 100)
        yi = kernel.evaluate(xi)

        mu = xi[np.argmax(yi)]
        U = count_filtered[count_filtered > mu].mean()
        sigma = (U - mu) * np.sqrt(np.pi / 2)

        zcount = (count - mu) / sigma

        score_list = [
            count[list(zcount).index(value)] if value > args.zscore else 'NaN'
            for value in zcount
        ]
        df_cells_log2_filtered[column] = pd.Series(score_list, index=rows)

    df_cells_log2_filtered.to_csv(args.output_file)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Gene expression filtration based on individual cell count and z-score.',
    )
    parser.add_argument('-i', '--input_file', required=True, help='Input CSV file with gene expression data.')
    parser.add_argument(
        '-zscore',
        '--zscore',
        required=True,
        type=int,
        default=-3,
        help='Z-score cut-off to filter lowly expressed genes',
    )
    parser.add_argument('-o', '--output_file', required=True, help='Output CSV file for filtered results.')
    return parser.parse_args()


def cli_main() -> None:
    main(parse_args())


if __name__ == '__main__':
    cli_main()
