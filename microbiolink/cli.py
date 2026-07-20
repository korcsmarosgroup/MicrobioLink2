"""Console-script entry points for MicrobioLink. Argparse boilerplate only."""

import argparse


def zscore_filter() -> int:
    """Filter a gene/protein count matrix by a z-score cutoff."""

    from . import zscore_filter as zscore_filter_module

    parser = argparse.ArgumentParser(
        description='Gene expression filtration based on individual cell count and z-score.',
    )
    parser.add_argument(
        '-i',
        '--input_file',
        required=True,
        help='Input CSV file with gene expression data.',
    )
    parser.add_argument(
        '-zscore',
        '--zscore',
        required=True,
        type=float,
        help='Z-score cut-off to filter lowly expressed genes.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required=True,
        help='Output CSV file for filtered results.',
    )
    args = parser.parse_args()

    zscore_filter_module.filter_count_matrix_file(
        args.input_file,
        zscore_threshold=args.zscore,
        output_file=args.output_file,
    )
    return 0
