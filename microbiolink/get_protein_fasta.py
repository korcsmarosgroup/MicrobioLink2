#!/usr/bin/env python

"""Download FASTA sequences from UniProt for any list of protein accessions."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from microbiolink.core.uniprot import DEFAULT_BATCH_SIZE
from microbiolink.core.uniprot import fetch_fasta_sequences
from microbiolink.core.uniprot import read_ids


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description = 'Download FASTA sequences for a list of UniProt accessions.',
    )
    parser.add_argument(
        '--id_list',
        required = True,
        help = 'Path to a file containing UniProt accessions.',
    )
    parser.add_argument(
        '--sep',
        required = True,
        help = 'Field separator in the identifier file.',
    )
    parser.add_argument(
        '--id_column',
        type = int,
        required = True,
        help = 'One-based column number containing the identifiers.',
    )
    parser.add_argument(
        '--output',
        required = True,
        help = 'Path to the output FASTA file.',
    )
    parser.add_argument(
        '--batch_size',
        type = int,
        default = DEFAULT_BATCH_SIZE,
        help = f'Number of accessions per request (default: {DEFAULT_BATCH_SIZE}).',
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Run the FASTA download workflow.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Exit code.
    """
    args = parse_args(argv)
    ids = read_ids(args.id_list, args.sep, args.id_column)

    with open(Path(args.output), 'w', encoding = 'utf-8') as output_file:
        for start in range(0, len(ids), args.batch_size):
            batch = ids[start : start + args.batch_size]
            output_file.write(fetch_fasta_sequences(batch))

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
