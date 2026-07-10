#!/usr/bin/env python

"""Download FASTA sequences for all bacterial proteins in a proteome or accession list."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from microbiolink.core.uniprot import DEFAULT_BATCH_SIZE
from microbiolink.core.uniprot import fetch_fasta_sequences
from microbiolink.core.uniprot import fetch_proteome_fasta
from microbiolink.core.uniprot import read_ids


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description = (
            'Download FASTA sequences for all bacterial proteins in a proteome '
            'or accession list.'
        ),
    )
    parser.add_argument(
        '--id_list',
        required = True,
        help = 'Path to the same ID list used for download_bacterial_proteins.',
    )
    parser.add_argument(
        '--sep',
        required = True,
        help = 'Field separator in the identifier file.',
    )
    parser.add_argument(
        '--id_type',
        choices = ['uniprot', 'UP'],
        required = True,
        help = 'Identifier type: UniProt accessions or UniProt proteomes.',
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
        help = f'Accessions per request for uniprot mode (default: {DEFAULT_BATCH_SIZE}).',
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Run the bacterial FASTA download workflow.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Exit code.
    """
    args = parse_args(argv)
    ids = read_ids(args.id_list, args.sep, args.id_column)

    with open(Path(args.output), 'w', encoding = 'utf-8') as output_file:
        if args.id_type == 'UP':
            for proteome_id in ids:
                output_file.write(fetch_proteome_fasta(proteome_id))
        else:
            for start in range(0, len(ids), args.batch_size):
                batch = ids[start : start + args.batch_size]
                output_file.write(fetch_fasta_sequences(batch))

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
