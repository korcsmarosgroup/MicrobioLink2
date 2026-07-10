#!/usr/bin/env python

"""Retrieve human protein FASTA sequences for expressed genes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from microbiolink.core.human_domains import fetch_protein_sequences
from microbiolink.core.human_domains import get_proteins


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description = 'Retrieve protein sequences for proteins from a gene list.',
    )
    parser.add_argument(
        '-genes',
        '--gene_expression',
        type = str,
        help = 'Path to the transcriptomics data',
    )
    parser.add_argument(
        '-id',
        '--id_type',
        choices = ['genesymbol', 'uniprot'],
        help = 'Type of gene identifier (genesymbol or uniprot)',
    )
    parser.add_argument(
        '-s',
        '--sep',
        help = 'Field separator in the protein list file',
    )
    parser.add_argument(
        '-lfl',
        '--location_filter_list',
        nargs = '+',
        default = None,
        help = (
            'Location filter list (options: plasma_membrane_transmembrane '
            'and/or plasma_membrane_peripheral and/or secreted)'
        ),
    )
    parser.add_argument(
        '-of',
        '--output_folder',
        default = '.',
        help = 'Output folder for result files',
    )
    parser.add_argument(
        '-oseq',
        '--output_sequences',
        default = 'protein_sequences.fasta',
        help = 'Output file for protein sequences',
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Run the human FASTA retrieval workflow.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Exit code.
    """
    args = parse_args(argv)

    proteins = get_proteins(
        args.gene_expression,
        args.id_type,
        args.sep,
        args.location_filter_list,
        args.output_folder,
    )
    print(proteins)

    batch_size = 100
    fasta_sequences = []

    for i in range(0, len(proteins), batch_size):
        batch_ids = proteins[i : i + batch_size]
        fasta_sequences.extend(fetch_protein_sequences(batch_ids))

    output_path = Path(args.output_folder) / args.output_sequences
    with open(output_path, 'w') as fasta_file:
        fasta_file.write(''.join(fasta_sequences))

    print(f'Protein sequences saved to {args.output_sequences}')

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
