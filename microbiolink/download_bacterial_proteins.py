#!/usr/bin/env python

"""Download bacterial proteins or proteomes from UniProt."""

from __future__ import annotations

import argparse
import sys

from microbiolink.core.uniprot import UNIPROT_BATCH_SIZE
from microbiolink.core.uniprot import download_protein_list_with_fields
from microbiolink.core.uniprot import download_proteome_with_fields
from microbiolink.core.uniprot import read_ids


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--id_list',
        required = True,
        help = 'Path to a file describing UniProt or proteome identifiers.',
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
        help = 'Path to the output file.',
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Run the UniProt download workflow.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Exit code.
    """
    args = parse_args(argv)
    ids = read_ids(args.id_list, args.sep, args.id_column)
    header_written = False

    with open(args.output, 'w', encoding = 'utf-8') as output_file:
        for start in range(0, len(ids), UNIPROT_BATCH_SIZE):
            batch = ids[start : start + UNIPROT_BATCH_SIZE]

            if args.id_type == 'UP':
                for proteome_id in batch:
                    lines = download_proteome_with_fields(proteome_id).rstrip().split('\n')

                    if not lines or not lines[0]:
                        continue

                    if not header_written:
                        output_file.write(lines[0] + '\tProteome_ID\n')
                        header_written = True

                    for row in lines[1:]:
                        if row:
                            output_file.write(f'{row}\t{proteome_id}\n')
                continue

            lines = download_protein_list_with_fields(batch).rstrip().split('\n')

            if not lines or not lines[0]:
                continue

            if not header_written:
                output_file.write(lines[0] + '\n')
                header_written = True

            for row in lines[1:]:
                if row:
                    output_file.write(row + '\n')

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
