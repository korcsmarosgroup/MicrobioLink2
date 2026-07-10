#!/usr/bin/env python

"""Download Pfam domain annotations for expressed human proteins."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from microbiolink.core.human_domains import _extract_swissprot_ids
from microbiolink.core.human_domains import read_expressed_genes
from microbiolink.core.human_domains import translate_symbol_to_uniprot
from microbiolink.core.uniprot import UNIPROT_BATCH_SIZE
from microbiolink.core.uniprot import download_protein_list_with_fields


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description = 'Download Pfam domain annotations for expressed human proteins.',
    )
    parser.add_argument(
        '--gene_expression',
        required = True,
        help = 'Path to z-score filtered gene expression CSV.',
    )
    parser.add_argument(
        '--id_type',
        choices = ['genesymbol', 'uniprot'],
        required = True,
        help = 'Type of gene identifier: genesymbol or uniprot.',
    )
    parser.add_argument(
        '--sep',
        required = True,
        help = 'Field separator in the gene expression file.',
    )
    parser.add_argument(
        '--output',
        required = True,
        help = 'Path to the output TSV file.',
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Run the human domain download workflow.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Exit code.
    """
    args = parse_args(argv)
    genes = read_expressed_genes(args.gene_expression, args.sep)

    if args.id_type == 'genesymbol':
        translation_dict = translate_symbol_to_uniprot(genes)
        uniprot_ids = _extract_swissprot_ids(translation_dict)
    else:
        uniprot_ids = genes

    header_written = False

    with open(Path(args.output), 'w', encoding = 'utf-8') as output_file:
        for start in range(0, len(uniprot_ids), UNIPROT_BATCH_SIZE):
            batch = uniprot_ids[start : start + UNIPROT_BATCH_SIZE]
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
