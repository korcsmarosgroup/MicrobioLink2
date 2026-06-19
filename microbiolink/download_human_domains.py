#!/usr/bin/env python

"""Download Pfam domain annotations for expressed human proteins."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from microbiolink.download_protein_domains import (
    UNIPROT_BATCH_SIZE,
    download_protein_list_with_fields,
)


def _is_expressed(value: str) -> bool:
    """Return True if an expression value is non-NaN and non-zero."""
    if value in ('NaN', '', 'nan'):
        return False
    try:
        return float(value) != 0.0
    except ValueError:
        return False


def read_expressed_genes(
    filename: str | Path,
    separator: str,
) -> list[str]:
    """Read expressed gene or protein IDs from a z-score filtered CSV.

    Includes any row where at least one expression column is non-NaN and
    non-zero, consistent with z_score_filter_terminal.py output format.

    Args:
        filename: Path to the z-score filtered CSV file.
        separator: Field separator used in the file.

    Returns:
        List of gene or protein identifiers for expressed entries.
    """
    genes: list[str] = []

    with open(Path(filename), encoding='utf-8-sig') as gene_file:
        next(gene_file, None)

        for line in gene_file:
            fields = line.strip().split(separator)
            if len(fields) < 2:
                continue
            gene = fields[0]
            if any(_is_expressed(v) for v in fields[1:]):
                genes.append(gene)

    return genes


def _extract_swissprot_ids(translation_dict: dict) -> list[str]:
    """Extract Swiss-Prot UniProt accessions from a MyGene translation dict.

    Args:
        translation_dict: Mapping of gene symbol to MyGene uniprot entry.

    Returns:
        Flat list of Swiss-Prot UniProt accession strings.
    """
    proteins: list[str] = []
    for entry in translation_dict.values():
        if 'Swiss-Prot' not in entry:
            continue
        value = entry['Swiss-Prot']
        if isinstance(value, str):
            proteins.append(value)
        elif isinstance(value, list):
            proteins.extend(value)
    return proteins


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description='Download Pfam domain annotations for expressed human proteins.',
    )
    parser.add_argument(
        '--gene_expression',
        required=True,
        help='Path to z-score filtered gene expression CSV.',
    )
    parser.add_argument(
        '--id_type',
        choices=['genesymbol', 'uniprot'],
        required=True,
        help='Type of gene identifier: genesymbol or uniprot.',
    )
    parser.add_argument(
        '--sep',
        required=True,
        help='Field separator in the gene expression file.',
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to the output TSV file.',
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
        from microbiolink.get_human_fasta import translate_symbol_to_uniprot
        translation_dict = translate_symbol_to_uniprot(genes)
        uniprot_ids = _extract_swissprot_ids(translation_dict)
    else:
        uniprot_ids = genes

    header_written = False

    with open(Path(args.output), 'w', encoding='utf-8') as output_file:
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
