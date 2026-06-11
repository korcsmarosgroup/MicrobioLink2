#!/usr/bin/env python

"""Predict host-microbe interactions from bacterial motifs matching human domains."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from microbiolink.DMI import (
    create_uniprot_motif_dict,
    parse_elm_regex,
    parse_motif_domain,
    parse_protein_domain,
    read_fasta_sequences,
)


def filter_cleavage_motifs(elm_regex: dict[str, str]) -> dict[str, str]:
    """Remove cleavage site motifs from an ELM regex dictionary.

    Cleavage (CLV_*) motifs mediate proteolytic digestion rather than
    domain binding and are not meaningful in a domain-motif interaction
    context.

    Args:
        elm_regex: Mapping of ELM identifier to regex pattern.

    Returns:
        Copy of the dict with all CLV_* entries removed.
    """
    return {
        name: pattern
        for name, pattern in elm_regex.items()
        if not name.startswith('CLV_')
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description=(
            'Predict interactions between microbial and human proteins based on '
            'domain-motif interactions (reversed: bacterial motifs, human domains).'
        ),
    )
    parser.add_argument(
        '-fasta',
        '--fasta_file',
        required=True,
        help='Path to the bacterial protein FASTA file.',
    )
    parser.add_argument(
        '-motif',
        '--elm_regex_file',
        required=True,
        help='Path to the ELM regex file.',
    )
    parser.add_argument(
        '-interaction',
        '--motif_domain_file',
        required=True,
        help='Path to the motif-domain interaction file.',
    )
    parser.add_argument(
        '-domain',
        '--human_domain_file',
        required=True,
        help='Path to the human protein domain file.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required=True,
        help='Path to the output file.',
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    """Run the reversed domain-motif interaction workflow.

    Args:
        args: Parsed argument namespace from parse_args().
    """
    bacterial_sequences = read_fasta_sequences(args.fasta_file)
    elm_regex = parse_elm_regex(args.elm_regex_file)
    elm_regex = filter_cleavage_motifs(elm_regex)
    motif_domain = parse_motif_domain(args.motif_domain_file)
    pfam_human = parse_protein_domain(args.human_domain_file)
    bacterial_motif = create_uniprot_motif_dict(bacterial_sequences, elm_regex)

    with open(Path(args.output_file), 'w', encoding='utf-8') as output_file:
        output_file.write(
            '# Bacterial Protein;Motif;Start;End;Human Domain;Human Protein\n',
        )

        for motif_name, motif_domains in motif_domain.items():
            motif_hits = [
                (bacterial_id, start, end)
                for bacterial_id, matches in bacterial_motif.items()
                for match_name, start, end in matches
                if match_name == motif_name
            ]

            for domain in motif_domains:
                for human_protein in pfam_human.get(domain, []):
                    for bacterial_id, start, end in motif_hits:
                        output_file.write(
                            f'{bacterial_id};{motif_name};{start};{end};'
                            f'{domain};{human_protein}\n',
                        )


if __name__ == '__main__':
    main(parse_args())
