#!/usr/bin/env python

"""Predict host-microbe interactions from bacterial motifs matching human domains."""

from __future__ import annotations

import argparse
from pathlib import Path

from microbiolink.core.dmi import predict_reverse_domain_motif_interactions
from microbiolink.core.dmi import resolve_dmi_resource_bundle_by_name


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description = (
            'Predict interactions between microbial and human proteins based on '
            'domain-motif interactions (reversed: bacterial motifs, human domains).'
        ),
    )
    parser.add_argument(
        '-fasta',
        '--fasta_file',
        required = True,
        help = 'Path to the bacterial protein FASTA file.',
    )
    parser.add_argument(
        '-motif',
        '--elm_regex_file',
        required = False,
        default = None,
        help = 'Path to the ELM regex file. Omit to use the packaged '
        'resource set selected by --resource_set.',
    )
    parser.add_argument(
        '-interaction',
        '--motif_domain_file',
        required = False,
        default = None,
        help = 'Path to the motif-domain interaction file. Omit to use the '
        'packaged resource set selected by --resource_set.',
    )
    parser.add_argument(
        '-domain',
        '--human_domain_file',
        required = True,
        help = 'Path to the human protein domain file.',
    )
    parser.add_argument(
        '--resource_set',
        choices = ['default', 'elm', '3did'],
        default = 'default',
        help = 'Packaged DMI resource set to use when --motif/--interaction '
        'are omitted.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required = True,
        help = 'Path to the output file.',
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    """Run the reversed domain-motif interaction workflow.

    Args:
        args: Parsed argument namespace from parse_args().
    """
    resource_bundle = None
    if args.elm_regex_file is None or args.motif_domain_file is None:
        resource_bundle = resolve_dmi_resource_bundle_by_name(args.resource_set)

    interactions = predict_reverse_domain_motif_interactions(
        bacterial_fasta_file = args.fasta_file,
        human_domain_file = args.human_domain_file,
        elm_regex_file = args.elm_regex_file,
        motif_domain_file = args.motif_domain_file,
        resource_bundle = resource_bundle,
    )

    with open(Path(args.output_file), 'w', encoding = 'utf-8') as output_file:
        output_file.write(
            '# Bacterial Protein;Motif;Start;End;Human Domain;Human Protein\n',
        )
        for interaction in interactions:
            output_file.write(
                f'{interaction.microbial_protein};{interaction.motif};'
                f'{interaction.start};{interaction.end};'
                f'{interaction.domain};{interaction.host_protein}\n',
            )


if __name__ == '__main__':
    main(parse_args())
