#!/usr/bin/env python

"""Predict host-microbe interactions from domain-motif matches."""

from __future__ import annotations

import argparse

from microbiolink.core.dmi import predict_domain_motif_interactions
from microbiolink.core.dmi import resolve_dmi_resource_bundle_by_name


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for standalone execution."""

    parser = argparse.ArgumentParser(
        description = (
            'Predict interactions between human and microbial proteins based '
            'on domain-motif interactions.'
        ),
    )
    parser.add_argument(
        '-fasta',
        '--fasta_file',
        required = True,
        help = 'Path to the human protein FASTA file.',
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
        '--bacterial_domain_file',
        required = True,
        help = 'Path to the bacterial protein domain file.',
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
    """Run the domain-motif interaction workflow."""

    resource_bundle = None
    if args.elm_regex_file is None or args.motif_domain_file is None:
        resource_bundle = resolve_dmi_resource_bundle_by_name(args.resource_set)

    interactions = predict_domain_motif_interactions(
        fasta_file = args.fasta_file,
        bacterial_domain_file = args.bacterial_domain_file,
        elm_regex_file = args.elm_regex_file,
        motif_domain_file = args.motif_domain_file,
        resource_bundle = resource_bundle,
    )

    with open(args.output_file, 'w', encoding = 'utf-8') as output_file:
        output_file.write(
            '# Human Protein;Motif;Start;End;Bacterial domain;Bacteria Protein\n',
        )

        for interaction in interactions:
            output_file.write(
                f'{interaction.human_protein};{interaction.motif};'
                f'{interaction.start};{interaction.end};'
                f'{interaction.bacterial_domain};{interaction.bacterial_protein}\n',
            )


if __name__ == '__main__':
    main(parse_args())
