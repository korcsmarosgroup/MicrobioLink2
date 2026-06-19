#!/usr/bin/env python

"""Predict host-microbe interactions from domain-motif matches."""

from __future__ import annotations

import argparse
from pathlib import Path

from microbiolink_api.dmi import BidirectionalDomainMotifInteraction
from microbiolink_api.dmi import DomainMotifInteraction
from microbiolink_api.dmi import extract_uniprot_id
from microbiolink_api.dmi import predict_bidirectional_domain_motif_interactions
from microbiolink_api.dmi import predict_domain_motif_interactions
from microbiolink_api.dmi import predict_reverse_domain_motif_interactions
from microbiolink_api.dmi import read_fasta_sequences


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for standalone execution."""

    parser = argparse.ArgumentParser(
        description = (
            'Predict interactions between human and microbial proteins based '
            'on domain-motif interactions.'
        ),
    )
    parser.add_argument(
        '--mode',
        choices = ['forward', 'reverse', 'both'],
        default = 'forward',
        help = 'Prediction direction. Defaults to forward.',
    )
    parser.add_argument(
        '-fasta',
        '--fasta_file',
        help = 'Path to the human protein FASTA file for forward mode.',
    )
    parser.add_argument(
        '--bacterial_fasta_file',
        help = 'Path to the bacterial protein FASTA file for reverse mode.',
    )
    parser.add_argument(
        '-motif',
        '--elm_regex_file',
        help = 'Optional path to a motif regex file.',
    )
    parser.add_argument(
        '-interaction',
        '--motif_domain_file',
        help = 'Optional path to a motif-domain interaction file.',
    )
    parser.add_argument(
        '-domain',
        '--bacterial_domain_file',
        help = 'Path to the bacterial protein domain file for forward mode.',
    )
    parser.add_argument(
        '--human_domain_file',
        help = 'Path to the human protein domain file for reverse mode.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required = True,
        help = 'Path to the output file.',
    )
    return parser.parse_args()


def _require_mode_arguments(
    mode: str,
    args: argparse.Namespace,
) -> None:
    """Validate the required arguments for the chosen mode."""

    required_by_mode = {
        'forward': ['fasta_file', 'bacterial_domain_file'],
        'reverse': ['bacterial_fasta_file', 'human_domain_file'],
        'both': [
            'fasta_file',
            'bacterial_domain_file',
            'bacterial_fasta_file',
            'human_domain_file',
        ],
    }

    missing_arguments = [
        argument_name
        for argument_name in required_by_mode[mode]
        if not getattr(args, argument_name, None)
    ]

    if missing_arguments:
        joined_arguments = ', '.join(f'`{argument_name}`' for argument_name in missing_arguments)
        raise ValueError(
            f'Mode `{mode}` requires the following arguments: '
            f'{joined_arguments}.',
        )


def _write_forward_interactions(
    interactions: list[DomainMotifInteraction],
    output_file: str | Path,
) -> None:
    """Write forward-mode interactions in the legacy output format."""

    with open(output_file, 'w', encoding = 'utf-8') as outfile:
        outfile.write(
            '# Human Protein;Motif;Start;End;Bacterial domain;Bacteria Protein\n',
        )

        for interaction in interactions:
            outfile.write(
                f'{interaction.human_protein};{interaction.motif};'
                f'{interaction.start};{interaction.end};'
                f'{interaction.bacterial_domain};'
                f'{interaction.bacterial_protein}\n',
            )


def _write_bidirectional_interactions(
    interactions: list[BidirectionalDomainMotifInteraction],
    output_file: str | Path,
) -> None:
    """Write reverse or bidirectional interactions in an explicit schema."""

    with open(output_file, 'w', encoding = 'utf-8') as outfile:
        outfile.write(
            '# Host Protein;Microbial Protein;Motif;Start;End;Domain;'
            'Motif Protein Side;Domain Protein Side;Resource\n',
        )

        for interaction in interactions:
            outfile.write(
                f'{interaction.host_protein};{interaction.microbial_protein};'
                f'{interaction.motif};{interaction.start};{interaction.end};'
                f'{interaction.domain};{interaction.motif_protein_side};'
                f'{interaction.domain_protein_side};{interaction.resource}\n',
            )


def main(args: argparse.Namespace) -> None:
    """Run the domain-motif interaction workflow."""

    _require_mode_arguments(args.mode, args)

    if args.mode == 'forward':
        interactions = predict_domain_motif_interactions(
            fasta_file = args.fasta_file,
            bacterial_domain_file = args.bacterial_domain_file,
            elm_regex_file = args.elm_regex_file,
            motif_domain_file = args.motif_domain_file,
        )
        _write_forward_interactions(interactions, args.output_file)
        return

    if args.mode == 'reverse':
        interactions = predict_reverse_domain_motif_interactions(
            bacterial_fasta_file = args.bacterial_fasta_file,
            human_domain_file = args.human_domain_file,
            elm_regex_file = args.elm_regex_file,
            motif_domain_file = args.motif_domain_file,
        )
        _write_bidirectional_interactions(interactions, args.output_file)
        return

    interactions = predict_bidirectional_domain_motif_interactions(
        human_fasta_file = args.fasta_file,
        bacterial_domain_file = args.bacterial_domain_file,
        bacterial_fasta_file = args.bacterial_fasta_file,
        human_domain_file = args.human_domain_file,
        elm_regex_file = args.elm_regex_file,
        motif_domain_file = args.motif_domain_file,
        mode = 'both',
    )
    _write_bidirectional_interactions(interactions, args.output_file)


if __name__ == '__main__':
    main(parse_args())
