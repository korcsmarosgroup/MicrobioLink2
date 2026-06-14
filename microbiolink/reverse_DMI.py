#!/usr/bin/env python

"""Predict host-microbe interactions from bacterial motifs matching human domains."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from microbiolink.DMI import (
    create_uniprot_motif_dict,
    parse_elm_regex,
    parse_motif_domain,
    parse_protein_domain,
    read_fasta_sequences,
)


@dataclass(frozen=True)
class ReverseDomainMotifInteraction:
    """Represent one predicted reverse domain-motif interaction."""

    bacterial_protein: str
    motif: str
    start: int
    end: int
    human_domain: str
    human_protein: str


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


def predict_reverse_domain_motif_interactions_from_data(
    bacterial_sequences: dict[str, str],
    elm_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    human_domain_table: dict[str, list[str]],
) -> list[ReverseDomainMotifInteraction]:
    """Predict reverse domain-motif interactions from in-memory inputs.

    Scans bacterial protein sequences for ELM motif patterns, then maps each
    hit through the motif-domain interaction table to the human proteins that
    carry a compatible Pfam domain.

    Args:
        bacterial_sequences: Mapping from FASTA header to bacterial protein
            sequence.
        elm_regex: Mapping from ELM motif identifier to regex pattern.
            Pass through filter_cleavage_motifs() first to exclude CLV_ entries.
        motif_domains: Mapping from ELM motif identifier to compatible Pfam
            domain identifiers.
        human_domain_table: Mapping from Pfam domain identifier to human
            UniProt accessions carrying that domain.

    Returns:
        List of predicted interactions, one entry per
        (bacterial protein, motif hit, human domain, human protein) combination.
    """
    bacterial_motif = create_uniprot_motif_dict(bacterial_sequences, elm_regex)
    interactions: list[ReverseDomainMotifInteraction] = []

    for motif_name, compatible_domains in motif_domains.items():
        motif_hits = [
            (bacterial_id, int(start), int(end))
            for bacterial_id, matches in bacterial_motif.items()
            for match_name, start, end in matches
            if match_name == motif_name
        ]

        for domain in compatible_domains:
            for human_protein in human_domain_table.get(domain, []):
                for bacterial_id, start, end in motif_hits:
                    interactions.append(
                        ReverseDomainMotifInteraction(
                            bacterial_protein=bacterial_id,
                            motif=motif_name,
                            start=start,
                            end=end,
                            human_domain=domain,
                            human_protein=human_protein,
                        ),
                    )

    return interactions


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
    elm_regex = filter_cleavage_motifs(parse_elm_regex(args.elm_regex_file))
    motif_domain = parse_motif_domain(args.motif_domain_file)
    human_domain_table = parse_protein_domain(args.human_domain_file)

    interactions = predict_reverse_domain_motif_interactions_from_data(
        bacterial_sequences=bacterial_sequences,
        elm_regex=elm_regex,
        motif_domains=motif_domain,
        human_domain_table=human_domain_table,
    )

    with open(Path(args.output_file), 'w', encoding='utf-8') as output_file:
        output_file.write(
            '# Bacterial Protein;Motif;Start;End;Human Domain;Human Protein\n',
        )
        for interaction in interactions:
            output_file.write(
                f'{interaction.bacterial_protein};{interaction.motif};'
                f'{interaction.start};{interaction.end};'
                f'{interaction.human_domain};{interaction.human_protein}\n',
            )


if __name__ == '__main__':
    main(parse_args())
