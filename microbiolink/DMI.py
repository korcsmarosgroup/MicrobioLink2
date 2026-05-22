#!/usr/bin/env python

"""Predict host-microbe interactions from domain-motif matches."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def extract_uniprot_id(fasta_header: str) -> str:
    """Extract the UniProt identifier from a FASTA header."""

    fields = fasta_header.split('|')
    if len(fields) < 2:
        raise ValueError(
            f'FASTA header does not contain a UniProt identifier: {fasta_header}',
        )
    return fields[1]


def read_fasta_sequences(filename: str | Path) -> dict[str, str]:
    """Read a FASTA file into a header-to-sequence mapping."""

    sequences: dict[str, str] = {}
    current_header: str | None = None
    current_fragments: list[str] = []

    with open(filename, encoding='utf-8') as fasta_file:
        for raw_line in fasta_file:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_fragments)

                current_header = line[1:]
                current_fragments = []
                continue

            current_fragments.append(line)

    if current_header is not None:
        sequences[current_header] = ''.join(current_fragments)

    return sequences


def parse_elm_regex(filename: str | Path) -> dict[str, str]:
    """Parse ELM motif definitions from a tab-separated file."""

    elm_regex: dict[str, str] = {}

    with open(filename, encoding='utf-8') as motif_table:
        next(motif_table, None)

        for line in motif_table:
            if not line or line.startswith('#'):
                continue

            fields = line.replace('"', '').strip().split('\t')
            if len(fields) > 4:
                elm_regex[fields[1]] = fields[4]

    return elm_regex


def parse_motif_domain(filename: str | Path) -> dict[str, list[str]]:
    """Parse motif-to-domain mappings."""

    motif_domain: dict[str, list[str]] = {}

    with open(filename, encoding='utf-8') as motif_domain_table:
        next(motif_domain_table, None)

        for line in motif_domain_table:
            fields = line.replace('"', '').strip().split('\t')
            if len(fields) <= 1:
                continue

            motif_domain.setdefault(fields[0], []).append(fields[1])

    return motif_domain


def parse_protein_domain(filename: str | Path) -> dict[str, list[str]]:
    """Parse bacterial proteins indexed by Pfam domain."""

    pfam_uniprot: dict[str, list[str]] = {}

    with open(filename, encoding='utf-8') as protein_domain:
        next(protein_domain, None)

        for line in protein_domain:
            fields = line.strip().split('\t')
            if len(fields) <= 1:
                continue

            for pfam in fields[1].split(';'):
                pfam_uniprot.setdefault(pfam, []).append(fields[0])

    return pfam_uniprot


def create_uniprot_motif_dict(
    human_sequences: dict[str, str],
    elm_regex: dict[str, str],
) -> dict[str, list[tuple[str, str, str]]]:
    """Create a mapping of UniProt IDs to matched motifs."""

    uniprot_motif: dict[str, list[tuple[str, str, str]]] = {}

    for header, sequence in human_sequences.items():
        uniprot_id = extract_uniprot_id(header)
        matches: list[tuple[str, str, str]] = []

        for motif_name, motif_pattern in elm_regex.items():
            for match in re.finditer(motif_pattern, sequence):
                matches.append(
                    (
                        motif_name,
                        str(match.start()),
                        str(match.end()),
                    ),
                )

        if matches:
            uniprot_motif[uniprot_id] = matches

    return uniprot_motif


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for standalone execution."""

    parser = argparse.ArgumentParser(
        description=(
            'Predict interactions between human and microbial proteins based '
            'on domain-motif interactions.'
        ),
    )
    parser.add_argument(
        '-fasta',
        '--fasta_file',
        required=True,
        help='Path to the human protein FASTA file.',
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
        '--bacterial_domain_file',
        required=True,
        help='Path to the bacterial protein domain file.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required=True,
        help='Path to the output file.',
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    """Run the domain-motif interaction workflow."""

    human_sequences = read_fasta_sequences(args.fasta_file)
    elm_regex = parse_elm_regex(args.elm_regex_file)
    motif_domain = parse_motif_domain(args.motif_domain_file)
    pfam_uniprot = parse_protein_domain(args.bacterial_domain_file)
    uniprot_motif = create_uniprot_motif_dict(human_sequences, elm_regex)

    with open(args.output_file, 'w', encoding='utf-8') as output_file:
        output_file.write(
            '# Human Protein;Motif;Start;End;Bacterial domain;Bacteria Protein\n',
        )

        for motif_name, motif_domains in motif_domain.items():
            motif_hits = [
                (uniprot_id, start, end)
                for uniprot_id, matches in uniprot_motif.items()
                for match_name, start, end in matches
                if match_name == motif_name
            ]

            for domain in motif_domains:
                for bacterial_protein in pfam_uniprot.get(domain, []):
                    for uniprot_id, start, end in motif_hits:
                        output_file.write(
                            f'{uniprot_id};{motif_name};{start};{end};'
                            f'{domain};{bacterial_protein}\n',
                        )


if __name__ == '__main__':
    main(parse_args())
