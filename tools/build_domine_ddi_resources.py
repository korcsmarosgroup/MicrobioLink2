#!/usr/bin/env python

"""Build packaged DOMINE DDI resources from the v2 INTERACTION table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_domine_pairs(
    interaction_file: Path,
) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """Load all and high-confidence Pfam pairs from DOMINE."""

    all_pairs: set[tuple[str, str]] = set()
    high_confidence_pairs: set[tuple[str, str]] = set()

    with open(interaction_file, encoding = 'utf-8') as infile:
        for line in infile:
            fields = line.strip().split('|')

            if len(fields) < 18:
                continue

            pfam_pair = tuple(sorted((fields[0], fields[1])))
            all_pairs.add(pfam_pair)

            if fields[17] == 'HC':
                high_confidence_pairs.add(pfam_pair)

    return all_pairs, high_confidence_pairs


def write_pairs(
    pfam_pairs: set[tuple[str, str]],
    output_file: Path,
) -> None:
    """Write a canonical two-column Pfam-Pfam table."""

    with open(output_file, 'w', encoding = 'utf-8', newline = '') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')

        for pfam_a, pfam_b in sorted(pfam_pairs):
            writer.writerow([pfam_a, pfam_b])


def main() -> None:
    """Build packaged DOMINE DDI TSV resources."""

    parser = argparse.ArgumentParser(
        description = 'Build DOMINE Pfam-Pfam TSV resources.',
    )
    parser.add_argument('--interaction_file', type = Path, required = True)
    parser.add_argument('--all_output_file', type = Path, required = True)
    parser.add_argument('--hc_output_file', type = Path, required = True)

    args = parser.parse_args()

    all_pairs, high_confidence_pairs = load_domine_pairs(args.interaction_file)
    write_pairs(all_pairs, args.all_output_file)
    write_pairs(high_confidence_pairs, args.hc_output_file)

    print(
        'Built '
        f'{len(all_pairs)} DOMINE pairs and '
        f'{len(high_confidence_pairs)} high-confidence DOMINE pairs.',
    )


if __name__ == '__main__':
    main()
