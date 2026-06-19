#!/usr/bin/env python

"""Build packaged 3did DMI resources from current 3did flat files."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
import re


PFAM_ACCESSION_PATTERN = re.compile(r'PF\d{5}')


def _build_pfam_name_mapping(flat_path: Path) -> dict[str, str]:
    """Map 3did Pfam names to Pfam accessions using 3did_flat.gz."""

    mapping: dict[str, str] = {}

    with gzip.open(flat_path, 'rt', encoding = 'utf-8') as flat_file:
        for line in flat_file:
            if not line.startswith('#=ID'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            domain_names = fields[1:3]
            accessions = PFAM_ACCESSION_PATTERN.findall(line)

            if len(accessions) != 2:
                continue

            for domain_name, accession in zip(domain_names, accessions):
                existing = mapping.get(domain_name)

                if existing is None:
                    mapping[domain_name] = accession
                    continue

                if existing != accession:
                    raise ValueError(
                        f'Pfam name {domain_name} mapped to multiple '
                        f'accessions: {existing}, {accession}',
                    )

    return mapping


def build_resources(
    ddi_flat_path: Path,
    dmi_flat_path: Path,
    regex_output_path: Path,
    motif_domain_output_path: Path,
) -> tuple[int, list[str]]:
    """Write packaged 3did regex and motif-domain tables."""

    pfam_name_mapping = _build_pfam_name_mapping(ddi_flat_path)
    regex_rows: dict[str, dict[str, str]] = {}
    motif_domain_rows: dict[tuple[str, str], dict[str, str]] = {}
    unresolved_domains: list[str] = []
    current_domain_name: str | None = None
    current_motif_name: str | None = None

    with gzip.open(dmi_flat_path, 'rt', encoding = 'utf-8') as dmi_file:
        for line in dmi_file:
            stripped_line = line.strip()

            if stripped_line.startswith('#=ID'):
                fields = stripped_line.split('\t')
                if len(fields) < 4:
                    current_domain_name = None
                    current_motif_name = None
                    continue

                current_domain_name = fields[1]
                current_motif_name = f'3DID_{fields[3]}'
                continue

            if not stripped_line.startswith('#=PT'):
                continue

            if current_domain_name is None or current_motif_name is None:
                continue

            pfam_accession = pfam_name_mapping.get(current_domain_name)
            if pfam_accession is None:
                unresolved_domains.append(current_domain_name)
                continue

            motif_pattern = stripped_line.split('\t', 1)[1].rsplit(' (', 1)[0]

            existing_regex_row = regex_rows.get(current_motif_name)
            if existing_regex_row is not None:
                if existing_regex_row['Regex'] != motif_pattern:
                    raise ValueError(
                        f'Motif {current_motif_name} has multiple patterns: '
                        f"{existing_regex_row['Regex']} vs {motif_pattern}",
                    )
            else:
                regex_rows[current_motif_name] = {
                    'Accession': current_motif_name,
                    'Identifier': current_motif_name,
                    'Description': f'3did structural motif for {current_domain_name}',
                    'Type': 'pattern',
                    'Regex': motif_pattern,
                }

            motif_domain_rows[(current_motif_name, pfam_accession)] = {
                'Motif': current_motif_name,
                'Domain': pfam_accession,
                'Interaction Domain Description': current_domain_name,
                'Interaction Domain Name': current_domain_name,
            }

    with open(regex_output_path, 'w', encoding = 'utf-8', newline = '') as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames = [
                'Accession',
                'Identifier',
                'Description',
                'Type',
                'Regex',
            ],
            delimiter = '\t',
        )
        writer.writeheader()
        writer.writerows(
            regex_rows[motif_name]
            for motif_name in sorted(regex_rows)
        )

    with open(
        motif_domain_output_path,
        'w',
        encoding = 'utf-8',
        newline = '',
    ) as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames = [
                'Motif',
                'Domain',
                'Interaction Domain Description',
                'Interaction Domain Name',
            ],
            delimiter = '\t',
        )
        writer.writeheader()
        writer.writerows(
            motif_domain_rows[key]
            for key in sorted(motif_domain_rows)
        )

    return len(regex_rows), sorted(set(unresolved_domains))


def main() -> None:
    """Build packaged 3did DMI resource tables."""

    parser = argparse.ArgumentParser(
        description = 'Build packaged 3did DMI TSV resources.',
    )
    parser.add_argument('--ddi_flat_path', type = Path, required = True)
    parser.add_argument('--dmi_flat_path', type = Path, required = True)
    parser.add_argument('--regex_output_path', type = Path, required = True)
    parser.add_argument(
        '--motif_domain_output_path',
        type = Path,
        required = True,
    )

    args = parser.parse_args()

    rule_count, unresolved_domains = build_resources(
        ddi_flat_path = args.ddi_flat_path,
        dmi_flat_path = args.dmi_flat_path,
        regex_output_path = args.regex_output_path,
        motif_domain_output_path = args.motif_domain_output_path,
    )

    print(
        f'Built {rule_count} 3did DMI rules. '
        f'Unresolved domains: {", ".join(unresolved_domains) or "none"}.',
    )


if __name__ == '__main__':
    main()
