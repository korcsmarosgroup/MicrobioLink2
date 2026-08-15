"""Console-script entry point for downloading Pfam domain annotations. Argparse boilerplate only."""

import argparse
from pathlib import Path

from ._common import (
    _add_human_identifier_arguments,
    _add_microbial_identifier_arguments,
    _resolve_human_identifiers,
    _resolve_microbial_identifiers,
)


def _build_domain_download_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the domain-download CLI command."""

    parser = argparse.ArgumentParser(
        description="Download Pfam domain annotations for human and/or microbial proteins.",
    )
    _add_human_identifier_arguments(parser)
    _add_microbial_identifier_arguments(parser)
    parser.add_argument(
        "-o",
        "--output_folder",
        required=True,
        help="Folder to write human_domains.tsv / microbial_domains.tsv into.",
    )
    return parser


def _write_domain_table(domains: dict[str, list[str]], output_path: Path) -> None:
    """Write a pfam_id -> uniprot_ids mapping as a Pfam/Entries TSV."""

    with open(output_path, "w") as output_file:
        output_file.write("Pfam\tEntries\n")
        for pfam_id, uniprot_ids in domains.items():
            entries_field = ";".join(uniprot_ids) + ";" if uniprot_ids else ""
            output_file.write(f"{pfam_id}\t{entries_field}\n")


def main() -> int:
    """Download Pfam domains for human and/or microbial proteins."""

    from ..workflow import domain_download

    args = _build_domain_download_parser().parse_args()
    human_identifiers = _resolve_human_identifiers(args)
    microbial_identifiers = _resolve_microbial_identifiers(args)

    results = domain_download.download_domains(
        human_identifiers=human_identifiers,
        human_id_type=args.human_id_type,
        microbial_identifiers=microbial_identifiers,
        microbial_id_type=args.microbial_id_type,
    )

    output_folder = Path(args.output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    filenames = {"human": "human_domains.tsv", "microbial": "microbial_domains.tsv"}
    for species, domains in results.items():
        _write_domain_table(domains, output_folder / filenames[species])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
