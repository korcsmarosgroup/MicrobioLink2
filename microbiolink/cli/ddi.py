"""Console-script entry point for domain-domain interaction prediction. Argparse boilerplate only."""

import argparse


from ._common import _read_domain_mapping


def _build_ddi_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the ddi CLI command."""

    parser = argparse.ArgumentParser(
        description="Predict domain-domain interactions between bacterial and human proteins.",
    )
    parser.add_argument(
        "-b",
        "--bacterial_domain_file",
        required=True,
        help="Bacterial Pfam/Entries domain TSV file (microbiolink-download-domains output).",
    )
    parser.add_argument(
        "-hu",
        "--human_domain_file",
        required=True,
        help="Human Pfam/Entries domain TSV file (microbiolink-download-domains output).",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for predicted domain-domain interactions.",
    )
    return parser


def main() -> int:
    """Predict domain-domain interactions between bacterial and human proteins."""

    from ..workflow import ddi as ddi_workflow

    args = _build_ddi_parser().parse_args()
    bacterial_domains = _read_domain_mapping(args.bacterial_domain_file)
    human_domains = _read_domain_mapping(args.human_domain_file)

    result = ddi_workflow.predict_domain_domain_interactions(
        bacterial_domains, human_domains
    )
    result.to_csv(args.output_file, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
