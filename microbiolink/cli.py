"""Console-script entry points for MicrobioLink. Argparse boilerplate only."""

import argparse
from pathlib import Path

import pandas as pd


def zscore_filter() -> int:
    """Filter a gene/protein count matrix by a z-score cutoff."""

    from .workflow import zscore_filter as zscore_filter_module

    parser = argparse.ArgumentParser(
        description="Gene expression filtration based on individual cell count and z-score.",
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Input CSV file with gene expression data.",
    )
    parser.add_argument(
        "-zscore",
        "--zscore",
        required=True,
        type=float,
        help="Z-score cut-off to filter lowly expressed genes.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for filtered results.",
    )
    args = parser.parse_args()

    zscore_filter_module.filter_count_matrix_file(
        args.input_file,
        zscore_threshold=args.zscore,
        output_file=args.output_file,
    )
    return 0


def _add_identifier_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing where identifiers come from."""

    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Identifier list file, or a count matrix if --from-count-matrix is set.",
    )
    parser.add_argument(
        "--from-count-matrix",
        action="store_true",
        help="Treat --input_file as a gene count matrix and extract gene symbols from it.",
    )
    parser.add_argument(
        "-sep",
        "--sep",
        default=",",
        help="Field separator for a plain identifier list (ignored with --from-count-matrix).",
    )
    parser.add_argument(
        "-col",
        "--id_column",
        type=int,
        default=1,
        help="One-based identifier column for a plain list (ignored with --from-count-matrix).",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Treat the identifier list as having no header row (ignored with --from-count-matrix).",
    )


def _add_membrane_filter_target_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the filter and output."""

    parser.add_argument(
        "-id",
        "--id_type",
        required=True,
        choices=["uniprot", "genesymbol", "proteome"],
        help="Identifier type in the input file.",
    )
    parser.add_argument(
        "-sp",
        "--species",
        required=True,
        choices=["human", "microbial"],
        help="Species the identifiers belong to.",
    )
    parser.add_argument(
        "-lfl",
        "--location_filters",
        required=True,
        nargs="+",
        help="Location categories (human) or location substrings (microbial) to keep.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for the filtered results.",
    )


def _build_membrane_filter_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the membrane-filter CLI command."""

    parser = argparse.ArgumentParser(
        description="Filter human or microbial proteins down to membrane/secreted proteins.",
    )
    _add_identifier_source_arguments(parser)
    _add_membrane_filter_target_arguments(parser)
    return parser


def _resolve_membrane_filter_identifiers(args: argparse.Namespace) -> list[str]:
    """Resolve the identifier list from CLI args, per --from-count-matrix."""

    from .utils import uniprot_client

    if args.from_count_matrix:
        count_matrix = pd.read_csv(args.input_file, index_col=0)
        return count_matrix.index.tolist()

    return uniprot_client.read_ids(
        args.input_file, args.sep, args.id_column, has_header=not args.no_header
    )


def membrane_filter() -> int:
    """Filter human or microbial proteins to membrane/secreted proteins."""

    from .workflow import membrane_filter as membrane_filter_module

    args = _build_membrane_filter_parser().parse_args()
    identifiers = _resolve_membrane_filter_identifiers(args)

    result = membrane_filter_module.filter_membrane_proteins(
        identifiers,
        id_type=args.id_type,
        species=args.species,
        location_filters=args.location_filters,
    )
    result.to_csv(args.output_file, index=False)
    return 0


def _add_human_identifier_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments for the human fasta identifiers (all optional)."""

    parser.add_argument(
        "-hi",
        "--human_input_file",
        help="Human identifier list file, or a count matrix if --human-from-count-matrix is set.",
    )
    parser.add_argument(
        "-hid",
        "--human_id_type",
        choices=["uniprot", "genesymbol"],
        help="Identifier type in --human_input_file.",
    )
    parser.add_argument(
        "--human-from-count-matrix",
        action="store_true",
        help="Treat --human_input_file as a gene count matrix and extract expressed gene symbols from it.",
    )
    parser.add_argument(
        "-hsep",
        "--human_sep",
        default=",",
        help="Field separator for a plain human identifier list.",
    )
    parser.add_argument(
        "-hcol",
        "--human_id_column",
        type=int,
        default=1,
        help="One-based identifier column for a plain human identifier list.",
    )
    parser.add_argument(
        "--human-no-header",
        action="store_true",
        help="Treat the human identifier list as having no header row.",
    )


def _add_microbial_identifier_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments for the microbial fasta identifiers (all optional)."""

    parser.add_argument(
        "-mi",
        "--microbial_input_file",
        help="Microbial identifier list file.",
    )
    parser.add_argument(
        "-mid",
        "--microbial_id_type",
        choices=["uniprot", "proteome"],
        help="Identifier type in --microbial_input_file.",
    )
    parser.add_argument(
        "-msep",
        "--microbial_sep",
        default=",",
        help="Field separator for the microbial identifier list.",
    )
    parser.add_argument(
        "-mcol",
        "--microbial_id_column",
        type=int,
        default=1,
        help="One-based identifier column for the microbial identifier list.",
    )
    parser.add_argument(
        "--microbial-no-header",
        action="store_true",
        help="Treat the microbial identifier list as having no header row.",
    )


def _build_fasta_download_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the fasta-download CLI command."""

    parser = argparse.ArgumentParser(
        description="Download human and/or microbial protein sequences as FASTA files.",
    )
    _add_human_identifier_arguments(parser)
    _add_microbial_identifier_arguments(parser)
    parser.add_argument(
        "-o",
        "--output_folder",
        required=True,
        help="Folder to write human_proteins.fasta / microbial_proteins.fasta into.",
    )
    return parser


def _resolve_human_identifiers(args: argparse.Namespace) -> list[str] | None:
    """Resolve the human identifier list from CLI args, or None if not supplied."""

    if args.human_input_file is None:
        return None

    from .utils import uniprot_client

    if args.human_from_count_matrix:
        count_matrix = pd.read_csv(args.human_input_file, index_col=0)
        expressed_mask = count_matrix.notna().any(axis=1) & (count_matrix != 0).any(
            axis=1
        )
        return count_matrix.index[expressed_mask].tolist()

    return uniprot_client.read_ids(
        args.human_input_file,
        args.human_sep,
        args.human_id_column,
        has_header=not args.human_no_header,
    )


def _resolve_microbial_identifiers(args: argparse.Namespace) -> list[str] | None:
    """Resolve the microbial identifier list from CLI args, or None if not supplied."""

    if args.microbial_input_file is None:
        return None

    from .utils import uniprot_client

    return uniprot_client.read_ids(
        args.microbial_input_file,
        args.microbial_sep,
        args.microbial_id_column,
        has_header=not args.microbial_no_header,
    )


def download_fasta() -> int:
    """Download human and/or microbial protein sequences as FASTA files."""

    from .workflow import fasta_download

    args = _build_fasta_download_parser().parse_args()
    human_identifiers = _resolve_human_identifiers(args)
    microbial_identifiers = _resolve_microbial_identifiers(args)

    fasta_download.download_fasta(
        args.output_folder,
        human_identifiers=human_identifiers,
        human_id_type=args.human_id_type,
        microbial_identifiers=microbial_identifiers,
        microbial_id_type=args.microbial_id_type,
    )
    return 0


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


def download_domains() -> int:
    """Download Pfam domains for human and/or microbial proteins."""

    from .workflow import domain_download

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


def _read_domain_mapping(filename: str) -> dict[str, list[str]]:
    """Read a Pfam/Entries domain TSV into a pfam_id -> uniprot_ids mapping."""

    domain_table = pd.read_csv(filename, sep="\t")
    return {
        row["Pfam"]: [entry for entry in row["Entries"].split(";") if entry]
        for _, row in domain_table.iterrows()
    }


def ddi() -> int:
    """Predict domain-domain interactions between bacterial and human proteins."""

    from .workflow import ddi as ddi_workflow

    args = _build_ddi_parser().parse_args()
    bacterial_domains = _read_domain_mapping(args.bacterial_domain_file)
    human_domains = _read_domain_mapping(args.human_domain_file)

    result = ddi_workflow.predict_domain_domain_interactions(
        bacterial_domains, human_domains
    )
    result.to_csv(args.output_file, index=False)
    return 0


def _build_dmi_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the dmi CLI command."""

    parser = argparse.ArgumentParser(
        description="Predict domain-motif interactions between bacterial and human proteins.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        choices=["forward", "reverse", "both"],
        help="forward: human motif -> bacterial domain. reverse: bacterial motif -> human domain.",
    )
    parser.add_argument(
        "-hf",
        "--human_fasta_file",
        help="Human protein FASTA file (microbiolink-download-fasta output). Required for forward/both.",
    )
    parser.add_argument(
        "-b",
        "--bacterial_domain_file",
        help="Bacterial Pfam/Entries domain TSV file (microbiolink-download-domains output). Required for forward/both.",
    )
    parser.add_argument(
        "-bf",
        "--bacterial_fasta_file",
        help="Bacterial protein FASTA file (microbiolink-download-fasta output). Required for reverse/both.",
    )
    parser.add_argument(
        "-hu",
        "--human_domain_file",
        help="Human Pfam/Entries domain TSV file (microbiolink-download-domains output). Required for reverse/both.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for predicted domain-motif interactions.",
    )
    return parser


def dmi() -> int:
    """Predict domain-motif interactions between bacterial and human proteins."""

    from .utils import fasta
    from .workflow import dmi as dmi_workflow

    args = _build_dmi_parser().parse_args()
    human_sequences = (
        fasta.read_fasta_sequences(args.human_fasta_file)
        if args.human_fasta_file
        else None
    )
    bacterial_domains = (
        _read_domain_mapping(args.bacterial_domain_file)
        if args.bacterial_domain_file
        else None
    )
    bacterial_sequences = (
        fasta.read_fasta_sequences(args.bacterial_fasta_file)
        if args.bacterial_fasta_file
        else None
    )
    human_domains = (
        _read_domain_mapping(args.human_domain_file) if args.human_domain_file else None
    )

    result = dmi_workflow.predict_domain_motif_interactions(
        args.mode,
        human_sequences=human_sequences,
        bacterial_domains=bacterial_domains,
        bacterial_sequences=bacterial_sequences,
        human_domains=human_domains,
    )
    result.to_csv(args.output_file, index=False)
    return 0


def _add_idr_filter_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the DMI table and its source FASTA files."""

    parser.add_argument(
        "--dmi_file",
        required=True,
        help="DMI predictions CSV file (microbiolink-dmi output).",
    )
    parser.add_argument(
        "--human_fasta_file",
        help="Human protein FASTA file (microbiolink-download-fasta output). Required if dmi_file has 'forward' rows.",
    )
    parser.add_argument(
        "--bacterial_fasta_file",
        help="Bacterial protein FASTA file (microbiolink-download-fasta output). Required if dmi_file has 'reverse' rows.",
    )


def _add_idr_filter_scoring_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the disorder/binding method, cutoffs, and device."""

    parser.add_argument(
        "--method",
        required=True,
        choices=["iupred", "aiupred"],
        help="Disorder/binding prediction method.",
    )
    parser.add_argument(
        "--disorder_cutoff",
        required=True,
        type=float,
        help="Minimum per-residue disorder score required across the whole motif window.",
    )
    parser.add_argument(
        "--binding_cutoff",
        required=True,
        type=float,
        help="Minimum per-residue binding score required across the whole motif window.",
    )
    parser.add_argument(
        "--force-cpu",
        dest="force_cpu",
        action="store_true",
        help="Force CPU inference (aiupred only; ignored for iupred).",
    )
    parser.add_argument(
        "--gpu-num",
        dest="gpu_num",
        type=int,
        default=0,
        help="GPU index to use (aiupred only; ignored for iupred).",
    )


def _build_idr_filter_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the idr-filter CLI command."""

    parser = argparse.ArgumentParser(
        description="Filter domain-motif interactions to those in a disordered, binding-prone region.",
    )
    _add_idr_filter_source_arguments(parser)
    _add_idr_filter_scoring_arguments(parser)
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for the filtered domain-motif interactions.",
    )
    return parser


def idr_filter() -> int:
    """Filter domain-motif interactions to those in a disordered, binding-prone region."""

    from .utils import fasta
    from .workflow import idr_filter as idr_filter_module

    args = _build_idr_filter_parser().parse_args()
    dmi_table = pd.read_csv(args.dmi_file)
    human_sequences = (
        fasta.read_fasta_sequences(args.human_fasta_file)
        if args.human_fasta_file
        else None
    )
    bacterial_sequences = (
        fasta.read_fasta_sequences(args.bacterial_fasta_file)
        if args.bacterial_fasta_file
        else None
    )

    result = idr_filter_module.filter_by_disorder(
        dmi_table,
        human_sequences=human_sequences,
        bacterial_sequences=bacterial_sequences,
        method=args.method,
        disorder_cutoff=args.disorder_cutoff,
        binding_cutoff=args.binding_cutoff,
        force_cpu=args.force_cpu,
        gpu_num=args.gpu_num,
    )
    result.to_csv(args.output_file, index=False)
    return 0
