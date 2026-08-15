"""Console-script entry point for functional enrichment analysis. Argparse boilerplate only."""

import argparse

import pandas as pd


def _add_enrichment_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the target table and the background gene universe."""

    parser.add_argument(
        "--analysis_level",
        required=True,
        choices=["HMI", "TieDIE"],
        help="HMI: enrich the distinct human targets of a DDI/DMI table. TieDIE: enrich every node of a Module 9 network.",
    )
    parser.add_argument(
        "--target_file",
        required=True,
        help="A Module 5/6/7/8 DDI/DMI table (HMI) or a Module 9 final network TSV (TieDIE).",
    )
    parser.add_argument(
        "--background_gene_list",
        required=True,
        help="Expressed-gene-symbol background universe (first column read).",
    )
    parser.add_argument(
        "--sep",
        required=True,
        help="Field separator for --background_gene_list.",
    )


def _build_enrichment_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the enrichment CLI command."""

    parser = argparse.ArgumentParser(
        description="Run functional over-representation (Enrichr) enrichment on human targets.",
    )
    _add_enrichment_source_arguments(parser)
    parser.add_argument(
        "--database",
        default="Reactome_2022",
        help="""Enrichr gene-set library. Default Reactome_2022. Supported gget shortcuts:
                'pathway' (KEGG_2021_Human), 'transcription' (ChEA_2016),
                'ontology' (GO_Biological_Process_2021), 'diseases_drugs' (GWAS_Catalog_2019),
                'celltypes' (PanglaoDB_Augmented_2021), 'kinase_interactions' (KEA_2015),
                or any library at https://maayanlab.cloud/Enrichr/#libraries""",
    )
    parser.add_argument(
        "--ranking",
        default="combined_score",
        choices=["combined_score", "adj_p"],
        help="Rank the plot by 'combined_score' (twin-axis) or 'adj_p'.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for the significant-enrichment results.",
    )
    parser.add_argument(
        "--output_image",
        required=True,
        help="Output PNG file for the top-20 enrichment plot.",
    )
    return parser


def main() -> int:
    """Run functional enrichment on the human targets of an HMI table or a TieDie network."""

    from ..utils import uniprot_client
    from ..workflow import enrichment as enrichment_workflow

    args = _build_enrichment_parser().parse_args()
    # HMI tables are comma-separated (Module 6/7/8); a Module 9 network is tab-separated.
    target_sep = "," if args.analysis_level == "HMI" else "\t"
    target_table = pd.read_csv(args.target_file, sep=target_sep)
    background_symbols = uniprot_client.read_ids(
        args.background_gene_list, args.sep, 1, has_header=False
    )

    results, figure = enrichment_workflow.run_enrichment_analysis(
        target_table,
        background_symbols,
        analysis_level=args.analysis_level,
        database=args.database,
        ranking=args.ranking,
    )

    results.to_csv(args.output_file, index=False)
    if figure is not None:
        figure.savefig(
            args.output_image, dpi=300, bbox_inches="tight", transparent=True
        )
    else:
        print(
            f"No terms passed the adj_p_val < {enrichment_workflow.SIGNIFICANCE_CUTOFF} "
            f"cutoff; wrote empty results to {args.output_file} and skipped the plot."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
