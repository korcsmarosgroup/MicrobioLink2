#!/usr/bin/env python

"""Human gene-identifier resolution workflow for domain/FASTA lookups."""

from __future__ import annotations

import sys
from pathlib import Path

import requests

from microbiolink.core.uniprot import fetch_fasta_sequences


PathLike = str | Path


def _is_expressed(value: str) -> bool:
    """Return True if an expression value is non-NaN and non-zero."""

    if value in ('NaN', '', 'nan'):
        return False
    try:
        return float(value) != 0.0
    except ValueError:
        return False


def read_expressed_genes(
    filename: PathLike,
    separator: str,
) -> list[str]:
    """Read expressed gene or protein IDs from a z-score filtered CSV.

    Includes any row where at least one expression column is non-NaN and
    non-zero, consistent with z_score_filter_terminal.py output format.

    Args:
        filename: Path to the z-score filtered CSV file.
        separator: Field separator used in the file.

    Returns:
        List of gene or protein identifiers for expressed entries.
    """

    genes: list[str] = []

    with open(Path(filename), encoding = 'utf-8-sig') as gene_file:
        next(gene_file, None)

        for line in gene_file:
            fields = line.strip().split(separator)
            if len(fields) < 2:
                continue
            gene = fields[0]
            if any(_is_expressed(v) for v in fields[1:]):
                genes.append(gene)

    return genes


def _extract_swissprot_ids(translation_dict: dict) -> list[str]:
    """Extract Swiss-Prot UniProt accessions from a MyGene translation dict.

    Args:
        translation_dict: Mapping of gene symbol to MyGene uniprot entry.

    Returns:
        Flat list of Swiss-Prot UniProt accession strings.
    """

    proteins: list[str] = []
    for entry in translation_dict.values():
        if 'Swiss-Prot' not in entry:
            continue
        value = entry['Swiss-Prot']
        if isinstance(value, str):
            proteins.append(value)
        elif isinstance(value, list):
            proteins.extend(value)
    return proteins


def translate_symbol_to_uniprot(symbol, species = 'human'):
    """Translate gene symbols to UniProt accessions via MyGene.

    Args:
        symbol: Gene symbol or list of gene symbols.
        species: Species name (default: 'human').

    Returns:
        Dict mapping query symbol to UniProt entry dict.
    """

    from mygene import MyGeneInfo

    mg = MyGeneInfo()
    target_genesymbols_translation = mg.querymany(
        symbol,
        scopes = 'symbol',
        fields = 'uniprot',
        species = 'human',
        returnall = True,
    )
    translation_dict = {
        entry['query']: entry['uniprot']
        for entry in target_genesymbols_translation['out']
        if 'uniprot' in entry
    }
    return translation_dict


def get_proteins(
    gene_expression_file,
    id_type,
    sep,
    location_filter_list,
    output_folder,
):
    """Resolve which proteins to fetch from a gene-expression file.

    Filters a raw (not-yet-z-score-filtered) gene-expression file down to
    expressed genes, optionally restricting to a subcellular location via
    the omnipath intercell table, and resolves the result to UniProt
    accessions.

    Args:
        gene_expression_file: Path to the raw gene-expression file.
        id_type: Identifier type: 'genesymbol' or 'uniprot'.
        sep: Field separator in the gene-expression file.
        location_filter_list: Optional list of subcellular location filters
            (omnipath intercell parent categories). When given, genes are
            additionally required to appear in the intercell table for one
            of these locations, and matching rows are written to
            `location_filtered_genes.csv` in `output_folder`.
        output_folder: Folder to write `location_filtered_genes.csv` into
            when `location_filter_list` is given.

    Returns:
        Deduplicated list of UniProt accessions.
    """

    proteins = []

    if location_filter_list:
        import omnipath as op

        pmtm = op.requests.Intercell.get(
            parent = location_filter_list,
            scope = ['generic', 'specific'],
            source = ['resource_specific', 'composite'],
            entity_type = 'protein',
        )

        with open(gene_expression_file) as gene_expression:
            gene_expression.readline()
            location_output_path = Path(output_folder) / 'location_filtered_genes.csv'
            with open(location_output_path, 'w') as location_output:
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) > 1 and _is_expressed(line[1]):
                        gene = line[0]
                        if id_type == 'genesymbol':
                            if gene in list(pmtm['genesymbol']):
                                uniprot = pmtm.loc[pmtm['genesymbol'] == gene, 'uniprot'].values
                                if len(uniprot) > 0:
                                    proteins.append(uniprot[0])
                                    location_output.write(','.join(line) + '\n')

                        elif id_type == 'uniprot':
                            if gene in list(pmtm['uniprot']):
                                proteins.append(gene)
                                location_output.write(','.join(line) + '\n')

    else:
        if id_type == 'uniprot':
            proteins = read_expressed_genes(gene_expression_file, sep)

        elif id_type == 'genesymbol':
            symbols = read_expressed_genes(gene_expression_file, sep)
            translation_dict = translate_symbol_to_uniprot(symbols)
            proteins = _extract_swissprot_ids(translation_dict)

    proteins = list(set(proteins))
    return proteins


def fetch_protein_sequences(uniprots):
    """Fetch FASTA sequences for a list of UniProt accessions.

    Delegates to fetch_fasta_sequences and preserves skip-and-continue
    behaviour on HTTP errors for backward compatibility.

    Args:
        uniprots: List of UniProt accession strings.

    Returns:
        List containing a single FASTA string, or empty list on failure.
    """

    try:
        return [fetch_fasta_sequences(uniprots)]
    except requests.HTTPError:
        print(
            f'Warning: failed to fetch sequences for batch of {len(uniprots)} accessions.',
            file = sys.stderr,
        )
        return []
