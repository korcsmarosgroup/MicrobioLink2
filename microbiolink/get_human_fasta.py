import argparse
import sys
from pathlib import Path

import requests

from microbiolink.get_protein_fasta import fetch_fasta_sequences


# Function to retrieve proteins (optionally based on their location) from gene list
def get_proteins(gene_expression_file, id_type, sep, location_filter_list, output_folder):
    proteins = []

    if location_filter_list:
        import omnipath as op

        # Fetch intercell data with optional parent parameter
        pmtm = op.requests.Intercell.get(
            parent=location_filter_list,
            scope=['generic', 'specific'],
            source=['resource_specific', 'composite'],
            entity_type='protein',
        )

        with open(gene_expression_file) as gene_expression:
            gene_expression.readline()
            location_output_path = Path(output_folder) / 'location_filtered_genes.csv'
            with open(location_output_path, 'w') as location_output:
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) > 1:
                        gene = line[0]
                        if line[1] != 'NaN':
                            expression = float(line[1])
                            if expression != 0.0:
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
        # Directly download all protein data without filtering
        if id_type == 'uniprot':
            with open(gene_expression_file) as gene_expression:
                gene_expression.readline()
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) > 1:
                        gene = line[0]
                        if line[1] != 'NaN':
                            expression = float(line[1])
                            if expression != 0.0:
                                proteins.append(gene)

        elif id_type == 'genesymbol':
            uniprots = []
            with open(gene_expression_file) as gene_expression:
                gene_expression.readline()
                symbols = []
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) > 1:
                        gene = line[0]
                        if line[1] != 'NaN':
                            expression = float(line[1])
                            if expression != 0.0:
                                symbols.append(gene)
            translation_dict = translate_symbol_to_uniprot(symbols)
            proteins.extend(translation_dict.values())

            uniprots.extend(translation_dict.values())

            proteins = []
            for protein in uniprots:
                for ids in protein:
                    if ids == 'Swiss-Prot':
                        if isinstance(protein[ids], str):
                            proteins.append(protein[ids])
                        elif isinstance(protein[ids], list):
                            for i in protein[ids]:
                                proteins.append(i)

    proteins = list(set(proteins))
    return proteins


def translate_symbol_to_uniprot(symbol, species='human'):
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
        scopes='symbol',
        fields='uniprot',
        species='human',
        returnall=True,
    )
    translation_dict = {
        entry['query']: entry['uniprot']
        for entry in target_genesymbols_translation['out']
        if 'uniprot' in entry
    }
    return translation_dict


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
            file=sys.stderr,
        )
        return []


def main():
    parser = argparse.ArgumentParser(
        description='Retrieve protein sequences for proteins from a gene list.',
    )
    parser.add_argument('-genes', '--gene_expression', type=str, help='Path to the transcriptomics data')
    parser.add_argument('-id', '--id_type', choices=['genesymbol', 'uniprot'], help='Type of gene identifier (genesymbol or uniprot)')
    parser.add_argument('-s', '--sep', help='Field separator in the protein list file')
    parser.add_argument('-lfl', '--location_filter_list', nargs='+', default=None, help='Location filter list (options: plasma_membrane_transmembrane and/or plasma_membrane_peripheral and/or secreted)')
    parser.add_argument('-of', '--output_folder', default='.', help='Output folder for result files')
    parser.add_argument('-oseq', '--output_sequences', default='protein_sequences.fasta', help='Output file for protein sequences')

    args = parser.parse_args()

    proteins = get_proteins(
        args.gene_expression,
        args.id_type,
        args.sep,
        args.location_filter_list,
        args.output_folder,
    )
    print(proteins)

    batch_size = 100
    fasta_sequences = []

    for i in range(0, len(proteins), batch_size):
        batch_ids = proteins[i : i + batch_size]
        fasta_sequences.extend(fetch_protein_sequences(batch_ids))

    output_path = Path(args.output_folder) / args.output_sequences
    with open(output_path, 'w') as fasta_file:
        fasta_file.write(''.join(fasta_sequences))

    print(f'Protein sequences saved to {args.output_sequences}')


if __name__ == '__main__':
    main()
