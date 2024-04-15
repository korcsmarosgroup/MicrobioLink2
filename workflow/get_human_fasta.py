import omnipath as op
import requests
import argparse
from mygene import MyGeneInfo
import gget

# Function to retrieve proteins (optionally based on there location) from gene list
def get_proteins(gene_expression_file, id_type, sep, pmtm, location_filter_list):
    proteins = []

    if location_filter_list:
        # Fetch intercell data with optional parent parameter
        pmtm = op.requests.Intercell.get(
            parent=location_filter_list,
            scope=['generic', 'specific'],
            source=['resource_specific', 'composite'],
            entity_type='protein',
        )

        with open(gene_expression_file) as gene_expression
            with open('location_filtered_genes.csv', 'w') as output:
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) >1:
                        gene = line[0]
                        expression = line[1]

                        if expression > 0.0 or expression < 0.0:
                            if id_type == 'genesymbol':
                               if gene in list(pmtm['genesymbol']):
                                    proteins.append(gene)
                            elif id_type == 'uniprot':
                                 if gene in list(pmtm['uniprot']):
                                     # Find the corresponding 'genesymbol' and add it to the list
                                     genesymbol = pmtm.loc[pmtm['uniprot'] == gene, 'genesymbol'].values
                                     # genesymbol = array(['CFC1', 'CFC1', 'CFC1', 'CFC1', 'CFC1'], dtype=object)
                                     if len(genesymbol) > 0:
                                         proteins.append(genesymbol[0])
                            output.write(",".join(line) + "\n")


    else:
        # Directly download all protein data without filtering
        if id_type == 'uniprot':
            uniprot_ids = []
            with open(gene_expression_file) as gene_expression:
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) >1:
                        gene = line[0]
                        expression = line[1]
                        if expression > 0.0 or expression < 0.0:
                            uniprot_ids.append(gene)

            translation_dict = translate_uniprot_to_symbols(uniprot_ids)
            proteins.extend(translation_dict.values())

        elif id_type == 'genesymbol':
            with open(gene_expression_file) as gene_expression:
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) >1:
                        gene = line[0]
                        expression = line[1]
                        if expression > 0.0 or expression < 0.0:
                            proteins.append(gene)

    return proteins

#Translate UniProt IDs to gene symbols.
def translate_uniprot_to_symbols(uniprot, species='human'):
    mg = MyGeneInfo()
    target_genesymbols_translation = mg.querymany(uniprot, scopes='uniprot', fields='symbol', species='human', returnall=True)

    translation_dict = {entry['query']: entry['symbol'] for entry in target_genesymbols_translation['out'] if 'symbol' in entry}
    return translation_dict

# Function to fetch protein sequences
def fetch_protein_sequences(gene_symbols):
    fasta_data = []

    for index, value in enumerate(gene_symbols):
        if index == len(gene_symbols) - 1:
            gene_symbols[index] = '%28gene%3A' + gene_symbols[index] + '%29%29'
        else:
            gene_symbols[index] = '%28gene%3A' + gene_symbols[index] + '%29+OR+'

    url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28' + "".join(gene_symbols) + '+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29'
    response = requests.get(url)
    print(url)

    if response.status_code == 200:
        fasta_data.append(response.text)
    else:
        print(f"Failed to fetch data for gene symbols: {gene_symbols}")

    return fasta_data

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Retrieve protein sequences for proteins from a gene list.')

    # Define command-line arguments
    parser.add_argument('-genes','--gene_expression', type=str, help='Path to the transcripomics data')
    parser.add_argument('-id','--id_type', choices=['genesymbol', 'uniprot'], help='Type of gene identifier (genesymbol or uniprot)')
    parser.add_argument('-s','--sep', help='Field separator in the protein list file')
    parser.add_argument('-lfl', '--location_filter_list', type=list, default=None, help='Location filter list (options:plasma_membrane_transmembrane and/or plasma_membrane_peripheral and/or secreted), (required format: []), (default: None)')
    parser.add_argument('-o', '--output', default='protein_sequences.fasta', help='Output file for protein sequences')


    # Parse the command-line arguments
    args = parser.parse_args()

    # Get proteins from the gene list
    proteins = get_proteins(args.gene_expression, args.id_type, args.sep, pmtm, args.location_filter_list)
    print(proteins)

    batch_size = 100
    fasta_sequences = []
    for i in range(0, len(proteins), batch_size):
        # Split the ids list into batches
        #0:100; 100:200; 200:300
        batch_ids = proteins[i:i + batch_size]

    # Fetch protein sequences for proteins
        fasta_sequences.extend(fetch_protein_sequences(batch_ids))

    # Save the protein sequences to the specified output file
    with open(args.output, "w") as fasta_file:
        fasta_file.write("".join(fasta_sequences))

    print(f"Protein sequences saved to {args.output}")

if __name__ == '__main__':
    main()
