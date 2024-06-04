import omnipath as op
import requests
import argparse
from mygene import MyGeneInfo
import gget
import os

# Function to retrieve proteins (optionally based on there location) from gene list
def get_proteins(gene_expression_file, id_type, sep, location_filter_list, output_folder):
    proteins = []

    if location_filter_list:
        # Fetch intercell data with optional parent parameter
        pmtm = op.requests.Intercell.get(
            parent=location_filter_list,
            scope=['generic', 'specific'],
            source=['resource_specific', 'composite'],
            entity_type='protein',
        )

        with open(gene_expression_file) as gene_expression:
            gene_expression.readline()
            with open(os.path.join(output_folder, 'location_filtered_genes.csv'), 'w') as location_output:
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) >1:
                        gene = line[0]
                        if line[1] !='NaN':
                            expression = float(line[1])
                            if expression != 0.0:
                                if id_type == 'genesymbol':
                                    if gene in list(pmtm['genesymbol']):
                                     # Find the corresponding 'uniprots' and add it to the list
                                     # uniprot = array([''], dtype=object)
                                        uniprot = pmtm.loc[pmtm['genesymbol'] == gene, 'uniprot'].values
                                        if len(uniprot) > 0:
                                            proteins.append(uniprot[0])
                                            location_output.write(",".join(line) + "\n")

                                elif id_type == 'uniprot':
                                    if gene in list(pmtm['uniprot']):
                                        proteins.append(gene)
                                        location_output.write(",".join(line) + "\n")

    else:
        # Directly download all protein data without filtering
        if id_type == 'uniprot':
            uniprot_ids = []
            with open(gene_expression_file) as gene_expression:
                gene_expression.readline()
                for line in gene_expression:
                    line = line.strip().split(sep)
                    if len(line) >1:
                        gene = line[0]
                        if line[1] !='NaN':
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
                    if len(line) >1:
                        gene = line[0]
                        if line[1] !='NaN':
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

    proteins=list(set(proteins))
    return proteins

#Translate UniProt IDs to gene symbols.
def translate_symbol_to_uniprot(symbol, species='human'):
    mg = MyGeneInfo()
    target_genesymbols_translation = mg.querymany(symbol, scopes='symbol', fields='uniprot', species='human', returnall=True)

    translation_dict = {entry['query']: entry['uniprot'] for entry in target_genesymbols_translation['out'] if 'uniprot' in entry}
    return translation_dict

# Function to fetch protein sequences
def fetch_protein_sequences(uniprots):
    fasta_data = []

    
    url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=accession:(' + "+OR+".join(uniprots) + ')'
    response = requests.get(url)
    print(url)

    if response.status_code == 200:
        fasta_data.append(response.text)
    else:
       print(f"Failed to fetch data for uniprots: {uniprots}")

    return fasta_data

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Retrieve protein sequences for proteins from a gene list.')

    # Define command-line arguments
    parser.add_argument('-genes','--gene_expression', type=str, help='Path to the transcripomics data')
    parser.add_argument('-id','--id_type', choices=['genesymbol', 'uniprot'], help='Type of gene identifier (genesymbol or uniprot)')
    parser.add_argument('-s','--sep', help='Field separator in the protein list file')
    parser.add_argument('-lfl', '--location_filter_list', nargs='+', default=None, help='Location filter list (options:plasma_membrane_transmembrane and/or plasma_membrane_peripheral and/or secreted), (required format: without '', separated by spaces), (default: None)')
    parser.add_argument('-of', '--output_folder', default='.', help='Output folder for result files')
    parser.add_argument('-oseq', '--output_sequences', default='protein_sequences.fasta', help='Output file for protein sequences')


    # Parse the command-line arguments
    args = parser.parse_args()

    # Get proteins from the gene list
    proteins = get_proteins(args.gene_expression, args.id_type, args.sep, args.location_filter_list, args.output_folder)
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
    with open(os.path.join(args.output_folder, args.output_sequences), 'w') as fasta_file:
        fasta_file.write("".join(fasta_sequences))

    print(f"Protein sequences saved to {args.output_sequences}")

if __name__ == '__main__':
    main()
