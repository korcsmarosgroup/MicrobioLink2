import argparse
import requests
import sys

# Creating a list of proteome ids to download from Uniprot
def read_ids(file, separator, id_column):
    with open(file, encoding='utf-8-sig') as id_list:
        id_list.readline()
        ids = []
        for line in id_list:
            line = line.strip().split(separator)

            #Python starts to count by 0, therefore the user-provided column number should be decreased by one
            id_column_new = int(id_column) - 1
            ids.append(line[id_column_new])

    return ids

# OPTION A: Download the whole proteome
def download_proteome(UP_id):
    url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_pfam%2Cgene_names&format=tsv&query=%28%28proteome%3A' + UP_id + '%29%29'
    r = requests.get(url)
    return r.text

# OPTION B: Download the list of proteins with Uniprot ID
def download_protein_list(Uniprot_ids):
    for index, value in enumerate(Uniprot_ids):
        print(Uniprot_ids)
        if index == len(Uniprot_ids) - 1:
            Uniprot_ids[index] = '%28accession%3A' + Uniprot_ids[index] + '%29%29'
        else:
            Uniprot_ids[index] = '%28accession%3A' + Uniprot_ids[index] + '%29+OR+'

    url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_pfam%2Cgene_names&format=tsv&query=%28' + "".join(Uniprot_ids)

    response = requests.get(url)
    if response.status_code == 200:
        return response.text

def parse_args(argv):
    """ Command line interface for the the module """
    parser = argparse.ArgumentParser()
    parser.add_argument("--id_list",
                        help="<Path to an existing FILE describing Uniprot or UP IDs>",
                        dest="id_list",
                        action="store",
                        required=True)

    parser.add_argument("--sep",
                        help="<Field separator>",
                        dest="sep",
                        action="store",
                        required=True)

    parser.add_argument("--id_type",
                        help="<Type of ID - Uniprot or UniprotProteome (UP)>",
                        choices=['Uniprot','UP'],
                        dest="id_type",
                        action="store",
                        required=True)


    parser.add_argument("--id_column",
                        help="<Column number for proteome/protein ID>",
                        dest="id_column",
                        action="store",
                        type=int,
                        required=True)

    parser.add_argument("--output",
                        help="<Path to the output file>",
                        dest="output",
                        action="store",
                        required=True)


    results = parser.parse_args(argv)
    return results


def main(argv):
    """ Main method and logic """

    # Read args
    args = parse_args(argv)

    # Download proteomes or list of proteins with their domains
    ids = read_ids(args.id_list, args.sep, args.id_column)

    batch_size = 1000

    with open(args.output, 'w') as output:
        header_added = False
        for i in range(0, len(ids), batch_size):
            # Split the ids list into batches
            #0:1000; 1000:2000; 2000:4000
            batch_ids = ids[i:i + batch_size]

            # Download protein details for the batch
            batch_results = []
            if args.id_type == 'UP':
                for id_ in batch_ids:
                    result = download_proteome(id_)
                    batch_results.append([result,id_])

            elif args.id_type == 'Uniprot':
                batch_results = []
                result = download_protein_list(batch_ids)
                batch_results.append(result)

            # Add the header only if it hasn't been added yet
            if header_added == False:
                header_added = True
                header = batch_results[0][0].split('\n', 1)[0]  # Extract the header from the first batch result
                output.write(header + '\n')

            # Write the batch results to the main output file, skipping the header in subsequent batches
            for result_list in batch_results:
                result = result_list[0].split('\n')
                print(result)
                for ids in result[1:]:
                    if ids:
                        output.write(ids + "\t" + result_list[1] + "\n")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
