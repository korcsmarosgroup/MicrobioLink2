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
        #print(Uniprot_ids)
        if index == len(Uniprot_ids) - 1:
            Uniprot_ids[index] = '%28accession%3A' + Uniprot_ids[index] + '%29%29'
        else:
            Uniprot_ids[index] = '%28accession%3A' + Uniprot_ids[index] + '%29+OR+'
    print(Uniprot_ids)
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
    args = parse_args(argv)
    ids = read_ids(args.id_list, args.sep, args.id_column)

    BATCH = 1000
    header_written = False

    with open(args.output, "w") as fout:
        for start in range(0, len(ids), BATCH):
            batch = ids[start: start+BATCH]

            if args.id_type == "UP":
                for up_id in batch:
                    tsv = download_proteome(up_id)
                    lines = tsv.rstrip().split("\n")

                    # write header (once) + extra column name
                    if not header_written:
                        fout.write(lines[0] + "\tProteome_ID\n")
                        header_written = True

                    for row in lines[1:]:
                        if row:
                            fout.write(f"{row}\t{up_id}\n")

            else:   # id_type == "Uniprot"
                tsv = download_protein_list(batch)
                lines = tsv.rstrip().split("\n")

                if not header_written:
                    fout.write(lines[0] + "\n")
                    header_written = True

                for row in lines[1:]:
                    if row:
                        fout.write(row + "\n")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
