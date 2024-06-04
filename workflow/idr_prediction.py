import argparse
import os
import sys
from iupred2a import read_seq, iupred, anchor2
from Bio import SeqIO



def parse_args(argv):
    """ Command line interface for the the module """
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmi_prediction",
                        help="<path to an existing FILE> [mandatory]",
                        dest="hmi_prediction",
                        action="store",
                        required=True)

    parser.add_argument("--fasta_file",
                        help="<path to an existing FILE> [mandatory]",
                        dest="fasta_file",
                        action="store",
                        required=True)

    parser.add_argument("--resources",
                        help="<path to resources FILES> [mandatory]",
                        dest="resources",
                        action="store",
                        required=True)

    parser.add_argument("--results",
                        help="<path to results> [mandatory]",
                        dest="results",
                        action="store",
                        required=True)

    parser.add_argument("--output",
                        help="<path to output> [mandatory]",
                        dest="output",
                        action="store",
                        required=True)

    results = parser.parse_args(argv)
    return results


def process_hmi(file):
    """ Opening hmi file and select the human interactors """
    with open(file) as hmi_table:
        human_proteins = []
        for line in hmi_table:
            line = line.strip().split(";")
            protein = line[0]
            if protein not in human_proteins:
                human_proteins.append(protein)
    return human_proteins


def get_interaction(file):
    """ Opening hmi file and select the human interactors """
    with open(file) as hmi_table:
        hmi = []
        for line in hmi_table:
            line = line.strip().split(";")
            hmi.append(line)
    return hmi


def get_motif(file):
    """ Get information about bacteria interacting human motifs """
    with open(file) as hmi_table:
        motif = {}
        for line in hmi_table:
            line = line.strip().split(";")
            if line[0] not in motif:
                motif[line[0]] = set()
            motif[line[0]].add((line[1], line[2], line[3]))

    return motif

def split_large_fasta_file(fasta_file, protein_list, folder):
    """Creating small fasta files from a large one - not downloading one by one"""
    sequence_folder = os.path.join(folder, "protein_sequences")
    os.makedirs(sequence_folder, exist_ok=True)
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        if seq.description.split('|')[1] in protein_list:
            filename = os.path.join(sequence_folder, seq.description.split('|')[1] + ".fasta")
            #filename = sequence_folder + "\\" + seq.description.split('|')[1] + ".fasta"
            with open(f"{filename}", "w") as small_fasta:
                SeqIO.write(seq, small_fasta, "fasta")

def run_iupred(folder, method_type='short', anchor=True):
    """ IUPred modelling with an optional usage of ANCHOR2 """
    motif_score = {}
    sequence_folder = os.path.join(folder, "protein_sequences")
    os.makedirs(sequence_folder, exist_ok=True)
    iupred_data_folder = os.path.join(folder, "iupred_data")
    os.makedirs(iupred_data_folder, exist_ok=True)
    #iupred_data_folder = os.path.dirname(iupred_data_folder)

    for sequence in os.listdir(sequence_folder):
        protein_name = str(sequence).split(".")[0]
        sequence = read_seq(os.path.join(sequence_folder, sequence))
        iupred2_result = iupred(iupred_data_folder, sequence, method_type)
        if anchor:
            if method_type == 'long':
                anchor2_res = anchor2(iupred_data_folder, sequence, iupred2_result[0])
            else:
                anchor2_res = anchor2(iupred_data_folder, sequence, iupred(iupred_data_folder, sequence, 'long')[0])
        if method_type == 'glob':
            motif_score[protein_name] = iupred2_result[1]
        for pos, residue in enumerate(sequence):
            output_list = [protein_name, str(pos + 1), residue, str(iupred2_result[0][pos])]
            if anchor:
                output_list = [protein_name, str(pos + 1), residue, str(iupred2_result[0][pos]), str(anchor2_res[pos])]
            if protein_name not in motif_score:
                motif_score[protein_name] = []
            motif_score[protein_name].append(output_list[1:])


    return motif_score


def motif_selection(aa_scores, motif_dictionary):
    """ Selection of disordered motifs """
    disordered_motifs = []
    for protein in aa_scores:
        if protein in motif_dictionary:
            for list_ in motif_dictionary[protein]:
                motif = range(int(list_[1]), int(list_[2])+1)
                motif_size = len(motif)
                dis_aa_count = 0
                for details in aa_scores[protein]:
                    if int(details[0]) in motif:
                        if (float(details[2]) > 0.5) and (float(details[3]) > 0.5):
                            dis_aa_count += 1
                if (dis_aa_count == motif_size) or (dis_aa_count == motif_size - 1) or (dis_aa_count == motif_size + 1):
                    disordered_motifs.append((protein, list_[0], list_[1], list_[2]))
    return disordered_motifs


def write_output(folder, idr_motifs, hmi, output, anchor=True):
    """ Writing output file """
    output_file = os.path.join(folder, output)
    written_entries = set()

    with open(output_file, 'w') as output_file:
        output_file.write("# Human Protein" + "\t" + "Motif" + "\t" +  "Start" + "\t" +  "End" +
                          "\t" +"Bacterial domain" + "\t" + "Bacterial protein" "\n")
        idr_motifs = list(set(idr_motifs))

        for motif in idr_motifs:

            for interaction in hmi:
                if motif[0] in interaction and motif[1] in interaction and motif[2] in interaction and motif[3] in interaction:
                    entry = "\t".join(interaction)
                    if entry not in written_entries:
                        output_file.write(entry + "\n")
                        written_entries.add(entry)

def delete_files_in_folder(folder_path):
    try:
        # Check if the folder exists
        if os.path.exists(folder_path):
            # Iterate over all files in the folder and delete them
            for filename in os.listdir(folder_path):
                file_path = os.path.join(folder_path, filename)
                if os.path.isfile(file_path):
                    os.remove(file_path)
            print(f"All files in {folder_path} have been deleted.")
        else:
            print(f"The folder {folder_path} does not exist.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def main(argv):
    """ Main method and logic """

    # Read args
    args = parse_args(argv)

    # Get human interactors from the HMI prediction
    human_proteins = process_hmi(args.hmi_prediction)
    hmi = get_interaction(args.hmi_prediction)

    # Get interacting motif information
    motif_dict = get_motif(args.hmi_prediction)

    # Download fasta files - 1 protein/file
    split_large_fasta_file(args.fasta_file, human_proteins, args.resources)

    # Assign IUPred and ANCHOR scores to AAs
    scores = run_iupred(args.resources, 'short')

    # Select disordered motifs
    idr_motifs = motif_selection(scores, motif_dict)

    # Write the output file
    write_output(args.results, idr_motifs, hmi, args.output)

    # Remove files from protein_sequences folder
    folder_to_delete_files = os.path.join(args.resources, "protein_sequences")
    delete_files_in_folder(folder_to_delete_files)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
