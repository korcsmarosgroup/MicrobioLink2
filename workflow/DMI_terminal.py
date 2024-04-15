import re
import argparse
from pyfasta import Fasta

def rename(fasta_key):
    """Extract UniProt ID from FASTA header."""
    fasta_key = fasta_key.split("|")
    fasta_key = fasta_key[1]
    return fasta_key

def parse_elm_regex(filename):
    """Parse the ELM regex from a file."""
    elm_regex = {}
    with open(filename, "r") as motif_table:
        motif_table.readline()
        for line in motif_table:
            if line[0] != '#':
                line = line.replace("\"", "")
                line = line.strip().split("\t")
                elm_regex[line[1]] = line[4]
    return elm_regex

def parse_motif_domain(filename):
    """Parse motif-domain interactions from a file."""
    motif_domain = {}
    with open(filename, "r") as motif_domain_table:
        motif_domain_table.readline()
        for line in motif_domain_table:
            line = line.replace("\"", "")
            line = line.strip("\n").split("\t")
            if len(line) > 1:
                if line[0] not in motif_domain:
                    motif_domain[line[0]] = []
                motif_domain[line[0]].append(line[1])
    return motif_domain


def parse_protein_domain(filename):
    """Parse protein-domain information from a file."""
    pfam_uniprot = {}
    with open(filename, "r") as protein_domain:
        protein_domain.readline()
        for line in protein_domain:
            line = line.strip().split("\t")
            if len(line) > 1:
                pfams = line[1].split(";")
                for pfam in pfams:
                    if pfam not in pfam_uniprot:
                        pfam_uniprot[pfam] = []
                    pfam_uniprot[pfam].append(line[0])
    return pfam_uniprot

def create_uniprot_motif_dict(human, elm_regex):
    """Create a dictionary mapping UniProt keys to motifs."""
    uniprot_motif = {}
    for key in human.keys():
        if rename(key) not in uniprot_motif:
            uniprot_motif[rename(key)] = []
        for motif in elm_regex:
            match = re.finditer(str(elm_regex[motif]), str(human[key]))
            for m in match:
                uniprot_motif[rename(key)].append((motif, str(m.start()), str(m.end())))
    return uniprot_motif

def main(args):

    # Load the human protein FASTA file
    human = Fasta(args.fasta_file)

    # Parse ELM regex
    elm_regex = parse_elm_regex(args.elm_regex_file)

    # Parse motif-domain interactions
    motif_domain = parse_motif_domain(args.motif_domain_file)

    # Parse protein domains
    pfam_uniprot = parse_protein_domain(args.bacterial_domain_file)

    # Create the uniprot_motif dictionary
    uniprot_motif = create_uniprot_motif_dict(human, elm_regex)


    # Create output file and write header
    with open(args.output_file, "w") as output:
        output.write('# Human Protein;Motif;Start;End;Bacterial domain;Bacteria Protein\n')
        for pfam, uniprot_list in pfam_uniprot.items():
            for uniprot in uniprot_list:
                for motif in motif_domain:
                    if pfam in motif_domain[motif]:
                        for uni, motif_list in uniprot_motif.items():
                            for motif_2 in motif_list:
                                if motif_2[0] == motif:
                                    output.write(f"{uni};{';'.join(motif_2)};{pfam};{uniprot}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict interaction between human and microbial proteins based on domain-motif interactions.")
    parser.add_argument("-fasta", "--fasta_file", help="Path to human protein FASTA file")
    parser.add_argument("-motif", "--elm_regex_file", help="Path to ELM regex file")
    parser.add_argument("-interaction", "--motif_domain_file", help="Path to motif-domain interaction file")
    parser.add_argument("-domain", "--bacterial_domain_file", help="Path to bacterial protein domain file")
    parser.add_argument("-o", "--output_file", help="Path to output file")
    args = parser.parse_args()
    main(args)
