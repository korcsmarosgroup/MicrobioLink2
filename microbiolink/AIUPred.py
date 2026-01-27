# Author: Lejla Gul, Toby Lawrence
# Date: Sep, 2025

# Aim: Filter host-micere protein-protein interactions, keeping only those where
# the target motif lies on residues that are both disordered
# (AIUPred‑disorder ≥ 0.60)and binding‑prone (AIUPred‑binding ≥ 0.60)

# Input: HMI table must be semicolon‑separated
# The script expects that AIUPred model files are present in
# resources/aiupred_data/

# Output: Filtered HMI table

#!/usr/bin/env python3

# ── standard library ──────────────────────────────────────────────────────────
import argparse, logging, os, sys, warnings

# ── third‑party ───────────────────────────────────────────────────────────────
from Bio import SeqIO
from iupred import (
    init_aiupred_models,
    predict_aiupred_binding,
    predict_aiupred_disorder,
)

warnings.filterwarnings("ignore")
THRESHOLD = 0.60  # score cutoff for both tracks

# ═════════════════════════════════════════════════════════════════════════════
# CLI helpers
# ═════════════════════════════════════════════════════════════════════════════

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

    parser.add_argument("-v", "--verbose",
                        help="Increase output verbosity",
                        action="store_true")

    parser.add_argument("--gpu",
                    help="Index of GPU to use, default=0",
                    default=0)

    parser.add_argument("--force-cpu",
                    help="Force the network to only utilize the CPU. Calculation will be very slow, not recommended",
                    action="store_true")

    return parser.parse_args(argv)

# ═════════════════════════════════════════════════════════════════════════════
# I/O helpers
# ═════════════════════════════════════════════════════════════════════════════

def process_hmi(hmi_path):
    """Return list of unique human protein IDs in the HMI table."""
    ids = []
    with open(hmi_path) as fh:
        for line in fh:
            pid = line.split(';', 1)[0].strip()
            if pid not in ids:
                ids.append(pid)
    return ids

def get_interaction(hmi_path):
    """ Opening hmi file and select the human interactors """
    with open(hmi_path) as hmi_table:
        hmi = []
        for line in hmi_table:
            line = line.strip().split(";")
            hmi.append(line)
    return hmi


def get_motif_dict(hmi_path):
    """Map protein‑ID → set{ (motif_seq, start, end) }."""
    motifs = {}
    with open(hmi_path) as fh:
        for line in fh:
            pid, mseq, start, end, *_ = line.strip().split(';')
            motifs.setdefault(pid, set()).add((mseq, start, end))
    return motifs


def split_large_fasta(multi_fa_path, proteins, out_folder):
    """Write <out_folder>/protein_sequences/<ID>.fasta files for requested IDs."""
    seq_dir = os.path.join(out_folder, 'protein_sequences')
    os.makedirs(seq_dir, exist_ok=True)

    want = set(proteins)
    for record in SeqIO.parse(multi_fa_path, 'fasta'):
        # UniProt ID is second field in description: >db|ID|...
        try:
            pid = record.description.split('|')[1]
        except IndexError:
            continue
        if pid in want:
            SeqIO.write(record, os.path.join(seq_dir, f'{pid}.fasta'), 'fasta')
            want.discard(pid)
        if not want:
            break  # all extracted


def read_seq(fasta_path):
    """Return the sequence of a 1‑sequence FASTA file as a string."""
    with open(fasta_path) as fh:
        return ''.join(line.strip() for line in fh if not line.startswith('>'))






# ═════════════════════════════════════════════════════════════════════════════
# Motif filtering
# ═════════════════════════════════════════════════════════════════════════════

def predict_tracks(seq_dict, force_cpu, gpu_num):
    """Return two dicts: disorder[pid], binding[pid]  (NumPy 1-D arrays)."""
    dis_e, dis_d, device = init_aiupred_models(
        'disorder', force_cpu=force_cpu, gpu_num=gpu_num)

    bin_e, bin_d, _ = init_aiupred_models(
        'binding', force_cpu=force_cpu, gpu_num=gpu_num)

    disorder_profiles, binding_profiles = {}, {}
    for pid, seq in seq_dict.items():
        disorder_profiles[pid] = predict_aiupred_disorder(
            seq, dis_e, dis_d, device, smoothing=True)

        binding_profiles[pid] = predict_aiupred_binding(
            seq, bin_e, bin_d, device,
            smoothing=True, binding=True)

    return disorder_profiles, binding_profiles

def keep_high_confidence(motif_dict, disorder, binding, thr=THRESHOLD):
    kept = []
    for pid, motifs in motif_dict.items():
        if pid not in disorder:
            continue
        d_prof = disorder[pid]
        b_prof = binding[pid]
        for mseq, start, end in motifs:

            s, e = int(start), int(end)  # inclusive indices
            length = len(d_prof[s:e+1])
            disordered_average = sum(d_prof[s:e+1]) / length
            binding_average = sum(b_prof[s:e+1]) / length
            combined_score = (binding_average + disordered_average) / 2

            if (d_prof[s:e+1] >= thr).all() and (b_prof[s:e+1] >= thr).all():
                kept.append((pid, mseq, s, e, length,disordered_average, binding_average, combined_score))

    return kept


# def write_output(rows, out_path):
#     os.makedirs(os.path.dirname(out_path), exist_ok=True)
#     with open(out_path, 'w') as fh:
#         fh.write('protein;motif_seq;start;end\n')
#         for r in rows:
#             fh.write(';'.join(map(str, r)) + '\n')

def write_output(folder, idr_motifs, hmi, output):
    """ Writing output file """
    output_file = os.path.join(folder, output)
    written_entries = set()

    with open(output_file, 'w') as output_file:
        output_file.write("Bacterial protein" + "\t" + "Bacterial domain" + "\t"
        + "Human Protein" + "\t" + "Motif" + "\t" +  "Start" + "\t" +  "End" +
        "\t" + "Motif length" + '\t' + 'Avg IUPred score'+ '\t' + 'Avg Binding score' + "\t" +'Combined score' + "\n")

        idr_motifs = list(set(idr_motifs))
        print(idr_motifs)

        for motif in idr_motifs:
            for interaction in hmi:
                if motif[0] in interaction and motif[1] in interaction and str(motif[2]) in interaction and str(motif[3]) in interaction:
                    #entry = "\t".join(interaction) + "\t" + str(motif[4]) + "\t" + str(motif[5]) + "\t" + str(motif[6]) + "\t" + str(motif[7])
                    entry = interaction[5] + "\t" + interaction[4] + "\t" + interaction[0] + "\t" + interaction[1] + "\t" + interaction[2] + "\t" + interaction[3] + "\t" + str(motif[4]) + "\t" + str(motif[5]) + "\t" + str(motif[6]) + "\t" + str(motif[7])

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

# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════

def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format='%(levelname)s: %(message)s')

    # 1) collect protein IDs & motif metadata
    logging.info('Reading HMI table %s', args.hmi_prediction)
    prot_ids   = process_hmi(args.hmi_prediction)
    hmi = get_interaction(args.hmi_prediction)
    motif_dict = get_motif_dict(args.hmi_prediction)

    # 2) write per‑protein FASTAs (if not already there)
    logging.info('Extracting %d host protein sequences', len(prot_ids))
    split_large_fasta(args.fasta_file, prot_ids, args.resources)

    seq_dir = os.path.join(args.resources, 'protein_sequences')

    seqs = {f[:-6]: read_seq(os.path.join(seq_dir, f))
            for f in os.listdir(seq_dir) if f.endswith('.fasta')}


    # 3) predict tracks
    logging.info('Running AIUPred on %d proteins', len(seqs))
    disorder, binding = predict_tracks(seqs, force_cpu=args.force_cpu, gpu_num=args.gpu)

    # 5) filter motifs
    logging.info('Filtering motifs with score ≥%.2f in both tracks', THRESHOLD)
    kept = keep_high_confidence(motif_dict, disorder, binding)
    logging.info('Kept %d motifs', len(kept))

    # 6) write output
    write_output(args.results, kept, hmi, args.output)
    logging.info('Results written to %s', args.output)


    # Remove files from protein_sequences folder
    folder_to_delete_files = os.path.join(args.resources, "protein_sequences")
    delete_files_in_folder(folder_to_delete_files)

if __name__ == '__main__':
    sys.exit(main())
