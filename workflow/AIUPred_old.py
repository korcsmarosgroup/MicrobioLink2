# Author: Lejla Gul, Toby Lawrence
# Date: May, 2025

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
import numpy as np
from Bio import SeqIO
import torch
from scipy.signal import savgol_filter
from torch.nn.functional import pad  # needed by custom binding()

# ── local ────────────────────────────────────────────────────────────────────
import aiupred_lib

warnings.filterwarnings("ignore")
AA_CODE = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']
WINDOW = 100
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

# def binding_transform_v2(folder, prediction, smoothing=True):
#     transform = {}
#     iupred_data_folder = os.path.join(folder, "aiupred_data")
#     with open(f'{iupred_data_folder}/binding_transform') as fn:
#         for line in fn:
#             key, value = line.strip().split()
#             transform[int(float(key) * 1000)] = float(value)
#     rounded_pred = np.rint(prediction * 1000)
#     transformed_pred = np.vectorize(transform.get)(rounded_pred)
#     if not smoothing:
#         return transformed_pred
#     pred = savgol_filter(transformed_pred, 11, 5)
#     pred[pred > 1] = 1
#     return pred

# ═════════════════════════════════════════════════════════════════════════════
# AIUPred wrappers
# ═════════════════════════════════════════════════════════════════════════════


# @torch.no_grad()
# def tokenize(sequence, device):
#     """
#     Tokenize an amino acid sequence. Non-standard amino acids are treated as X
#     :param sequence: Amino acid sequence in string
#     :param device: Device to run on. CUDA{x} or CPU
#     :return: Tokenized tensors
#     """
#     return torch.tensor([AA_CODE.index(aa) if aa in AA_CODE else 20 for aa in sequence], device=device)

def initialise_models(resources_folder, which, device, gpu_num, force_cpu):
    data_dir = os.path.join(resources_folder, 'aiupred_data')
    if which == 'disorder':
        emb = aiupred_lib.TransformerModel()
        emb.load_state_dict(torch.load(os.path.join(data_dir, 'embedding_disorder.pt'),
                                       map_location=device))
        dec = aiupred_lib.DecoderModel()
        dec.load_state_dict(torch.load(os.path.join(data_dir, 'disorder_decoder.pt'),
                                       map_location=device))
    else:
        emb = aiupred_lib.BindingTransformerModel()
        emb.load_state_dict(torch.load(os.path.join(data_dir, 'embedding_binding.pt'),
                                       map_location=device))
        dec = aiupred_lib.BindingDecoderModel()
        dec.load_state_dict(torch.load(os.path.join(data_dir, 'binding_decoder.pt'),
                                       map_location=device))

    emb.to(device).eval()
    dec.to(device).eval()
    return emb, dec

def binding_transform_v2(resources_folder, prediction, smoothing=True):
    tf_path = os.path.join(resources_folder, 'aiupred_data', 'binding_transform')
    table = {}
    with open(tf_path) as fh:
        for line in fh:
            k, v = line.strip().split()
            table[int(float(k) * 1000)] = float(v)

    rounded = np.rint(prediction * 1000)
    transformed = np.vectorize(table.get)(rounded)
    if not smoothing:
        return transformed
    pred = savgol_filter(transformed, 11, 5)
    pred[pred > 1] = 1
    return pred

def predict_binding_v2(resources_folder, seq, embed_model, dec_model, device,
                    smoothing=True, binding=True):
    """AIUPred‑binding per‑residue profile (0‑1).
    Re‑implements library call so we can inject our own transform path."""
    tokens = aiupred_lib.tokenize(seq, device)
    padded = pad(tokens, (WINDOW//2, WINDOW//2), value=0)
    unfold = padded.unfold(0, WINDOW+1, 1)

    emb = embed_model(unfold, embed_only=True)
    emb_pad = pad(emb, (0, 0, 0, 0, WINDOW//2, WINDOW//2), value=0)
    emb_unf = emb_pad.unfold(0, WINDOW+1, 1)
    decoder_in = emb_unf.permute(0, 2, 1, 3)
    pred = dec_model(decoder_in).detach().cpu().numpy()

    if binding:
        print(binding_transform_v2(resources_folder, pred, smoothing=smoothing))
        return binding_transform_v2(resources_folder, pred, smoothing=smoothing)
    if smoothing and len(seq) > 10:
        return savgol_filter(pred, 11, 5)
    return pred




# ═════════════════════════════════════════════════════════════════════════════
# Motif filtering
# ═════════════════════════════════════════════════════════════════════════════

def predict_tracks(folder, seq_dict, device, gpu_num, force_cpu):
    """Return two dicts: disorder[pid], binding[pid]  (NumPy 1-D arrays)."""
    dis_e, dis_d = initialise_models(
        folder, 'disorder', device=device, gpu_num=gpu_num, force_cpu=force_cpu)

    # binding models
    bin_e, bin_d = initialise_models(
        folder, 'binding',  device=device, gpu_num=gpu_num, force_cpu=force_cpu)

    disorder_profiles, binding_profiles = {}, {}
    print(seq_dict.items())
    print(len(seq_dict.items()))
    for pid, seq in seq_dict.items():
        disorder_profiles[pid] = aiupred_lib.predict_disorder(
            seq, dis_e, dis_d, device, smoothing=True)

        binding_profiles[pid] = predict_binding_v2(
            folder, seq, bin_e, bin_d, device,
            smoothing=True, binding=True)      # 0–1 range

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
            if (d_prof[s:e+1] >= thr).all() and (b_prof[s:e+1] >= thr).all():
                kept.append((pid, mseq, s, e))
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
        output_file.write("# Human Protein" + "\t" + "Motif" + "\t" +  "Start" + "\t" +  "End" +
                          "\t" +"Bacterial domain" + "\t" + "Bacterial protein" "\n")

        idr_motifs = list(set(idr_motifs))

        for motif in idr_motifs:
            for interaction in hmi:
                if motif[0] in interaction and motif[1] in interaction and str(motif[2]) in interaction and str(motif[3]) in interaction:
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
    print(len(seqs))

    # 3) device selection
    device = torch.device('cpu') if args.force_cpu or not torch.cuda.is_available() \
             else torch.device(f'cuda:{args.gpu}')
    logging.info('Using device %s', device)

    # 4) predict tracks
    logging.info('Running AIUPred on %d proteins', len(seqs))
    disorder, binding = predict_tracks(args.resources, seqs, device, gpu_num=0, force_cpu=False)

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
