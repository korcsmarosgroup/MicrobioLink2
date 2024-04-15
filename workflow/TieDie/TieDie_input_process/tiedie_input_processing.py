import argparse
import omnipath as op
import pandas as pd
import numpy as np
import csv
import math

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process input files.')
    parser.add_argument('--transcriptomics_file', required=True, help='Path to transcriptomics file')
    parser.add_argument("--value_column", help="Column number for expression values", action="store", type=int, required=True)
    parser.add_argument("--sep_transcriptomics", help="Separator for the transcriptomic file", dest="sep_transcriptomics", action="store", required=True)
    parser.add_argument('--endpoint_file', required=True, help='Path to endpoint file')
    parser.add_argument('--endpoint_pvalue_column', required=False, help='Column number for the adjusted p-value/FDR value in the differentially expressed gene file', type=int)
    parser.add_argument('--endpoint_value_column', required=False, help='Column number for log2FC/Expression value in the endpoint file', type=int)
    parser.add_argument("--sep_endpoint", help="Separator for the endpoint file", required=True)
    parser.add_argument('--hmi_prediction_file', required=True, help='Path to hmi prediction output file')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory (resource)')
    parser.add_argument("--upstream_input_filename", type=str, default="upstream.input", help="Custom name for upstream output file")
    parser.add_argument("--downstream_input_filename", type=str, default="downstream.input", help="Custom name for downstream output file")
    parser.add_argument("--pathway_input_filename", type=str, default="pathway.sif", help="Custom name for pathway output file")

    return parser.parse_args()


def read_transcriptomics_file(file_path, sep, value_column):
    genes = []
    with open(file_path) as transcriptomics:
        transcriptomics.readline()
        for line in transcriptomics:
            line = line.strip().split(sep)
            value_column_new = int(value_column) - 1
            value = line[value_column_new]
            # Check if the value is a valid float and not equal to 0.0
            if value and not math.isnan(float(value)) and float(value) != 0.0:
                genes.append(line[0])
    return genes


def filter_ppi_by_genes(ppi, genes):
    return ppi[ppi['source_genesymbol'].isin(genes) & ppi['target_genesymbol'].isin(genes)]


def filter_tf_tg(tf_tg, genes, endpoint_genes):
    contextulised_tf_tg = tf_tg[tf_tg['source_genesymbol'].isin(genes) & tf_tg['target_genesymbol'].isin(endpoint_genes[endpoint_genes.columns[0]])]
    return contextulised_tf_tg.drop_duplicates()


def get_regulators_info(contextualised_tf_tg):
    regs = contextualised_tf_tg.pivot_table(columns=['source'], aggfunc='size')
    regs = pd.DataFrame(regs).reset_index().rename(columns={"source": "TF_name", 0: "num_degs"})
    return regs

def save_results(contextualised_tf_tg, regs, output_dir):
    contextualised_tf_tg.to_csv(f"{output_dir}/contextualised_regulator-target_network.txt", sep="\t", index=False)
    regs.to_csv(f"{output_dir}/contextualised_regulators_of_targets.txt", sep="\t", index=False)

def create_sif(ppi, output_dir, pathway_input_filename):
    ppis2 = ppi[['source', 'consensus_stimulation', 'target']].copy()
    ppis2.loc[ppis2['consensus_stimulation'] == 1, 'consensus_stimulation'] = 'stimulates>'
    ppis2.loc[ppis2['consensus_stimulation'] == 0, 'consensus_stimulation'] = 'inhibits>'
    ppis2 = ppis2.drop_duplicates().rename(columns={'consensus_stimulation': 'direction'})
    ppis2.to_csv(f"{output_dir}/{pathway_input_filename}", sep='\t', index=None, header=False)

def read_hmi_file(file_path):
    hbps = pd.read_csv(file_path, sep="\t", index_col=False)
    return hbps[['# Human Protein', 'Bacterial protein']].copy()

def process_upstream_input(hbps, output_dir, upstream_input_filename):
    if 'sign' in hbps.columns:
        hbps_f = hbps[(hbps['sign'] == '-') | (hbps['sign'] == '+')]
        if len(hbps) != len(hbps_f):
            print("WARNING: Some of the bacterial-human binding protein interactions were discarded as the values in the 'sign' column were not '+' or '-'.")
        if len(hbps_f) == 0:
            print("ERROR: bacterial-human binding protein interactions do not have the correct values in 'sign' column. They should be '+' or '-'.")
        hbps2 = hbps[['# Human Protein', 'sign']].copy().rename(columns={'sign': 'direction'})
        hbps2 = hbps2.groupby(['# Human Protein', 'direction']).size().reset_index(name='n')[['# Human Protein', 'n', 'direction']]
    else:
        hbps2 = hbps[['# Human Protein']].copy().groupby('# Human Protein').size().reset_index(name='n')
        hbps2['direction'] = '-'
    hbps2.to_csv(f"{output_dir}/{upstream_input_filename}", sep='\t', index=None, header=False)

def process_downstream_input(contextualised_tf_tg, genes, output_dir, downstream_input_filename, endpoint_value_column):
    expression = str(genes.columns[int(endpoint_value_column)-1])
    tfs2 = contextualised_tf_tg.merge(genes, on='target_genesymbol')
    tfs2_filtered = tfs2[['source', 'target', 'consensus_stimulation', expression]].copy()
    tfs2_filtered['exp_sign'] = np.where(tfs2_filtered['consensus_stimulation'] == True, tfs2_filtered[expression], -1*(tfs2_filtered[expression]))
    tfs3 = tfs2_filtered[['source', 'target', 'exp_sign']].copy()
    tfs3_1 = tfs3.groupby('source')
    tfs4 = tfs3_1.count()
    tfs4['sumof'] = tfs3_1.exp_sign.sum()
    tfs4['final_val'] = tfs4['sumof'] / tfs4['exp_sign']
    tfs4['sign'] = np.where(tfs4['final_val'] >= 0, '+', '-')
    tfs4 = tfs4.drop(['exp_sign', 'sumof', 'target'], axis=1)

    tfs4.to_csv(f"{output_dir}/{downstream_input_filename}", sep='\t', header=False)

def main():
    args = parse_arguments()

    ppis = op.interactions.OmniPath.get(genesymbols=1)
    tf_tg = op.interactions.Transcriptional.get(databases='CollecTRI', genesymbols=1)

    genes = read_transcriptomics_file(args.transcriptomics_file, args.sep_transcriptomics, args.value_column)
    contextulised_ppi = filter_ppi_by_genes(ppis, genes)

    endpoint_genes = pd.read_csv(args.endpoint_file, sep=args.sep_endpoint)
    endpoint_genes = endpoint_genes.rename(columns={endpoint_genes.columns[0]: 'target_genesymbol'})
    if args.endpoint_pvalue_column:
        padj_column = int(args.endpoint_pvalue_column) - 1
        endpoint_genes = endpoint_genes[endpoint_genes.iloc[:, padj_column] < 0.05]

    contextulised_tf_tg = filter_tf_tg(tf_tg, genes, endpoint_genes)
    regs = get_regulators_info(contextulised_tf_tg)

    save_results(contextulised_tf_tg, regs, args.output_dir)
    create_sif(contextulised_ppi, args.output_dir, args.pathway_input_filename)

    hbps = read_hmi_file(args.hmi_prediction_file)
    process_upstream_input(hbps, args.output_dir, args.upstream_input_filename)

    process_downstream_input(contextulised_tf_tg, endpoint_genes, args.output_dir, args.downstream_input_filename, args.endpoint_value_column)

if __name__ == "__main__":
    main()
