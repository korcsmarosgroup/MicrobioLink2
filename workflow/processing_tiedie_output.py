# Author : Lejla Gul
# Date : January 2024
#
# Script to combine TieDie outputs into 1 network table and one node table
#
# Input: Receptor-tf network output from tiedie
#        Bacteria - human binding protein network (columns 'bacterial_protein', 'human_protein')
#        Contextualised regulatory network (TFs - DEGs) (incl. columns 'DEG',"	TF", "consensus_stimulation")
#        Heats values output from TieDie
#        Bacterial proteins gene symbol to uniprot conversion table
#        Enfpoint gene (differential gene expression) table
#
# Output: Network file where each line represents an interaction
#         Node table where each lines represents one node in the network

import pandas as pd
import numpy as np
from mygene import MyGeneInfo
import argparse

def load_data(tiedie_network, tiedie_heats,hmi,tf_tg_network,endpoint):
    rec_tf = pd.read_csv(tiedie_network, header=None, names=["Source.node","Relationship","Target.node"], sep = "\t")

    heats = pd.read_csv(tiedie_heats, sep="=")
    heats = heats.reset_index()
    heats['index'] = heats['index'].apply(lambda x: x.replace(' ', ''))

    hbps = pd.read_csv(hmi, sep = "\t")
    hbps = hbps[['# Human Protein', 'Bacterial protein']]
    hbps['sign'] = '+'

    bacteria = hbps[['Bacterial protein']].copy()
    bacteria['Annotation'] = hbps['Bacterial protein']

    tf_deg = pd.read_csv(tf_tg_network, sep = "\t")
    tf_deg['consensus_stimulation'] = tf_deg['consensus_stimulation'].replace({True: 'stimulates>', False: 'inhibits>'})

    endpoint = pd.read_csv(endpoint)
    endpoint = endpoint[endpoint['adj.P.Val'] < 0.05]
    endpoint = endpoint[['Gene', 'logFC']]

    return rec_tf, heats, hbps, bacteria, tf_deg, endpoint

def process_edges(rec_tf, hbps, tf_deg, network_output):
    rec_tf['layer'] = 'bindingprot-tf'
    rec_tf = rec_tf[['Target.node','Source.node','Relationship','layer']]

    if 'sign' in hbps.columns:
        conditions = [hbps['sign'] == '-', hbps['sign'] == '+']
        choices = ['inhibits>', 'stimulates>']
        hbps['Relationship'] = np.select(conditions, choices, default='unknown')
        hbps2 = hbps.drop(columns=['sign'])

        hbps2['layer'] = 'bacteria-bindingprot'
        hbps2 = hbps2.rename(columns={'Bacterial protein':'Source.node', '# Human Protein': 'Target.node'})

        hbps2_filtered = hbps2[~hbps2['Source.node'].isin(rec_tf['Source.node']) & hbps2['Target.node'].isin(rec_tf['Source.node'])]
    else:
        hbps['Relationship'] = 'unknown'
        hbps['layer'] = 'bacteria-bindingprot'
        hbps2 = hbps.rename(columns={'bacterial_protein':'Baterial protein', 'human_protein':'# Human Protein'})

        hbps2_filtered = hbps2[~hbps2['Source.node'].isin(rec_tf['Source.node']) & hbps2['Target.node'].isin(rec_tf['Source.node'])]

    tf_deg_filtered = tf_deg[tf_deg['source'].isin(rec_tf['Target.node'])]
    tf_deg_filtered['layer'] = "tf-deg"
    tf_deg2 = tf_deg_filtered[['source','target','consensus_stimulation','layer']].copy()
    tf_deg2 = tf_deg2.rename(columns={'source':'Source.node','target':'Target.node','consensus_stimulation':'Relationship'})

    dfs = [hbps2_filtered, rec_tf, tf_deg2]

    whole_net = pd.concat(dfs)

    whole_net.to_csv(network_output, sep="\t", index=None)
    return whole_net

def process_node_table(whole_net):
    hmi = whole_net[(whole_net['layer'] == 'bacteria-bindingprot')]
    hmi.insert(4,'bacteria_layer', 'bacteria')
    bacteria = hmi[['Source.node','bacteria_layer']].copy()
    bacteria = bacteria.drop_duplicates()

    hmi2 = whole_net[(whole_net['layer'] == 'bacteria-bindingprot')]
    hmi2.insert(4,'bindingprot_layer', 'bindingprot')
    binding_prot = hmi2[['Target.node','bindingprot_layer']].copy()
    binding_prot = binding_prot.drop_duplicates()

    ppi1 = whole_net[(whole_net['layer'] == 'bindingprot-tf')]
    ppi1.insert(4,'ppi1_layer', 'bindingprot and/or protein')
    source_prot = ppi1[['Source.node','ppi1_layer']].copy()
    source_prot= source_prot.drop_duplicates()

    ppi2 = whole_net[(whole_net['layer'] == 'bindingprot-tf')]
    ppi2.insert(4,'ppi2_layer', 'protein and/or tf')
    target_prot = ppi2[['Target.node','ppi2_layer']].copy()
    target_prot = target_prot.drop_duplicates()

    tf_tg = whole_net[(whole_net['layer'] == 'tf-deg')]
    tf_tg.insert(4,'tf_layer', 'tf')
    tf = tf_tg[['Source.node','tf_layer']].copy()
    tf = tf.drop_duplicates()

    tf_tg2 = whole_net[(whole_net['layer'] == 'tf-deg')]
    tf_tg2.insert(4,'deg_layer', 'deg')
    deg = tf_tg2[['Target.node','deg_layer']].copy()
    deg = deg.drop_duplicates()

    source_df = pd.merge(whole_net[['Source.node']], bacteria, how='left', left_on='Source.node', right_on='Source.node')

    binding_prot = binding_prot.rename(columns={'Target.node':'Source.node'})
    source_df = pd.merge(source_df, binding_prot, how='left', left_on='Source.node', right_on='Source.node')

    source_df = pd.merge(source_df, source_prot, how='left', left_on='Source.node', right_on='Source.node')

    target_prot = target_prot.rename(columns={'Target.node':'Source.node'})
    source_df = pd.merge(source_df, target_prot, how='left', left_on='Source.node', right_on='Source.node')

    source_df = pd.merge(source_df, tf, how='left', left_on='Source.node', right_on='Source.node')

    deg = deg.rename(columns={'Target.node':'Source.node'})
    source_df = pd.merge(source_df, deg, how='left',left_on='Source.node', right_on='Source.node')


    source_df['ppi_layer'] = source_df['ppi1_layer'].astype(str) + ',' + source_df['ppi2_layer'].astype(str)
    source_df = source_df.drop(['ppi1_layer', 'ppi2_layer'], axis=1)
    source_df = source_df.fillna('NA')

    bacteria = bacteria.rename(columns={'Source.node':'Target.node'})
    target_df = pd.merge(whole_net[['Target.node']], bacteria, how='left', left_on='Target.node', right_on='Target.node')

    binding_prot = binding_prot.rename(columns={'Source.node':'Target.node'})
    target_df = pd.merge(target_df, binding_prot, how='left', left_on='Target.node', right_on='Target.node')

    source_prot = source_prot.rename(columns={'Source.node':'Target.node'})
    target_df = pd.merge(target_df, source_prot, how='left', left_on='Target.node', right_on='Target.node')

    target_prot = target_prot.rename(columns={'Source.node':'Target.node'})
    target_df = pd.merge(target_df, target_prot, how='left', left_on='Target.node', right_on='Target.node')

    tf = tf.rename(columns={'Source.node':'Target.node'})
    target_df = pd.merge(target_df, tf, how='left', left_on='Target.node', right_on='Target.node')

    deg = deg.rename(columns={'Source.node':'Target.node'})
    target_df = pd.merge(target_df, deg, how='left',left_on='Target.node', right_on='Target.node')

    target_df['ppi_layer'] = target_df['ppi1_layer'].astype(str) + ',' + target_df['ppi2_layer'].astype(str)
    target_df = target_df.drop(['ppi1_layer', 'ppi2_layer'], axis=1)
    target_df = target_df.fillna('NA')
    target_df = target_df.rename(columns={'Target.node':'Source.node'})
    merged_df = pd.concat([source_df, target_df], ignore_index=True)
    merged_df['all_nodes'] = merged_df[['bacteria_layer', 'bindingprot_layer', 'ppi_layer', 'tf_layer', 'deg_layer']].apply(lambda x: ','.join(filter(lambda val: val != 'NA', x)), axis=1)
    merged_df = merged_df.rename(columns={'Source.node':'node'})
    merged_df = merged_df[['node','all_nodes', 'bacteria_layer', 'bindingprot_layer', 'ppi_layer', 'tf_layer', 'deg_layer']]
    merged_df = merged_df.drop_duplicates()

    return merged_df

def batch_uniprot_to_gene_symbol(uniprot_ids):
    try:
        mg = MyGeneInfo()
        result = mg.querymany(uniprot_ids, scopes='uniprot', fields='symbol', species='human', returnall=True)
        mapping_dict = {entry['query']: entry.get('symbol', None) for entry in result['out']}
        return mapping_dict
    except Exception as e:
        print(f"Error translating UniProt IDs to gene symbols: {e}")
        return None


def fetch_gene_symbols(merged_df):
    uniprot_ids = merged_df['node'].tolist()
    batch_size = 100
    batches = [uniprot_ids[i:i + batch_size] for i in range(0, len(uniprot_ids), batch_size)]
    gene_symbol_mapping = {}

    for batch in batches:
        batch_mapping = batch_uniprot_to_gene_symbol(batch)
        if batch_mapping:
            gene_symbol_mapping.update(batch_mapping)

    merged_df['gene_symbol'] = merged_df['node'].map(gene_symbol_mapping)
    merged_df['gene_symbol'] = merged_df['gene_symbol'].fillna(merged_df['node'])

    return merged_df

def main():
    parser = argparse.ArgumentParser(description='Process network data and create node table.')
    parser.add_argument('--tiedie_file', type=str, default='tiedie.cn.sif', help='Path to tiedie.cn.sif file')
    parser.add_argument('--heats_file', type=str, default='heats.NA', help='Path to heats.NA file')
    parser.add_argument('--hmi_file', type=str, help='Path to the host-microbe interaction file')
    parser.add_argument('--tf_tg_file', type=str, help='Path to contextualised regulatory network file')
    parser.add_argument('--endpoint_file', type=str, help='Path to the endpoint gene list/DEG file')
    parser.add_argument('--network_output', type=str, help='Output path and name for the network file')
    parser.add_argument('--node_attr_output', type=str, help='Output path and name for the node attribute file')

    args = parser.parse_args()

    rec_tf, heats, hbps, _, tf_deg, exp = load_data(args.tiedie_file, args.heats_file, args.hmi_file, args.tf_tg_file, args.endpoint_file)

    whole_net = process_edges(rec_tf, hbps, tf_deg, args.network_output)

    merged_df = process_node_table(whole_net)

    merged_df = fetch_gene_symbols(merged_df)

    heats = heats.rename(columns={'index':'node'})

    merged_df = pd.merge(merged_df, heats, how='left', on='node')

    exp = exp.rename(columns={'Gene':'gene_symbol'})

    merged_df = pd.merge(merged_df, exp, how='left', on='gene_symbol')

    merged_df.to_csv(args.node_attr_output, sep="\t", index=None)

if __name__ == "__main__":
    main()
