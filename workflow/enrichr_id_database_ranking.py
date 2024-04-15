import argparse
from mygene import MyGeneInfo
import gget
import matplotlib.pyplot as plt  # Import Matplotlib
import numpy as np

#Parse command-line arguments.
def parse_arguments():
    parser = argparse.ArgumentParser(description="Your script description here.")
    parser.add_argument("--background_gene_list", required=True, help="Path to the file background gene list if applicable.")
    parser.add_argument("--target_gene_list", required=True, help="Path to the target gene list file.")
    parser.add_argument("--output_image", required=True, help="Path to the output image.")
    parser.add_argument("--output_file", required=True, help="Path to the output file.")
    parser.add_argument("--database", default="Reactome_2022", help="""Database to use as a reference for the enrichment analysis. 
                                Database to use as reference for the enrichment analysis.
                        Supported shortcuts (and their default database):
                        'pathway' (KEGG_2021_Human)
                        'transcription' (ChEA_2016)
                        'ontology' (GO_Biological_Process_2021)
                        'diseases_drugs' (GWAS_Catalog_2019)
                        'celltypes' (PanglaoDB_Augmented_2021)
                        'kinase_interactions' (KEA_2015)
                        or any database listed under Gene-set Library at: https://maayanlab.cloud/Enrichr/#libraries""")
    parser.add_argument("--ranking", default="combined_score", help="Plot will be ranked based on \"combined_score\" or \"adj_p\".")
                                
                                
    return parser.parse_args()

#Read the background gene list from a file.
def read_background_gene_list(file_path):
    background_genes = []
    with open(file_path, "r") as background_gene_list:
        for line in background_gene_list:
            line = line.strip().split(',')
            background_genes.append(line[0])
    return background_genes

#Read the target gene list from a file.
def read_target_gene_list(file_path):
    target_geneset = []
    with open(file_path, "r") as target_gene_list:
        for line in target_gene_list:
            line = line.strip().split("\t")
            target_geneset.extend([line[0], line[2]])
    target_geneset = list(set(target_geneset))
    return target_geneset

#Translate UniProt IDs to gene symbols.
def translate_uniprot_to_symbols(uniprot, species='human'):
    mg = MyGeneInfo()
    target_genesymbols_translation = mg.querymany(uniprot, scopes='uniprot', fields='symbol', species='human', returnall=True)
    
    translation_dict = {entry['query']: entry.get('symbol', None) for entry in target_genesymbols_translation['out']}
    return translation_dict

def main():
    args = parse_arguments()

    background_expressed_gene_symbols = read_background_gene_list(args.background_gene_list)
    target_geneset = read_target_gene_list(args.target_gene_list)

    target_genesymbols = translate_uniprot_to_symbols(target_geneset)
    only_target_genesymbols = list(target_genesymbols.values())

    if args.ranking == "combined_score":

        enrichr_df_bkg = gget.enrichr(only_target_genesymbols, database=args.database, background_list=background_expressed_gene_symbols)

        #Filtering results based on adjusted p-values.
        filtered_enrichr_df_bkg = enrichr_df_bkg[enrichr_df_bkg["adj_p_val"] < 0.05]

        #Writing output file.
        filtered_enrichr_df_bkg.to_csv(args.output_file, index=False)
    
        #Create bar plot ranked by combined score
        filtered_enrichr_df_bkg = filtered_enrichr_df_bkg.sort_values('combined_score', ascending = False)
        filtered_enrichr_inf = filtered_enrichr_df_bkg[filtered_enrichr_df_bkg['combined_score'].astype(str) != 'inf']
    
        top_20_functions = filtered_enrichr_inf.head(20).sort_values('combined_score')
    
        fig, ax1 = plt.subplots()

        ax1.barh(y= top_20_functions['path_name'], width= top_20_functions['combined_score'])
        ax1.set_xlabel('Combined score')
        ax1.set_title('Top 20 enriched pathways by combined score')

        #Add adjusted p-values to the plot
        ax2 = ax1.twiny()
        ax2.scatter(x= -np.log10(top_20_functions['adj_p_val']), y= top_20_functions['path_name'], color= '#FF8800')
        ax2.set_xlabel('-log 10(Adjusted p-value)', color= '#FF8800')
        ax2.tick_params('x', colors= '#FF8800')

        fig.savefig(args.output_image, dpi= 300,bbox_inches= 'tight', transparent= True)

    elif args.ranking == "adj_p":

        enrichr_df_bkg = gget.enrichr(only_target_genesymbols, database=args.database, background_list=background_expressed_gene_symbols, plot=True)
                
        #Filtering results based on adjusted p-values.        
        filtered_enrichr_df_bkg = enrichr_df_bkg[enrichr_df_bkg["adj_p_val"] < 0.05]

        #Writing output file.
        filtered_enrichr_df_bkg.to_csv(args.output_file, index=False)

        #Create bar plot ranked by adjusted p-value (default in gget enrichr).
        plt.savefig(args.output_image, dpi=300, bbox_inches="tight", transparent=True)

if __name__ == "__main__":
    main()