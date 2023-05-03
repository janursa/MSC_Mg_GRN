import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_model_name_2_method_and_DE_type, TOP_LINKS_DIR, F_selected_models, GRN_DIR, change_protnames_to_genenames_in_links
from common_tools.links import choose_top_quantile
from uncertainity_analysis.ua_aux import NetworkNoiseAnalysis

def identify_top_links_around_target_genes(model_name, studies, top_quantile: float) -> None:
    """ Extract top links around the target genes
    """
    method, DE_type = F_model_name_2_method_and_DE_type(model_name)
    # - get the target genes
    genes_under_study = NetworkNoiseAnalysis.retrieve_robust_genes(model_name)
    # - get the links for ctr and sample
    links_ctr, links_sample = [pd.read_csv(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv')
                               for study in studies]
    links_ctr, links_sample = [change_protnames_to_genenames_in_links(links)
                               for links in [links_ctr, links_sample]]
    # - determine top links around each target gene
    links_short_stack = []
    for gene in genes_under_study:
        # - ctr and sample links around target gene
        gene_links_ctr, gene_links_sample = [links.loc[(links['Target'] == gene) | (links['Regulator'] == gene), :]
                                             for links in [links_ctr, links_sample]]
        # - calculate difference between ctr and sample
        links_diff = gene_links_ctr.copy()
        links_diff.drop('Weight', axis=1, inplace=True)
        weight_diff = abs(gene_links_ctr['Weight'] - gene_links_sample['Weight'])
        links_diff['Difference'] = weight_diff
        links_short = choose_top_quantile(links_diff, top_quantile, col_name='Difference')
        links_short_stack.append(links_short)
    # - concatenate all the links into one
    top_links = pd.concat(links_short_stack, axis=0, ignore_index=True)
    # - output
    top_links.to_csv(Path(TOP_LINKS_DIR) / f'top_links_around_target_genes_{model_name}.csv', index=False)

def identify_top_links(model_name, studies, top_quantile: float) -> None:
    """ Extract top links in the network
    """
    method, DE_type = F_model_name_2_method_and_DE_type(model_name)
    # - get the links for ctr and sample
    links_ctr, links_sample = [pd.read_csv(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv')
                               for study in studies]
    links_ctr, links_sample = [change_protnames_to_genenames_in_links(links)
                               for links in [links_ctr, links_sample]]
    # - determine top links in the network
    links_diff = links_ctr.copy()
    links_diff.drop('Weight', axis=1, inplace=True)
    weight_diff = abs(links_ctr['Weight'] - links_sample['Weight'])
    links_diff['Difference'] = weight_diff
    top_links = choose_top_quantile(links_diff, top_quantile, col_name='Difference')
    # - output
    top_links.to_csv(Path(TOP_LINKS_DIR) / f'top_links_{model_name}.csv', index=False)

if __name__ == '__main__':
    # - parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--top_quantile_target_genes', type=float, default=0.9,
                        help='Top quantile of the top links associated with target genes')
    parser.add_argument('--top_quantile', type=float, default=0.9,
                        help='Top quantile of top links')
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    top_quantile_links_around_target_genes = args.top_quantile_target_genes
    top_quantile = args.top_quantile

    if not os.path.isdir(TOP_LINKS_DIR):
        os.makedirs(TOP_LINKS_DIR)

    for model_name in F_selected_models():
        # - select target links around the target genes (those with the biggest changes from ctr to sample)
        identify_top_links_around_target_genes(model_name, studies, top_quantile=top_quantile_links_around_target_genes)
        # - select target links in the network (those with the biggest changes from ctr to sample)
        identify_top_links(model_name, studies, top_quantile=top_quantile)