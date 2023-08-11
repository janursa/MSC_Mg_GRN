import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

from proteomics_MSC.imports import DIVERGENCE_ANALYSIS_DIR, GRN_DIR
from proteomics_MSC.common_tools import F_model_name_2_method_and_DE_type, F_selected_models, change_protnames_to_genenames_in_links
from proteomics_MSC.common_tools.links import choose_top_quantile
from proteomics_MSC.uncertainity_analysis.ua_aux import NetworkNoiseAnalysis


import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)
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
    to_save = Path(DIVERGENCE_ANALYSIS_DIR) / f'top_links_{model_name}.csv'
    top_links.to_csv(to_save, index=False)
    print(f'output -> {to_save}')

if __name__ == '__main__':
    studies = config['studies']
    top_quantile_links_around_target_genes = config['top_quantile_links_change']
    top_quantiles = config['top_quantile_links_change']

    if not os.path.isdir(DIVERGENCE_ANALYSIS_DIR):
        os.makedirs(DIVERGENCE_ANALYSIS_DIR)

    for model_name in F_selected_models():
        # - select target links around the target genes (those with the biggest changes from ctr to sample)
        # identify_top_links_around_target_genes(model_name, studies, top_quantile=top_quantile_links_around_target_genes)
        # - select target links in the network (those with the biggest changes from ctr to sample)
        identify_top_links(model_name, studies, top_quantile=top_quantiles)