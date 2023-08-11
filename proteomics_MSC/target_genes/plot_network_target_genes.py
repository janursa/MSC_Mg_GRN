
import numpy as np
import requests
import sys
import os
from pathlib import Path
import argparse
import pandas as pd
from typing import *


from proteomics_MSC.common_tools.plot_network import PlotCytoscape, LegendPlot, create_nodes_edges
from proteomics_MSC.imports import GRN_DIR, TARGETS_DIR
from proteomics_MSC.common_tools import change_protnames_to_genenames_in_links, F_selected_models, F_model_name_2_method_and_DE_type, F_DE_genenames
from proteomics_MSC.common_tools.links import  choose_top_count, choose_top_quantile

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)


if __name__ == '__main__':
    studies = config['studies']
    # - for each selected model
    for model_name in F_selected_models():
        # - get the links and the gene names
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        gene_names = F_DE_genenames()[DE_type]
        # - get divergence score
        to_load = Path(TARGETS_DIR) / model_name / 'links_divergence_scores.csv'
        divergence_df = pd.read_csv(to_load)
        for study in studies:
            links = pd.read_csv(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv', index_col=False)
            links = change_protnames_to_genenames_in_links(links)
            # merge the links with the divergence scores
            links_merged = pd.merge(links, divergence_df, on=['Regulator','Target'])
            # create nodes edges
            nodes, edges = create_nodes_edges(links_merged, gene_names)
            # - define node color manually
            nodes['color'] = ['#D1B3FF' if (count >= 4) else 'white' for count in
                              nodes['connections']]
            # - send them to cytoscape
            obj = PlotCytoscape(nodes, edges, node_size_attribute='active_sum', edge_color_attribute='divergence_score',
                                edge_size_attribute='weight')
            obj.create_network(name=f'targets_{model_name}_{study}')
            obj.make_changes(style_name=f'style_targets_{model_name}_{study}')





