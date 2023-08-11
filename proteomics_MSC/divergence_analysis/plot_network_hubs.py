import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path


from proteomics_MSC.imports import MAIN_DIR, DIVERGENCE_ANALYSIS_DIR, DIVERGENCE_ANALYSIS_ROBUSTNESS_DIR, GRN_DIR
from proteomics_MSC.common_tools import F_DE_genenames, F_selected_models, change_protnames_to_genenames_in_links, F_model_name_2_method_and_DE_type
from proteomics_MSC.common_tools.links import choose_top_quantile
from proteomics_MSC.common_tools.plot_network import create_nodes_edges, PlotCytoscape
def assign_node_color(tag, nodes, target_genes):
    """We assign node color customly as it is different from one plot to another"""
    if tag == 'target_genes':
        nodes['color'] = ['cyan' if (node in target_genes) else 'white' for node in nodes['node']]


    return nodes

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    top_quantile = 0.9
    for model_name in F_selected_models():
        # - load divergence score
        to_load = Path(DIVERGENCE_ANALYSIS_ROBUSTNESS_DIR) / model_name / 'links_divergence_scores.csv'
        links_divergence_scores = pd.read_csv(to_load)
        # - select top quantile
        divergence_df = choose_top_quantile(links_divergence_scores, top_quantile,
                                                      col_name='divergence_score')
        gene_names = pd.concat([divergence_df['Regulator'], divergence_df['Target']]).unique()
        # - create nodes and edges
        divergence_df['Weight'] = 1
        nodes, edges = create_nodes_edges(divergence_df, gene_names)
        nodes['active_sum'] = 30
        # - define node color manually
        nodes['color'] = ['#D1B3FF' if (count >= 4) else 'white' for count in
                          nodes['connections']]
        # - send them to cytoscape
        obj = PlotCytoscape(nodes, edges, node_size_attribute='active_sum', edge_color_attribute='divergence_score',
                            edge_size_attribute='weight')
        obj.create_network(name=f'hubs_{model_name}')
        obj.make_changes(style_name=f'style_hub_{model_name}')


    # - plot the legends
    if False:
        def save_fig(fig, name, to_save_dir):
            plt.tight_layout()
            if not os.path.isdir(to_save_dir):
                os.makedirs(to_save_dir)
            fig.savefig(to_save_dir / f'{name}.pdf', bbox_inches='tight')
            fig.savefig(to_save_dir / f'{name}.png', dpi=300, transparent=True)


        # - plot the edge color map
        edge_colormap = LegendPlot.create_edge_color_legend()
        save_fig(edge_colormap, 'edge_colormap', )
        # - plot node size, regulatory importance
        node_size_legend = LegendPlot.create_node_size_legend()
        save_fig(node_size_legend, 'node_size_legend', )
        # - plot edge width legend, weight
        edge_width_legend = LegendPlot.create_edge_weight_legend()
        save_fig(edge_width_legend, 'edge_width_legend', )