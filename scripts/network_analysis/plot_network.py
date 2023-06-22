
import numpy as np
import requests
import sys
import os
from pathlib import Path
import argparse
import pandas as pd
from typing import *

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from network_analysis.vn_aux import PlotCytoscape, PlotIgplot, LegendPlot
from imports import change_protnames_to_genenames_in_links, NETWORK_ANALYSIS_DIR, TOP_LINKS_WITH_DISTANCE_SCORES_DIR, F_selected_models, F_model_name_2_method_and_DE_type, F_DE_genenames, GRN_DIR
from uncertainity_analysis.ua_aux import NoiseAnalysis
from common_tools.links import  choose_top_count, choose_top_quantile
from common_tools.role_analysis import RolePlot

# def delete_single_connections(nodes, edges):
#     """delete those nodes with only one connection"""
#     removed_nodes = nodes[nodes['connections'] < 2]
#     removed_nodes = removed_nodes['node'].tolist()
#     # - filter nodes
#     filtered_nodes = nodes[nodes['connections'] >= 2]
#     filtered_edges = edges[
#         (~edges['source'].isin(removed_nodes)) &
#         (~edges['target'].isin(removed_nodes))
#         ]
#     return filtered_nodes, filtered_edges
def manual_filter(nodes, edges, model_name):
    genes_to_go = []
    # if model_name == 'late_KNN_portia':
    #     genes_to_go = ['UCHL3', 'NEDD8', 'AP2M1', 'LRP1', 'CMPK1', 'SKP1', 'OGDH', 'ATP6V1B2', 'PXN', 'COPA', 'YWHAG']
    # else:
    #     genes_to_go = ['NUFIP2', 'KRT8', 'PFDN4']
    nodes = nodes.loc[~nodes['node'].isin(genes_to_go)].reset_index()
    edges = edges.loc[~edges['source'].isin(genes_to_go)].reset_index()
    edges = edges.loc[~edges['target'].isin(genes_to_go)].reset_index()
    # print()
    return nodes, edges
def add_connection_counts(nodes, edges):
    """Counts number of times a node is repeated in the network"""
    source_counts = edges.groupby("source").size()
    target_counts = edges.groupby("target").size()

    # Combine the counts from both columns
    combined_counts = pd.concat([source_counts, target_counts], axis=1).fillna(0)
    combined_counts.columns = ["source_count", "target_count"]

    # Calculate the total number of connections for each gene
    combined_counts["connections"] = combined_counts["source_count"] + combined_counts["target_count"]

    # Reset the index so that the gene names become a regular column
    combined_counts.reset_index(inplace=True)
    combined_counts.rename(columns={"index": "node"}, inplace=True)

    # Merge the connections column with the nodes DataFrame
    nodes = nodes.merge(combined_counts[["node", "connections"]], on="node")
    return nodes, edges
def assign_node_color(tag, nodes, target_genes):
    """We assign node color customly as it is different from one plot to another"""
    if tag == 'target_genes':
        nodes['color'] = ['cyan' if (node in target_genes) else 'white' for node in nodes['node']]
    else:
        nodes['color'] = ['#D1B3FF' if (count >= 4) else 'white' for count in
                          nodes['connections']]
    return nodes

def plot_network(model_name, studies, top_links_n: int | float, tag: str = 'top_links', target_genes=None):
    """
        Plot the network
    """
    assert (tag in ['top_links', 'top_links_same_weight', 'target_genes'])
    divergence_scores = NoiseAnalysis.retrieve_divergence_scores_df(model_name)
    if tag == 'target_genes':
        divergence_scores_short = filter_for_around_target_genes(divergence_scores, target_genes, top_links_n)
    else:
        divergence_scores_short = choose_top_quantile(divergence_scores, top_links_n,
                                                      col_name='divergence_score')
    # - plot for each study
    for i_study, study in enumerate(studies):
        # - create nodes and edges
        nodes, edges = create_nodes_edges(model_name, study, divergence_scores_short)
        if tag == 'top_links_same_weight':
            edges['weight'] = 1
            nodes['active_sum'] = 30
        # - manual filter
        nodes, edges = manual_filter(nodes, edges, model_name)
        # - define node color manually
        nodes = assign_node_color(tag, nodes, target_genes)
        # - send them to cytoscape
        if plot_cytoscape:
            obj = PlotCytoscape(nodes, edges, node_size_attribute='active_sum', edge_color_attribute='divergence_score', edge_size_attribute='weight')
            obj.create_network(name=f'{tag}_{model_name}_{study}')
            obj.make_changes(style_name=f'style_{tag}_{model_name}_{study}')
        # - plot igraph
        if plot_igraph:
            fig = PlotIgplot().plot(nodes, edges, model_name, study)
            to_save_dir = Path(NETWORK_ANALYSIS_DIR) / tag / 'igraphs'
            if not os.path.isdir(to_save_dir):
                os.makedirs(to_save_dir)
            fig.savefig(to_save_dir / f'{tag}_{model_name}_{study}.pdf', bbox_inches='tight')
def filter_for_around_target_genes(links: pd.DataFrame, target_genes: List[str],
                                         top_n_links: int) -> pd.DataFrame:
    """
        Filters the links around target genes to only a few
    """
    # - top n connections around the target genes to show
    links_to_show = []
    for gene in target_genes:
        links_gene = links.loc[(links['Regulator'] == gene) | (links['Target'] == gene),
                          :].reset_index(drop=True)
        links_gene = choose_top_quantile(links_gene, top_n_links,
                                                 col_name='divergence_score')
        links_to_show.append(links_gene)
    links_to_show = pd.concat(links_to_show).reset_index(drop=True)
    links_to_show.drop_duplicates(inplace=True)
    return links_to_show

def create_nodes_edges_basics(links, gene_names, map_names):
    """
        Format links into nodes and edges. Basic attributes.
    """
    nodes_active_sum = [sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in gene_names]
    nodes = pd.DataFrame({'node': gene_names, 'active_sum': nodes_active_sum})

    links.rename(columns=map_names, inplace=True)
    return nodes, links
def create_nodes_edges(model_name, study, divergence_scores_short) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """ Create nodes and edges

    """
    method, DE_type = F_model_name_2_method_and_DE_type(model_name)
    gene_names = F_DE_genenames()[DE_type]
    # - get links
    links = pd.read_csv(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv', index_col=False)
    # - change the names from protnames to genenames
    links = change_protnames_to_genenames_in_links(links)
    # - select only
    links_short = links.merge(divergence_scores_short, on=['Regulator', 'Target'])
    # - create the ExtractNodesEdges obj
    map_names = {'Regulator':'source', 'Target':'target', 'Weight': 'weight'}
    nodes, edges = create_nodes_edges_basics(links_short, gene_names, map_names)
    # - determine the number of interactions
    nodes, edges = add_connection_counts(nodes, edges)

    return nodes, edges
if __name__ == '__main__':
    # - parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--igraph', type=bool, default=False, help="Plot using Igraph")
    parser.add_argument('--cytoscape', type=bool, default=True, help="Plot using Cytoscape")
    parser.add_argument('--legends', type=bool, default=True, help="Plot legends")
    parser.add_argument('--top_n_links_around_target_genes', type=int, default=2,
                        help='Number of links connected to each target proteins to be visualized')
    parser.add_argument('--top_links_n', type=float, default=0.9,
                        help='Top percentile of the significant links to be visualized')
    args, remaining_args = parser.parse_known_args()
    studies = args.studies
    plot_igraph = args.igraph
    plot_cytoscape = args.cytoscape
    plot_legends = args.cytoscape
    top_n_links_around_target_genes = args.top_n_links_around_target_genes
    top_links_n = args.top_links_n
    # - target genes
    target_genes_models = {'early_MinProb_portia': ['ACO1', 'HDGFL2', 'MYL1'],
                    'late_KNN_portia':['MDH2', 'TRIM28', 'MYL1', 'GLS']}
    # - for each selected model
    for model_name in F_selected_models():
        if False:
            plot_network(model_name, studies, top_links_n, "top_links")
        if True:
            plot_network(model_name, studies, top_links_n, "target_genes", target_genes_models[model_name])
        if False: # - to identify hops, we dont need to plot ctr vs mg
            plot_network(model_name, studies, top_links_n, "top_links_same_weight")

    # - plot the legends
    if plot_legends:
        # - plot the edge color map
        LegendPlot.create_edge_color_legend()
        # - plot node size, regulatory importance
        LegendPlot.create_node_size_legend()
        # - plot edge width legend, weight
        LegendPlot.create_edge_weight_legend()
