"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them.

"""
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import igraph as ig
import json
from typing import List, Dict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import GRN_VISUALIZE_DIR, VSA_DIR, VSA_NOISE_DIR, F_protnames_to_genenames, F_selected_models
from common_tools.links import  choose_top_count, normalize_links
from common_tools.VSA import RolePlot
from VSA.analyse_roles import retreive_links_with_genenames

def manual_adjustments_2_graphs(model_name, study, vertex_sizes, edge_weights):
    #- adjust vertex size (node size)
    if model_name == 'early_MinProb_portia':
        if study == 'ctr':
            scale, base = .4, .4
        if study == 'mg':
            scale, base = .3, .3
    if model_name == 'late_KNN_portia':
        if study == 'ctr':
            scale, base = .3, .35
        if study == 'mg':
            scale, base = .3, .4
    def standardize(vector, scale, base):
        vector = vector - min(vector)
        vector = vector/max(vector)
        return scale*vector + base

    vertex_sizes = standardize(vertex_sizes, scale, base)
    #- adjust edge weight
    base_edge_weight = 1
    edge_weights = 4 * (edge_weights - min(edge_weights)) + base_edge_weight
    return vertex_sizes, edge_weights
def ig_plot(ax, nodes, edges, model_name, study, target_genes):
    """Plot protein connections using igplot"""
    node_names = nodes['Protein'].tolist()
    F_normalize = lambda vector: vector / np.max(vector)
    #- node colors are only for target proteins
    # node_colors = [RolePlot.roles_colors[role] if (gene in target_genes) else 'white' for gene, role in zip(nodes['Protein'], nodes['Role'])]
    target_genes_colors = {gene:color for gene, color in zip(target_genes, colors[0:len(target_genes)])}
    node_colors = [target_genes_colors[gene] if (gene in target_genes) else 'white' for gene in nodes['Protein']]

    edge_weights = F_normalize(edges['Weight'].to_numpy(float))
    vertex_sizes = F_normalize(nodes['Weight'].to_numpy(float))
    #- manual adjustement to the sizes
    vertex_sizes, edge_weights = manual_adjustments_2_graphs(model_name, study, vertex_sizes, edge_weights)
    #- convert gene-gene weight into edges numbers
    edges_list = [[node_names.index(reg), node_names.index(targ)] for reg, targ in edges.loc[:,['Regulator','Target']].to_numpy(str)]
    # - Construct a graph
    g = ig.Graph(n=len(node_names), edges=edges_list, directed=True)
    layouts = ["auto", "kamada_kawai","fruchterman_reingold"]
    ig.plot(
        g,
        target=ax,
        layout = layouts[0],
        vertex_size = vertex_sizes,
        vertex_color = node_colors,
        vertex_frame_width = 1,
        vertex_frame_color = 'grey',
        vertex_label = node_names,
        vertex_label_size = 9,
        edge_width = edge_weights,
        edge_color = 'dimgrey',
        edge_arrow_size=.005,
        edge_curved=True
    )

    # - color legend: roles
    if False:
        title = 'Regulatory role'
        handles = RolePlot.create_role_legends(ax)
        l1=ax.legend(handles=handles, loc='upper center', title=title,
                        bbox_to_anchor=(.2, -.1), prop={'size': 10}, title_fontproperties={'size': 9,'weight':'bold'}, frameon=False
                        )
        ax.add_artist(l1)

    # - vertex size legend
    if False:
        size_title = 'Protein importance'
        n_classes = 4
        sizes = (vertex_sizes-min(vertex_sizes))/(max(vertex_sizes)-min(vertex_sizes)) + .1
        sizes = sizes/max(sizes)
        sizes_classes = np.linspace(min(sizes), max(sizes), n_classes)
        handles = []
        for i, size in enumerate(sizes_classes):
            adj_size = 30 * (size) / np.std(sizes)
            handles.append(ax.scatter([], [], marker='o', label=round(size, 1), color='black',
                                      s=adj_size, alpha=1))
        l2 = ax.legend(loc='upper center', bbox_to_anchor=(.5, -.1), handles=handles, title=size_title, fancybox=False,
                       frameon=False, prop={'size': 10}, title_fontproperties={'size': 9,'weight':'bold'} )
        ax.add_artist(l2)
    # - edge weight legend
    if False:
        size_title = 'Regulatory effect'
        n_classes = 4
        sizes = edge_weights / max(edge_weights)
        sizes_classes = np.linspace(min(sizes), max(sizes), n_classes)
        handles = []
        for i, size in enumerate(sizes_classes):
            adj_size = 2*size / np.std(sizes)
            line, = ax.plot([], [], label=round(size, 1), color='black',
                                      linewidth=adj_size, alpha=1)
            line.set_solid_capstyle('round')
            handles.append(line)
        l3 = ax.legend(loc='upper center', bbox_to_anchor=(.8, -.1), handles=handles, title=size_title, fancybox=False,
                       frameon=False, prop={'size': 10},title_fontproperties={'size': 9,'weight':'bold'})
def manual_filter(links, model_name, study):
    to_go = []
    for item in to_go:
        links = links.drop(links[(links['Regulator']==item[0]) & (links['Target']==item[1])].index)
    return links
def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data
def filter_links_for_target_genes(links, target_genes):
    return links.loc[(links['Target'].isin(target_genes)) | (links['Regulator'].isin(target_genes)), :]
def choose_top_n_links_connected_2_target_genes(links, target_genes, top_n):
    links_target_gene_shortlisted_all = []
    for target_gene in target_genes:
        links_target_gene = links.loc[(links['Target']==target_gene) | (links['Regulator']==target_gene), :]
        links_target_gene_shortlisted = choose_top_count(links_target_gene, top_n)
        links_target_gene_shortlisted_all.append(links_target_gene_shortlisted)
    return pd.concat(links_target_gene_shortlisted_all, axis=0, ignore_index=True)
def create_nodes(edges, model_name, study):
    #- extract gene names
    gene_names = edges.loc[:, ['Regulator', 'Target']].to_numpy(str)
    temp_ = []
    [temp_.extend(names) for names in gene_names]
    gene_names = list(set(temp_))
    #- calculate active sum
    nodes_active_sum = [sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in gene_names]
    # - read vsa results and add it to nodes
    vsa_roles = pd.read_csv(Path(VSA_DIR) / f'vsa_{model_name}_{study}.csv', index_col=False)
    vsa_roles_short = vsa_roles[vsa_roles['Entry'].isin(gene_names)]
    gene_roles = [vsa_roles_short.query(f"Entry == '{prot}'")['Role'].iloc[0] for prot in gene_names]
    # - save nodes and edges to file
    nodes = pd.DataFrame({'Protein': gene_names, 'Weight': nodes_active_sum, 'Role': gene_roles})
    # edges.to_csv(Path(GRN_VISUALIZE_DIR) / f'edges_{model_name}_{study}.csv', index=False)
    # nodes.to_csv(Path(GRN_VISUALIZE_DIR) / f'nodes_{model_name}_{study}.csv', index=False)
    return nodes


if __name__ == '__main__':
    # - parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--top_n_links', type=int, default=2, help='Number of links connected to each target proteins to be visualized')
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    top_n_links = args.top_n_links
    #- dir
    if not os.path.isdir(GRN_VISUALIZE_DIR):
        os.makedirs(GRN_VISUALIZE_DIR)
    #- load target genes
    with open(Path(VSA_NOISE_DIR)/ "target_genes.json", 'r') as file:
        target_genes_all = json.load(file)

    # - load names of selected models
    selected_models = F_selected_models()
    # - read the links and store it for selected models
    links_all: Dict[str, List[pd.DataFrame]] = retreive_links_with_genenames(selected_models,
                                                                             F_protnames_to_genenames(), studies)
    #- standardize links
    links_all = {key: [normalize_links(links) for links in links_stack] for key, links_stack in links_all.items()}
    colors = ['lightgreen', 'lightblue', 'orange', 'yellow', 'lightolive']
    #- for each selected model, draw the network
    for model_name in selected_models:
        target_genes= target_genes_all[model_name]
        #- plot for each study seperately
        for i_study, study in enumerate(studies):
            #- retireve links
            links = links_all[model_name][i_study]
            #- filter it
            # links = filter_links_for_target_genes(links, target_genes) #- choose only those connected with target genes
            # links = manual_filter(links, model_name, study) #- remove to make plot look better
            # links = choose_top_count(links, top_n_links) #- choose top n links
            links = choose_top_n_links_connected_2_target_genes(links, target_genes, top_n_links)

            #- plot
            edges = links
            nodes = create_nodes(edges, model_name, study)
            fig, ax = plt.subplots(1, 1, figsize=(6,6), tight_layout=True)
            ig_plot(ax, nodes, edges, model_name, study, target_genes)
            fig.savefig(Path(GRN_VISUALIZE_DIR) / f'GRN_{model_name}_{study}.pdf', bbox_inches='tight')
            fig.savefig(Path(GRN_VISUALIZE_DIR) / f'GRN_{model_name}_{study}.png', dpi=300, transparent=True)




"""
igraph params:

layout: the layout to use for the graph. This can be a precomputed layout, or a string specifying one of several built-in layout algorithms (default is 'auto')
vertex_color: the color to use for vertices (default is 'white')
vertex_size: the size of the vertices in pixels (default is 10)
vertex_shape: the shape of the vertices (default is 'circle')
vertex_frame_color: the color of the border around each vertex (default is 'black')
vertex_label: the label to display for each vertex (default is None)
vertex_label_color: the color of the vertex labels (default is 'black')
vertex_label_size: the size of the vertex labels in points (default is 12)
edge_color: the color to use for edges (default is 'black')
edge_width: the width of the edges in pixels (default is 1)
edge_arrow_size: the size of the arrowheads on directed edges (default is 1)
edge_arrow_width: the width of the arrowheads on directed edges (default is 1)
edge_curved: whether to draw curved edges instead of straight lines (default is False)
edge_label: the label to display for each edge (default is None)
edge_label_color: the color of the edge labels (default is 'black')
edge_label_size: the size of the edge labels in points (default is 12)
margin: the size of the margin around the plot in pixels (default is 10)
background: the background color of the plot (default is 'white')
vertex_label_dist: the distance between the vertex label and the vertex in pixels (default is 0)
vertex_label_angle: the angle of rotation for the vertex label in degrees (default is 0)
vertex_label_family: the font family for the vertex label (default is 'sans-serif')
vertex_label_font: the font style for the vertex label (default is 'normal')
vertex_label_rect: whether to draw a rectangle behind the vertex label (default is False)
edge_label_dist: the distance between the edge label and the edge in pixels (default is 0)
edge_label_angle: the angle of rotation for the edge label in degrees (default is 0)
edge_label_family: the font family for the edge label (default is 'sans-serif')
edge_label_font: the font style for the edge label (default is 'normal')
edge_label_rect: whether to draw a rectangle behind the edge label (default is False)
bbox: a tuple specifying the size of the plot in pixels (default is None)
"""