"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them.
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
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import copy
import igraph as ig
import matplotlib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, F_DE_data, F_DE_protiens, ENRICH_DIR, VSA_DIR
from scripts.utils import calibration, serif_font, make_title_pretty
from scripts.utils.links import retrieve_scores, choose_top_quantile, choose_top_count
from scripts.utils.VSA import RolePlot

def ig_plot(ax, nodes, edges, model_name, study):
    node_names = list(nodes['Protein'].to_numpy())
    F_normalize = lambda vector: vector / np.max(vector)
    vertex_colors = [RolePlot.roles_colors[role] for role in nodes['Role']]
    edge_weights = F_normalize(edges['Weight'].to_numpy(float))
    vertex_sizes = F_normalize(nodes['Weight'].to_numpy(float))

    base_vertex_size = 0.3
    base_edge_weight = 1
    if model_name == 'day1_21_KNN_portia':
        if study == 'ctr':
            base_vertex_size = .15
        if study == 'mg':
            base_vertex_size = .2
    if model_name == 'day1_11_KNN_RF':
        if study == 'ctr':
            base_vertex_size = .15
        if study == 'mg':
            base_vertex_size = .3
    vertex_sizes = vertex_sizes/2 + base_vertex_size

    edge_weights =  4*(edge_weights - min(edge_weights))+  base_edge_weight


    edges_list = [[node_names.index(reg), node_names.index(targ)] for reg, targ in edges.loc[:,['Regulator','Target']].to_numpy(str)]
    # - Construct a graph
    g = ig.Graph(n=len(node_names), edges=edges_list, directed=True)
    arrow_shapes = ["triangle", "vee", "tee", "square", "circle", "diamond"]
    layouts = ["auto", "kamada_kawai","fruchterman_reingold"]
    ig.plot(
        g,
        target=ax,
        layout = layouts[0],
        vertex_size = vertex_sizes,
        vertex_color = vertex_colors,
        vertex_frame_width = .2,
        vertex_frame_color = vertex_colors,
        vertex_label = node_names,
        vertex_label_size = 9,
        edge_width = edge_weights,
        edge_arrow_size=.005  ,
        arrow_shape = arrow_shapes[3],
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


    # grad_colors = ig.GradientPalette("red", "green", len(protnames))
    # - create the layout
    # coords = []
    # for i in range(len(g.vs)):
    #     coords.append([i,i])
    # layout_obj = ig.Layout(coords)

def compare_to_golden_links(links_top, golden_links):
    l_regs = links_top['Regulator'].to_numpy(str)
    l_targs = links_top['Target'].to_numpy(str)

    s_regs = golden_links['Regulator'].to_numpy(str)
    s_targs = golden_links['Target'].to_numpy(str)

    n = 0
    for reg, target in zip(s_regs, s_targs):
        if ((l_regs == reg) & (l_targs == target)).any():
            n += 1

    return n
def manual_filter(links, model_name, study):
    if model_name == 'day1_21_KNN_portia':
        if study == 'ctr':
            to_go = [
                ['Q7Z417', 'P40926'],
                ['Q7Z417', 'Q14011'],
                ['P53621', 'Q7Z417'],
                ['Q7Z417', 'P53621'],
                ['Q7Z417','P13010'],
                ['Q7Z417','P62328'],
                ['Q7Z417','Q9H0U4'],
                ['Q7Z417','Q15233'],
                ['P40926', 'Q7Z417'],
                ['P62328','Q7Z417'],
                ['Q15843', 'Q7Z417'],
                ['Q7Z417','Q15843'],
                ['Q14011', 'Q7Z417'],
                ['Q15233','Q7Z417'],
                ['P13010','Q7Z417']



            ]
        if study == 'mg':
            to_go = []
    if model_name == 'day1_11_KNN_RF':
        if study == 'ctr':
            to_go = [
                    ['O00629', 'Q07866'],
                    #  ['P67936','Q14011'],
                     # ['P67936', 'Q9UKY7'],
                #      ['P67936', 'O00629'],
                #     ['P67936','Q96C90'],
                # ['P67936','P67809'],
                # ['P67936','Q9UN86']

                     ]
        if study == 'mg':
            to_go = [['Q3SX28','P67936']]
    for item in to_go:
        links = links.drop(links[(links['Regulator']==item[0]) & (links['Target']==item[1])].index)

        # links.loc[(links['Target'] == item[0]) & (links['Regulator'] == item[1]), 'Weight'] = 0
    return links
def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data
def determine_target_genes(model_name):
    if model_name == 'day1_11_KNN_RF':
        target_genes = ['P67936', 'Q07866', 'Q9UKY7', 'Q07866']
    if model_name == 'day1_21_KNN_portia':
        target_genes = ['P56192', 'Q13263', 'Q7Z417']
    return target_genes
def filter_target_genes(links, target_genes):
    return links.loc[(links['Target'].isin(target_genes)) | (links['Regulator'].isin(target_genes)), :]
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
    #- dir
    GRN_VISUALIZE_DIR = Path(GRN_DIR)/'visualize'
    if not os.path.isdir(GRN_VISUALIZE_DIR):
        os.makedirs(GRN_VISUALIZE_DIR)
    #- settings
    top_n_links = 10
    figsize = (6,6)
    studies = ['ctr','mg']
    methods = ['RF', 'ridge', 'portia']
    selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    # selected_models = ['day1_11_KNN_RF']

    model_i = 0
    for idx, (DE_type, DE_proteins) in enumerate(F_DE_protiens().items()):
        for method in methods:
            model_name = '_'.join([DE_type,method])
            if model_name not in selected_models: #only selected models
                continue
            target_genes = determine_target_genes(model_name)
            #- plot for each study seperately
            for i_study, study in enumerate(studies):
                #- retireve links
                links = pd.read_csv(Path(GRN_DIR)/method/f'links_{DE_type}_{study}.csv', index_col=False)
                #- filter it
                links = filter_target_genes(links, target_genes) #- choose only those connected with target genes
                links = manual_filter(links, model_name, study) #- remove to make plot look better
                links = choose_top_count(links, top_n_links) #- choose top n links


                #- plot
                edges = links
                nodes = create_nodes(edges, model_name, study)
                fig, ax = plt.subplots(1, 1, figsize=figsize, tight_layout=True)
                ig_plot(ax, nodes, edges, model_name, study)
                fig.savefig(Path(GRN_VISUALIZE_DIR) / f'GRN_{model_name}_{study}.pdf', bbox_inches='tight')
                fig.savefig(Path(GRN_VISUALIZE_DIR) / f'GRN_{model_name}_{study}.png', dpi=300, transparent=True)

            model_i+=1



