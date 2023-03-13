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
    edges_scale_factor = 4
    if model_name == 'day1_21_KNN_portia':
        if study == 'ctr':
            base_vertex_size = .15
        if study == 'mg':
            base_vertex_size = .2
    if model_name == 'day1_11_KNN_RF':
        if study == 'ctr':
            base_vertex_size = .15
            edges_scale_factor = 2
        if study == 'mg':
            base_vertex_size = .3
            edges_scale_factor = 3
    vertex_sizes = vertex_sizes/2 + base_vertex_size

    edge_weights =  edges_scale_factor*(edge_weights - min(edge_weights))+  base_edge_weight


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
        edge_arrow_size=.01  ,
        arrow_shape = arrow_shapes[3],
        edge_curved=True
    )

    # - color legend: roles
    if True:
        title = 'Regulatory role'
        handles = RolePlot.create_role_legends(ax)
        l1=ax.legend(handles=handles, loc='upper center', title=title,
                        bbox_to_anchor=(.2, -.1), prop={'size': 10}, title_fontproperties={'size': 9,'weight':'bold'}, frameon=False
                        )
        ax.add_artist(l1)

    # - vertex size legend
    if True:
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
    if True:
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
            to_go = [['Q14767','Q7Z417']]
        if study == 'mg':
            to_go = []
    if model_name == 'day1_11_KNN_RF':
        if study == 'ctr':
            # to_go = [['Q14011', 'P31949'],['Q14011','P14174'],['P31949','P67936'], ['Q8IVL6','P14174'],['P62328','P31949']]

            to_go = [['Q99613', 'P13667']]

        if study == 'mg':
            to_go = []
    for item in to_go:
        links.loc[(links['Regulator']==item[0]) & (links['Target']==item[1]), 'Weight'] = 0
        links.loc[(links['Target'] == item[0]) & (links['Regulator'] == item[1]), 'Weight'] = 0
if __name__ == '__main__':
    #- dir
    GRN_VISUALIZE_DIR = Path(GRN_DIR)/'visualize'
    if not os.path.isdir(GRN_VISUALIZE_DIR):
        os.makedirs(GRN_VISUALIZE_DIR)
    #- settings
    top_n_links_across_groups = 3
    figsize = (7,7)
    studies = ['mg']
    methods = ['RF', 'ridge', 'portia']
    selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    # selected_models = ['day1_11_KNN_RF']
    target_proteins = [['Q02790', 'Q07866', 'Q99613', 'P14174'], ['P02652','P00568','Q02218']]

    model_i = 0
    for idx, (DE_type, DE_proteins) in enumerate(F_DE_protiens().items()):
        # protnames = F_DE_protiens()[DE_type]
        for method in methods:
            model_name = '_'.join([DE_type,method])
            if model_name not in selected_models: #only selected models
                continue

            target_proteins_model = target_proteins[model_i]

            for i_study, study in enumerate(['ctr', 'mg']):

                #- retireve links
                links = pd.read_csv(Path(GRN_DIR)/method/f'links_{DE_type}_{study}.csv', index_col=False)
                #- filter links by choosing only target protein's edges
                links = links.loc[(links['Target'].isin(target_proteins_model)) | (links['Regulator'].isin(target_proteins_model)), :]
                manual_filter(links, model_name, study) #- to keep the network clean
                #- choose top n links for each target protein
                links_stack = []
                for i, prot in enumerate(target_proteins_model):
                    # links.loc[(links['Target'] == prot) | (links['Regulator'] == prot), 'class'] = i
                    class_links = links.loc[(links['Target'] == prot) | (links['Regulator'] == prot), :]
                    links_stack.append(choose_top_count(class_links, top_n_links_across_groups))

                edges = pd.concat(links_stack)
                edges.drop_duplicates(inplace=True)

                #- shotlisted node names and AS
                nodes_names = edges.loc[:,['Regulator', 'Target']].to_numpy(str)
                temp_ = []
                [temp_.extend(names) for names in nodes_names]
                nodes_names = list(set(temp_))
                nodes_active_sum= [sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in nodes_names]
                #- read vsa results and add it to nodes
                vsa_roles = pd.read_csv(Path(VSA_DIR) / f'vsa_{model_name}_{study}.csv', index_col=False)
                vsa_roles_short = vsa_roles[vsa_roles['Entry'].isin(nodes_names)]
                nodes_roles = [vsa_roles_short.query(f"Entry == '{prot}'")['Role'].iloc[0] for prot in nodes_names]
                #- save nodes and edges to file
                nodes = pd.DataFrame({'Protein': nodes_names, 'Weight': nodes_active_sum, 'Role': nodes_roles})
                edges.to_csv(Path(GRN_VISUALIZE_DIR) / f'edges_{model_name}_{study}.csv', index=False)
                nodes.to_csv(Path(GRN_VISUALIZE_DIR) / f'nodes_{model_name}_{study}.csv', index=False)
                #- plot
                fig, ax = plt.subplots(1, 1, figsize=figsize, tight_layout=True)
                ig_plot(ax, nodes, edges, model_name, study)
                fig.savefig(Path(GRN_VISUALIZE_DIR) / f'GRN_{model_name}_{study}.pdf', bbox_inches='tight')
                fig.savefig(Path(GRN_VISUALIZE_DIR) / f'GRN_{model_name}_{study}.png', dpi=300, transparent=True)

                # extra_artists_stack.extend(extra_artists)

            # if model_i==0:
            #     plt.show()
            # plt.tight_layout(pad=2, w_pad=5, h_pad=3)
            # fig.savefig(Path(GRN_VISUALIZE_DIR)/ f'GRN_{model_name}.pdf',bbox_extra_artists=(ll,), bbox_inches='tight')

            model_i+=1

    # golden_links = pd.read_csv(Path(ENRICH_DIR)/f'network_{DE_type}.csv', index_col=False)


