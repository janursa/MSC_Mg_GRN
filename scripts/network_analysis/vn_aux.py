import requests
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
import igraph as ig
import json
from typing import List, Dict
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import NETWORK_ANALYSIS_DIR, VSA_DIR
from common_tools import flatten
from uncertainity_analysis.ua_aux import NoiseAnalysis
from common_tools.role_analysis import RolePlot

# target_node_colors = ["#B0C4DE", "#C1FFC1", "#D8BFD8", "#87CEFA", "#98FB98", "#DDA0DD"]
target_node_colors = RolePlot.roles_colors

class PlotCytoscape:
    edge_color_map = mcolors.LinearSegmentedColormap.from_list("black_to_red", ["yellow", "red"])
    def __init__(self, nodes, edges, node_size_attribute=None, edge_color_attribute=None, edge_size_attribute=None):
        assert (len(nodes) > 0)
        assert (len(edges) > 0)
        # - scale node size
        if node_size_attribute is not None:
            scale, min_value = 20, 30
            nodes["size"] = self.normalize(nodes[node_size_attribute], scale=scale, min_value=min_value)
        # - scale edge width
        if edge_size_attribute is not None:
            scale, min_value = 3, 1
            edges["width"] = self.normalize(edges[edge_size_attribute], scale=scale,
                                                         min_value=min_value)
        # - assign color to nodes
        # if node_color_attribute is not None:
        #     nodes['color'] = self.assign_node_colors(nodes[node_color_attribute])
        # - assign color to edges
        if edge_color_attribute is not None:
            edges['edge_color'] = self.assign_edge_colors(edges[edge_color_attribute])

        self.base = 'http://localhost:' + str(1234) + '/v1/'
        self.ID = None
        self.nodes = self.extract_dict(nodes.rename(columns={'node': 'id'}))
        self.edges = self.extract_dict(edges)
    # @staticmethod
    # def assign_node_colors(values):
    #
    #     # colors = [target_node_colors[color_code] if (color_code != 0) else 'white' for color_code in values]
    #     colors = [target_node_colors[color_code] for color_code in values]
    #     return colors

    @classmethod
    def assign_edge_colors(cls, values):
        # - scale edge color
        scale, min_value = 1, 0  # fixed for color mapping
        values = cls.normalize(values, scale=scale,
                                        min_value=min_value)
        colors = [cls.edge_color_map(value) for value in values]
        colors = ['#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255)) for r, g, b, _ in colors]
        return colors
    @staticmethod
    def normalize(values, scale=1, min_value=0):
        min_val = min(values)
        max_val = max(values)
        if min_val == max_val:
            return values
        normalized_values = [min_value+ (val - min_val) / (max_val - min_val)*scale for val in values]
        return normalized_values
    @staticmethod
    def extract_dict(df) -> List[Dict]:
        """
        Puts the df into a list of dict where it presents each row of the df
        """
        nodes_dict = df.to_dict('records')
        list_of_data = []
        for node in nodes_dict:
            data = {'data': {}}
            for col in df.columns:
                data['data'][col] = node[col]
            list_of_data.append(data)
        return list_of_data
    def _delete_network(self, name):
        """
            If the name exist, delete the network before creating it
        """
        response = requests.get(self.base + 'networks.names')
        networks = response.json()
        network_suid = None
        for network in networks:
            if network['name'] == name:
                network_suid = network['SUID']
                break

        # Delete the network if it exists
        if network_suid is not None:
            response = requests.delete(self.base + f'networks/{network_suid}')
            if response.status_code == 200:
                print(f"Network '{name}' deleted successfully")
            else:
                print(f"Failed to delete network '{name}'")
                print(response.status_code)
                print(response.text)
        else:
            print(f"No network found with name '{name}'")
    def create_network(self, name='mynetwork'):
        self._delete_network(name)
        HEADERS = {'Content-Type': 'application/json'}
        elements = {"nodes": self.nodes, "edges": self.edges}
        network = {
            'data': {
                'name': name
            },
            'elements': elements,
        }

        response = requests.post(self.base + 'networks?collection=My%20Collection', data=json.dumps(network),
                                 headers=HEADERS)
        print(response.text)
        res_dict = response.json()
        new_suid = res_dict['networkSUID']
        self.ID = new_suid
        print(f"Network '{name}' created successfully")
    @staticmethod
    def check(response, code, name):
        if response.status_code == code:
            print(f"{name} successful")
        else:
            print(f"Failed to {name}")
            print(response.status_code)
            print(response.text)
    @staticmethod
    def passthrough_map(style_data, attribute, attribute_type, visual_property):
        """
            Used to map the given attribute to the visual property. The column attribute should be present.
        """
        mapping = {
            "mappingType": "passthrough",
            "mappingColumn": attribute,
            "mappingColumnType": attribute_type,
            "visualProperty": visual_property
        }
        style_data['mappings'].append(mapping)
    @staticmethod
    def single_value_map(style_data, visual_property, value):
        """
            Used to set a single value to the given visual property
        """
        mapping = {
            "visualProperty": visual_property,
            "value" : value
        }
        style_data['defaults'].append(mapping)

    def make_changes(self, style_name="CustomStyle"):
        # - apply layout
        response = requests.get(self.base + f'apply/layouts/force-directed/{self.ID}')
        self.check(response, 200, 'directed force')
        # - get the Directed style
        response = requests.get(self.base + 'styles/Directed')
        self.check(response, 200, 'directed')
        directed_style = response.json()
        # - Create a new style
        custom_style = directed_style
        custom_style['title'] = style_name
        # - make changes
        self.passthrough_map(custom_style, "color", 'String', 'NODE_FILL_COLOR')
        self.passthrough_map(custom_style, "size", 'Float', 'NODE_SIZE')
        self.passthrough_map(custom_style, "width", 'Float', 'EDGE_WIDTH')
        self.passthrough_map(custom_style, "edge_color", 'String', 'EDGE_STROKE_UNSELECTED_PAINT')
        self.passthrough_map(custom_style, "edge_color", 'String', 'EDGE_TARGET_ARROW_UNSELECTED_PAINT')
        self.passthrough_map(custom_style, "edge_color", 'String', 'EDGE_TARGET_ARROW_UNSELECTED_PAINT')
        self.single_value_map(custom_style, "NODE_LABEL_COLOR", "#000000")
        self.single_value_map(custom_style, "NODE_BORDER_WIDTH", 1)
        self.single_value_map(custom_style, "NODE_LABEL_FONT_SIZE", 8)
        # - delete the previous style with the same name
        response = requests.delete(self.base + f'styles/{style_name}')
        self.check(response, 200, 'delete styles')
        # - add style
        response = requests.post(self.base + 'styles', json=custom_style)
        self.check(response, 201, 'custom style added')
        # - apply style
        response = requests.get(self.base + f'apply/styles/{style_name}/{self.ID}')
        self.check(response, 200, 'apply style')

class PlotIgplot:
    @classmethod
    def plot(cls, nodes: pd.DataFrame, edges: pd.DataFrame, model_name: str, study: str) -> None:
        """
        Plots gene network using igplot
        """
        # - manual adjustement to the sizes
        norm_func = lambda vector: vector / np.max(vector)
        vertex_sizes = norm_func(nodes['active_sum'].to_numpy(float))
        edge_weights = norm_func(edges['weight'])
        # - convert gene-gene weight into edges numbers
        node_names = nodes['node'].values.tolist()
        edges_list = [[node_names.index(reg), node_names.index(targ)] for reg, targ in
                      edges.loc[:, ['source', 'target']].to_numpy(str)]
        # - Construct a graph
        g = ig.Graph(n=len(node_names), edges=edges_list, directed=True)
        layouts = ["auto", "kamada_kawai", "fruchterman_reingold"]

        # - create a color map for edge color
        colormap = mcolors.LinearSegmentedColormap.from_list(
            "black_to_red", ["black", "red"]
        )
        norm = mcolors.Normalize(vmin=min(edges['divergence_score']), vmax=max(edges['divergence_score']))
        edge_colors = [colormap(norm(value)) for value in edges['divergence_score']]
        # - nodes colors only for target
        node_colors = nodes['color']
        # - create the plot
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        ig.plot(
            g,
            target=ax,
            layout=layouts[0],
            vertex_size=vertex_sizes,
            vertex_color=node_colors,
            vertex_frame_width=1,
            vertex_frame_color='grey',
            vertex_label=node_names,
            vertex_label_size=9,
            edge_width=edge_weights,
            edge_color=edge_colors,
            edge_arrow_size=.005,
            edge_curved=True
        )
        # - add the colormap to the graph
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        cax = fig.add_axes([0.85, 0.3, 0.025, 0.4])
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])  # This is required to make the colorbar show up
        fig.colorbar(sm, cax=cax, label="Significance", aspect=20)
        fig.subplots_adjust(right=0.8)

        return fig





class LegendPlot():
    title_size = 9
    size = 10
    @staticmethod
    def __plot_size(tag):
        if tag == 'vertical':
            return (2, 2)
    @staticmethod
    def define_size_classes(ax):
        vertex_sizes = np.linspace(0, 1, 100)
        n_classes = 4
        sizes = (vertex_sizes - min(vertex_sizes)) / (max(vertex_sizes) - min(vertex_sizes)) + .1
        sizes = sizes / max(sizes)
        sizes_classes = np.linspace(min(sizes), max(sizes), n_classes)
        return sizes, sizes_classes
    @staticmethod
    def __save_fig(fig, name):
        plt.tight_layout()
        to_save_dir = Path(NETWORK_ANALYSIS_DIR) / 'legends'
        if not os.path.isdir(to_save_dir):
            os.makedirs(to_save_dir)
        fig.savefig(to_save_dir / f'{name}.pdf', bbox_inches='tight')
        fig.savefig(to_save_dir / f'{name}.png', dpi=300, transparent=True)
    @classmethod
    def create_legends(cls, ax, handles, title):
        ax.legend(loc='upper center', bbox_to_anchor=(.5, .7), handles=handles, title=title,
                  fancybox=False,
                  frameon=False, prop={'size': cls.size}, title_fontproperties={'size': cls.title_size, 'weight': 'bold'})
    @classmethod
    def create_node_size_legend(cls):
        title = 'Active Sum'
        fig, ax = plt.subplots(figsize=cls.__plot_size('vertical'))
        sizes, sizes_classes = cls.define_size_classes(ax)
        handles = []
        for i, size in enumerate(sizes_classes):
            adj_size = 30 * (size) / np.std(sizes)
            handles.append(ax.scatter([], [], marker='o', label=round(size, 1), color='black',
                                      s=adj_size, alpha=1))

        cls.create_legends(ax, handles, title)
        cls.posprocess(ax)
        cls.__save_fig(fig, 'node_size_legend')
    @staticmethod
    def posprocess(ax):
        ax.set_axis_off()
    @classmethod
    def create_edge_weight_legend(cls):
        title = 'Regulatory effect'
        fig, ax = plt.subplots(figsize=cls.__plot_size('vertical'))
        sizes, sizes_classes = cls.define_size_classes(ax)
        handles = []
        for i, size in enumerate(sizes_classes):
            adj_size = 2 * size / np.std(sizes)
            line, = ax.plot([], [], label=round(size, 1), color='black',
                            linewidth=adj_size, alpha=1)
            line.set_solid_capstyle('round')
            handles.append(line)

        cls.create_legends(ax, handles, title)
        cls.posprocess(ax)
        cls.__save_fig(fig, 'edge_width_legend')
    @classmethod
    def create_edge_color_legend(cls):
        title = ''
        sample_array = np.linspace(0, 1, 100).reshape(10, 10)
        fig, ax = plt.subplots(figsize=(1.1, 1.4))
        cmap = PlotCytoscape.edge_color_map
        cbar = plt.colorbar(plt.imshow(sample_array, cmap=cmap, aspect='auto'),
                            cax=ax, orientation='vertical', ticks=np.round(np.linspace(.1, 1, 4), decimals=1))
        cbar.ax.tick_params(labelsize=13)
        cls.__save_fig(fig, 'edge_colormap')
