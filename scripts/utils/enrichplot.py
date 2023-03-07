"""
    Sets of functions useful for ploting enriched terms
"""
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from scripts.imports import OUTPUT_DIR
from scripts.utils import comic_font, serif_font

def mscatter(x,y,ax=None, m=None, **kw):
    """ a custom plot built on scatter plot to enable multi marker visualization
    """
    import matplotlib.markers as mmarkers
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc
def plot_enrich(df_stack, tags, size_tag, color_tag, xlabel, marker_types, figsize,legend_color=True,
        legend_size=True, legend_marker=True, title='', scale_factor=.1):
    #-----------define prop---------------
    import matplotlib
    # comic_font()
    serif_font()
    matplotlib.rcParams.update({'font.size': 14})
    #-----------assignment---------------
    def flatten(data, tag):
        return [j for sub in [item[tag].values.tolist() for item in data] for j in sub]
    x = flatten(df_stack, xlabel)
    terms = flatten(df_stack,'Description' )
    y = list(range(len(terms), 0, -1)) # numbers for each term
    sizes = np.array(flatten(df_stack, size_tag))
    colors = flatten(df_stack, color_tag)
    counts = [len(data['Description'].values) for data in df_stack]
    markers = [marker_types[i] for sub in [np.full(shape=n, fill_value=i, dtype=int) for i,n in enumerate(counts)] for i in sub]
    #-----------scale----------------------
    adj_sizes = scale_factor * len(y) * figsize[1] * sizes/np.std(sizes)

    #-----------plot-----------------------
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    sc = mscatter(x, y, c=colors, s=adj_sizes,
        m=markers, ax=ax, 
        cmap='autumn'
        # cmap='Spectral'
        )
    ax.set_yticks(y)
    ax.set_yticklabels(terms)
    ax.set_xticks([min(x), int((max(x)+min(x))/2), max(x)])
    ax.set_xlabel('Protein count')
    ax.set_xmargin(0.2)
    ax.set_title(title)
    #- marker legend
    if legend_marker:
        handles = []
        for i, tag in enumerate(tags):
            handles.append(ax.scatter([],[],marker=marker_types[i], label=tag, color='black', s=100))
        l2 = plt.legend(loc='right', bbox_to_anchor=(-.5,-.3), handles=handles, title='Enriched Term')
        ax.add_artist(l2)
    #- size legend 
    if legend_size:
        handles = []
        n_classes = 4
        sizes_classes = np.linspace(min(sizes), max(sizes), n_classes)

        for i, size in enumerate(sizes_classes):
            adj_size = 50*size/np.std(sizes)
            handles.append(ax.scatter([],[], marker='o', label=round(size,1), color = 'black',
                                      s = adj_size,alpha = 1))
        l1 = ax.legend(loc='right', bbox_to_anchor=(.2,-.3), handles=handles, title=size_tag, fancybox=False, frameon=False)

        ax.add_artist(l1)
    #- color legend: FDR
    if legend_color:
        PCM=ax.get_children()[0]
        CB = plt.colorbar(PCM, ax=ax, 
                    aspect=3, 
                    location='right',
                    anchor=(0, 0),
                    # extend='both', 
                    ticks=[min(colors), max(colors)], 
                    # format="%4.2e",
                    shrink=.3
                    )
        CB.ax.set_title(color_tag)
        # CB.ax.tick_params(labelsize=fontsize['fontsize'])
    return fig
