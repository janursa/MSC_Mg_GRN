"""
    Sets of functions useful for ploting enriched terms
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
import utils
from ._imports import *
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
def plot(datas, tags, size_tag, color_tag, xlabel, marker_types, figsize,legend_color=True, 
        legend_size=True, legend_marker=True, title=''):
    #-----------define prop---------------
    utils.comic_font()

    scale_scatter_size = 60
    
    #-----------calculations---------------
    x = [j for sub in [data[xlabel].values.tolist() for data in datas] for j in sub]
    terms =[j for sub in [data['Description'].values.tolist() for data in datas] for j in sub]

    y = list(range(len(terms), 0, -1)) # numbers for each term
    sizes = np.array([j for sub in [data[size_tag].values.tolist() for data in datas] for j in sub])
    colors = [j for sub in [data[color_tag].values.tolist() for data in datas] for j in sub]
    counts = [len(data['Description'].values.tolist()) for data in datas]
    markers = [marker_types[i] for sub in [np.full(shape=n, fill_value=i, dtype=int) for i,n in enumerate(counts)] for i in sub]

    #-----------plot-----------------------
    fig, ax = plt.subplots(1,1, figsize=figsize, tight_layout=True)
    sc = mscatter(x, y, c=colors, s=scale_scatter_size*sizes,
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
        for i, marker in enumerate(marker_types):
            handles.append(ax.scatter([],[],marker=marker, label=tags[i], color='black'))
        l2 = plt.legend(handles=handles, 
            bbox_to_anchor=(2.6,1), 
            # title='Enriched Term'
            )
        ax.add_artist(l2)
    #- size legend 
    if legend_size:
        handles = []
        labels = []
        n_classes = 4
        sizes_classes = np.arange(min(sizes), max(sizes), (max(sizes)-min(sizes))/n_classes)
        for i, size in enumerate(sizes_classes):
            handles.append(ax.scatter([],[], marker='o', label=round(size,1),
                color='black', 
                s= size*scale_scatter_size,
                alpha=1
                ))
        l1 = plt.legend(bbox_to_anchor=(1, 1), handles=handles, 
                        title=size_tag, fancybox=False, frameon=False)

        ax.add_artist(l1)
    #- color legend
    if legend_color:
        PCM=ax.get_children()[0]
        CB = plt.colorbar(PCM, ax=ax, 
                    aspect=3, 
                    location='right',
                    anchor=(1, 0), 
                    # extend='both', 
                    ticks=[min(colors), max(colors)], 
                    # format="%4.2e",
                    shrink=1
                    )
        CB.ax.set_title(color_tag)
        # CB.ax.tick_params(labelsize=fontsize['fontsize'])
    return fig
