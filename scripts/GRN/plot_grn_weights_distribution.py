"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, F_DE_protnames
from scripts.utils import  serif_font

def plot_weight(links, color, name, dist_key='WeightPool'):
    serif_font()
    links = links.sort_values('Weight', ascending=False).reset_index(drop=True)

    top_n = 5
    nrows = 1
    ncols = 5
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(ncols * 2.2, nrows * 2))
    for idx in range(top_n):
        # i = int(idx / (nrows - 1))
        # j = idx % ncols
        # ax = axes[i][j]
        ax = axes[idx]
        ax.hist(links[dist_key][idx],
                bins=30,
                alpha=1,
                histtype='stepfilled',  # 'bar', 'barstacked', 'step', 'stepfilled'
                color=color,
                ec='lightgreen',
                rwidth=1,
                # label=name
                )
        ax.set_xlabel('Interaction strength')
        if idx==0:
            ax.set_ylabel('Model count')
        else:
            ax.set_ylabel('')
            ax.set_yticks([])
        ax.set_title(links['Regulator'][idx]+'-'+links['Target'][idx])
        # ax.set_ylim([0,35])
        ax.set_xlim([0.1, 0.45])
        if idx == top_n-1:
            handles = []
            for i, color in enumerate([color]):
                handles.append(ax.scatter([], [], marker='o', label=name, s=50,
                                          edgecolor='black', color=color, linewidth=.2))

            ax.legend(handles=handles,
                      bbox_to_anchor=(1, .8), prop={'size': 12}, loc='right',
                      frameon = False
                      )
        # ax.set_xmargin(.15)
    return fig
if __name__ == '__main__':
    method = 'RF'
    colors = ['lightblue', 'pink']
    for DE_type, _ in F_DE_protnames().items():
        for idx, study in enumerate(['ctr','mg']):
            links_pool = pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{DE_type}_{study}.csv'))
            fig = plot_weight(links_pool, color=colors[idx], name=study, dist_key='WeightPool')
            fig.savefig(os.path.join(GRN_DIR, f'weights_dist_{DE_type}_{study}.pdf'))
            fig.savefig(os.path.join(GRN_DIR, f'weights_dist_{DE_type}_{study}.png'), dpi=300, transparent=True)


