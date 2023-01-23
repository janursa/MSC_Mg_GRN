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

from imports import GRN_DIR

def plot_weight(links, dist_key='WeightPool'):
    links_s = links.sort_values('Weight', ascending=False).reset_index(drop=True)

    top_n = 12
    nrows = 4
    ncols = 3
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(ncols * 3, nrows * 2))
    for idx in range(top_n):
        i = int(idx / (nrows - 1))
        j = idx % ncols
        ax = axes[i][j]
        ax.hist(links_s[dist_key][idx],
                bins=20,
                alpha=0.5,
                histtype='stepfilled',  # 'bar', 'barstacked', 'step', 'stepfilled'
                color='lightgreen',
                ec='black',
                rwidth=.9
                )
        ax.set_xlabel('Interaction strength')
        ax.set_ylabel('Model count')
        #         title = links['Regulator'][idx]+'-->'+links['Target'][idx]
        title = idx
        ax.set_title(title)
        ax.set_ymargin(.15)
        ax.set_xmargin(.15)
    return fig
if __name__ == '__main__':
    method = 'RF'
    links_pool_ctr, links_pool_mg = (pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{study}.csv')) for study in ['ctr', 'mg'])

    fig = plot_weight(links_pool_ctr, dist_key='WeightPool')


    plt.show()
