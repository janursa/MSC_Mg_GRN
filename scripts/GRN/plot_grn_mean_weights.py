"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, F_DE_protiens
from scripts.utils.links import plot_mean_weights

def retreive_grn(method, DE_type, studies):
    links_stack = []
    for study in studies:
        links_stack.append(pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False))

    return links_stack
if __name__ == '__main__':
    studies = ['ctr', 'mg']
    methods = ['RF', 'Portia', 'Ridge']

    for DE_type, _ in F_DE_protiens().items():
        links_combined = []
        for method in methods:
            links_combined.append(retreive_grn(method, DE_type, studies))
        colors = ['lightblue', 'pink']

        fig = plot_mean_weights(links_combined, methods, colors, studies)
        fig.savefig(os.path.join(GRN_DIR, f'mean_weights_{DE_type}.pdf'))

    # plt.show()
