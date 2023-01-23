"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import GRN_DIR
from utils.links import read_write_links, plot_mean_weights

def retreive_grn(method):
    links_ctr= read_write_links(method=method, study='ctr', mode='read', output_dir=GRN_DIR)
    links_mg= read_write_links(method=method, study='mg', mode='read', output_dir=GRN_DIR)
    return links_ctr, links_mg
if __name__ == '__main__':

    links_ctr_rf, links_mg_rf = retreive_grn('RF')
    links_ctr_ridge, links_mg_ridge = retreive_grn('ridge')
    links_ctr_portia, links_mg_portia = retreive_grn('portia')

    links_combined = [[links_ctr_rf, links_mg_rf],
                      [links_ctr_portia, links_mg_ridge],
                      [links_ctr_ridge, links_mg_portia]]

    # - normalize
    labels = ['RF', 'Portia',
              'Ridge']
    colors = ['lightblue', 'pink']

    fig = plot_mean_weights(links_combined, labels, colors)
    fig.savefig(os.path.join(GRN_DIR, 'mean_weights.pdf'))

    # plt.show()
