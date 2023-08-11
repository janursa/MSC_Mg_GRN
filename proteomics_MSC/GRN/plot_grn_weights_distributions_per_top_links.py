"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import matplotlib.pyplot as plt
import argparse
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from _imports import GRN_DIR, F_DE_protnames, study_colors, PLOT_WEIGHTS_DIR
from common_tools.links import plot_grn_weights_distributions_per_top_links

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--GRN_method', type=str, default='RF')
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    method = args.GRN_method


    for DE_type, _ in F_DE_protnames().items():
        for idx, study in enumerate(['ctr','mg']):
            links_pool = pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{DE_type}_{study}.csv'))
            fig = plot_grn_weights_distributions_per_top_links(links_pool, color=study_colors[idx], name=study, dist_key='WeightPool')
            fig.savefig(os.path.join(PLOT_WEIGHTS_DIR, f'weights_dist_{DE_type}_{study}.pdf'))
            fig.savefig(os.path.join(PLOT_WEIGHTS_DIR, f'weights_dist_{DE_type}_{study}.png'), dpi=300, transparent=True)


