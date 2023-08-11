"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt


from proteomics_MSC.imports import MODELSELECTION_DIR
from proteomics_MSC.common_tools import make_title_pretty
from proteomics_MSC.common_tools.model_selection import retreieve_scores, is_single_value, ViolinPlot


import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    n_repeat = config['n_random_links']

    #- retreieve scores
    scores = retreieve_scores()
    # - plot violin
    methods_names = [make_title_pretty(score.split('_')[2]) for score in scores.keys()] # only method name
    scores_list = [item['epr'] for item in scores.values()]
    sig_flags = [item['sig_flag'] for item in scores.values()]
    percentile_ranks = [item['percentile_rank'] for item in scores.values()]

    scores_dist = np.asarray([np.repeat(score, n_repeat) if (is_single_value(score)) else score for score in scores_list])
    assert all(len(item) == len(scores_dist[0]) for item in scores_dist), "All items in scores_list must have the same length"


    fig, ax = plt.subplots(nrows=1, ncols=1, tight_layout=True, figsize=(6, 2.2))
    ViolinPlot.plot(ax, scores_dist, percentile_ranks, sig_flags, methods_names)
    to_save = os.path.join(MODELSELECTION_DIR, f'violinplot.png')
    fig.savefig(to_save, dpi=300, transparent=True)
    print(f'output -> {to_save}')
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot.pdf'))

