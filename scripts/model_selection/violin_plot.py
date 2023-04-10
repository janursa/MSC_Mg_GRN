"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Callable, Optional, TypeAlias, Any



SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from utils import make_title_pretty
from imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, CALIBRATION_DIR, RANDOM_MODELS_DIR
from utils import calibration
from utils.model_selection import lineplot, violinplot, create_random_links, calculate_early_precision

score_type: TypeAlias = Dict[str, Any]  # template to store score. score_name, e.g. ep: score value


def is_single_value(value):
    collection_types = (list, tuple, set, dict, np.ndarray)
    return not isinstance(value, collection_types)
def violinplot_all_models(scores, n_repeat, methods_preferred_names):
    methods_names = [ methods_preferred_names[score.split('_')[2]] for score in scores.keys()] # only method name

    ncols = 1
    nrows = 1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(7.5 * ncols, 3 * nrows))

    scores_list = [item['epr'] for item in scores.values()] 
    sig_flags = [item['sig_flag'] for item in scores.values()]
    percentages = [item['percentage'] for item in scores.values()]

    sig_signs = [r'$*$' if flag else '' for flag in sig_flags]
    percentages_methods = [f'{item}%' if item else '' for item in percentages]

    scores_dist = [np.repeat(score, n_repeat) if (is_single_value(score)) else score for score in scores_list]
    violinplot(ax=ax, data_stack=np.asarray(scores_dist), x_labels=methods_names, sig_signs=sig_signs, percentages=percentages_methods, title='')

    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot_all_models.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot_all_models.pdf'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_random_links', default=1000, type=int, help='Number of randomly generated noisy links')
    args, remaining_args = parser.parse_known_args()

    n_repeat = args.n_random_links

    #- to be displayed on the graph
    methods_preferred_names = {'RF':'RF', 'ridge':'Ridge', 'portia':'Portia', 'baseline':'Baseline'}

    #- retreieve scores
    with open(f'{MODELSELECTION_DIR}/scores.json', 'r') as f:
        # load the JSON data from the file
        scores = json.load(f)

    violinplot_all_models(scores, n_repeat, methods_preferred_names)

