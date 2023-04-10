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
from imports import MODELSELECTION_DIR, make_title_pretty
from md_aux import retreieve_scores
from md_aux import violinplot

def is_single_value(value):
    collection_types = (list, tuple, set, dict, np.ndarray)
    return not isinstance(value, collection_types)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_random_links', default=1000, type=int, help='Number of randomly generated noisy links')
    args, remaining_args = parser.parse_known_args()

    n_repeat = args.n_random_links

    #- retreieve scores
    scores = retreieve_scores()
    # - plot violin
    methods_names = [make_title_pretty(score.split('_')[2]) for score in scores.keys()] # only method name
    scores_list = [item['epr'] for item in scores.values()]
    sig_flags = [item['sig_flag'] for item in scores.values()]
    percentile_ranks = [item['percentile_rank'] for item in scores.values()]

    sig_signs = [r'$*$' if flag else '' for flag in sig_flags]
    percentile_ranks = [f'{item}%' if item else '' for item in percentile_ranks]

    scores_dist = [np.repeat(score, n_repeat) if (is_single_value(score)) else score for score in scores_list]

    fig = violinplot(scores_dist, percentile_ranks, sig_signs, methods_names)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'violinplot.pdf'))

