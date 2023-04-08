"""
    Create baseline models (random models) for each data types (DE_type)
"""
import sys
import os
import numpy as np
import pandas as pd
import tqdm
import argparse
from pathlib import Path
from typing import List

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from imports import ENRICH_DIR, F_DE_data, GRN_DIR, RANDOM_MODELS_DIR
from utils.model_selection import create_random_links, calculate_early_precision


def calculate_ep_for_random_links(links_stack: List[pd.DataFrame], data_type: str) -> List[List[float]]:
    """calculate early precision ratio for random links"""
    ep_scores_stack: List[List[float]] = []  # 1000 (replica)*10(cut-offs)
    # - ep scores for random
    for i in tqdm.tqdm(range(len(links_stack)), desc=f'Progress for {data_type}'):
        scores = calculate_early_precision(links_stack[i], golden_links,
                                           top_quantiles)  # 10 scores one for each top quantile
        ep_scores_stack.append(scores)
    return ep_scores_stack



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--GRN_methods', nargs='+', default=['RF', 'ridge', 'portia'])
    parser.add_argument('--target_study', type=str, default='ctr')
    parser.add_argument('--n_random_links', default=1000, type=int, help='Number of randomly generated noisy links')
    args, remaining_args = parser.parse_known_args()

    methods = args.GRN_methods
    n_repeat = args.n_random_links
    top_quantiles = np.linspace(.75, .9, 10)
    target_study = args.target_study

    if not os.path.isdir(RANDOM_MODELS_DIR):
        os.makedirs(RANDOM_MODELS_DIR)

    for DE_type in F_DE_data().keys():
        # - retrieve links
        links_methods = []
        for method in methods:
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{target_study}.csv'), index_col=False)
            links_methods.append(links)
        # - calculate random links
        links_random = create_random_links(links_methods, n_repeat)
        # - get the golden links
        golden_links = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
        # - calculate epr (10) and epr_mean for random links
        ep_scores_series = calculate_ep_for_random_links(links_random, DE_type)  # 1000*10
        assert (np.asarray(ep_scores_series).shape == (n_repeat,len(top_quantiles)))
        ep_scores = np.mean(ep_scores_series, axis=1)  # 1000
        assert (np.asarray(ep_scores).shape == (n_repeat,))

        np.savetxt(f'{Path(RANDOM_MODELS_DIR)}/ep_scores_{DE_type}.csv', np.asarray(ep_scores),
                   delimiter=',')
