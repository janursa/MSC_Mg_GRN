import sys
import os
import numpy as np
import pandas as pd
import tqdm
from typing import List

from proteomics_MSC.imports import STRING_DIR, GRN_DIR
from proteomics_MSC.common_tools import F_DE_data
from proteomics_MSC.common_tools.model_selection import save_baseline_scores, create_random_links, calculate_early_precision

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

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
    methods = config['GRN_methods']
    n_repeat = config['n_random_links']
    target_study = config['model_selection_study']
    top_quantiles = config['top_quantiles']

    for DE_type in F_DE_data().keys():
        # - retrieve links
        links_methods = []
        for method in methods:
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{target_study}.csv'), index_col=False)
            links_methods.append(links)
        # - calculate random links
        links_random = create_random_links(links_methods, n_repeat)
        # - get the golden links
        golden_links = pd.read_csv(os.path.join(STRING_DIR, f'links_{DE_type}.csv'), index_col=False)
        # - calculate epr (10) and epr_mean for random links
        ep_scores_series = calculate_ep_for_random_links(links_random, DE_type)  # 1000*10
        assert (np.asarray(ep_scores_series).shape == (n_repeat,len(top_quantiles)))
        ep_scores = np.mean(ep_scores_series, axis=1)  # 1000
        assert (np.asarray(ep_scores).shape == (n_repeat,))

        save_baseline_scores(ep_scores, ep_scores_series, DE_type)
