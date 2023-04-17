import argparse
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt

import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from imports import  F_selected_models, RANDOM_REGULATORS_DIR
from common_tools.model_selection import determine_sig_flag, retreieve_scores, is_single_value, ViolinPlot
from rr_aux import retrieve_random_scores
def calculate_scores_vs_random(selected_models):
    # - retreive ep scores for the real models
    scores = {key: value['ep_series'] for key, value in retreieve_scores().items()}

    percentile_ranks = [] # for both models
    ep_scores = []
    sig_flags = []
    for selected_model in selected_models:
        # - retreive ep scores from the random models
        ep_scores_random = retrieve_random_scores(selected_model)
        # - calculate what percentage of random values are bigger than random score
        ep_score = np.mean(scores[selected_model])
        count = sum(1 for x in ep_scores_random if
                    x > ep_score)  # Count how many elements in 'a' are bigger than 'b'
        percentile_rank = int((count / len(ep_scores_random)) * 100)  # Calculate the percentile_rank
        sig_flag = determine_sig_flag(ep_scores_random, ep_score)

        # - append scores: random, model
        ep_scores.extend([ep_scores_random, ep_score])
        percentile_ranks.extend([None, percentile_rank])
        sig_flags.extend([None, sig_flag])
    return percentile_ranks, ep_scores, sig_flags

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_repeat', type=int, default=100, help="How many random selection")
    args, remaining_args = parser.parse_known_args()
    n_repeat = args.n_repeat

    # - load names of selected models
    selected_models = F_selected_models()
    # - get ep scores and data for plot
    percentile_ranks, ep_scores, sig_flags = calculate_scores_vs_random(selected_models)
    # - repeat scores to create distribution similar to random
    scores_dist = np.asarray([np.repeat(score, n_repeat) if is_single_value(score) else score for score in ep_scores])
    assert all(len(item) == len(scores_dist[0]) for item in scores_dist), "All items in scores_list must have the same length"

    # - pile the scores of random and real scores
    methods_names = np.tile(['Random', 'DE'], 2)
    # - plot
    colors = np.repeat(ViolinPlot.group_colors, 2)
    x_ticks = [1, 2, 3.5, 4.5]
    x_margin = .05
    fig, ax = plt.subplots(nrows=1, ncols=1, tight_layout=True, figsize=(2.2, 2.2))
    ViolinPlot.plot(ax, scores_dist, percentile_ranks, sig_flags, methods_names, yaxis_name='Early precision (%)', violin_colors=colors, x_ticks=x_ticks, x_margin=x_margin)
    fig.savefig(os.path.join(RANDOM_REGULATORS_DIR, f'violinplot.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(RANDOM_REGULATORS_DIR, f'violinplot.pdf'))





