"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import tqdm
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Callable, Optional, TypeAlias, Any



SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, CALIBRATION_DIR, RANDOM_MODELS_DIR, top_quantiles, make_title_pretty
from md_aux import lineplot, retreieve_scores

def calculate_epr_series(model_names:List[str], ep_series_scores:Dict[str, List[float] | List[List[float]]]) -> List[List[float]]:
    """Calculate epr for a series of top quantiles (10 numbers)
        For each top quantile, we divide ep value to the ep value averaged
        over n baseline values.
    """
    # extract ep series for the baseline model, the one created on the same dataset of each model
    epr_series = []
    for model_name in model_names:
        data_type = '_'.join(model_name.split('_')[0:2])
        baseline_name = '_'.join([data_type, 'baseline'])
        epr_series_random = ep_series_scores[baseline_name]
        epr_series_mean = np.mean(epr_series_random, axis=0)
        epr_series.append(ep_series_scores[model_name]/epr_series_mean)
    return epr_series


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args, remaining_args = parser.parse_known_args()


    # - retreive scores
    scores = retreieve_scores()
    # - get the shortlisted model
    shortlisted_modelnames = list(np.genfromtxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), delimiter=",", dtype=str))
    # - extract ep series for shortlisted models
    ep_scores = [value['ep_series'] for key,value in scores.items() if key in shortlisted_modelnames]
    # - calculate epr series
    ep_series_scores = {key:value['ep_series'] for key,value in scores.items()}
    epr_series = calculate_epr_series(shortlisted_modelnames, ep_series_scores)
    # - get epr mean score, average over all top quantiles, to put into the graph
    epr_mean_scores = [np.round(scores[name]['epr'],2) for name in shortlisted_modelnames]
    # - add random score and name to the data
    assert (len(epr_series) == len(shortlisted_modelnames))
    epr_mean_scores = [1.0] + epr_mean_scores
    models_names = ['baseline'] + shortlisted_modelnames
    models_names = [f'{name} ({score})' for name, score in zip(models_names, epr_mean_scores)]
    epr_baseline = np.repeat(1, len(top_quantiles))
    data_stack = np.vstack([epr_baseline, epr_series])
    # - plot
    models_names = [make_title_pretty(name) for name in models_names]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.2))
    lineplot(ax=ax, x_data=top_quantiles, data_stack=data_stack,
              line_names=models_names, title='', yticks= None)
    ax.legend(loc='upper center', bbox_to_anchor=(1.54, 1), fancybox=False, frameon=False, title='Model (Mean EPR)', title_fontproperties={'weight':'bold'} )

    fig.savefig(os.path.join(MODELSELECTION_DIR, f'lineplot.png'), bbox_inches="tight", dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'lineplot.pdf'), bbox_inches="tight")

