"""
    Auxillary functions for model selection
"""
import sys
import os
import numpy as np
import json
from pathlib import Path
from typing import List, TypeAlias, Dict, Any, Tuple
import random
import pandas as pd
from typing import List, Dict


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from imports import MODELSELECTION_DIR, RANDOM_MODELS_DIR
from utils import serif_font, flatten
from utils.links import choose_top_quantile, normalize_links


score_type: TypeAlias = Dict[str, Any]  # template to store score. score_name, e.g. ep: score value

def save_scores(scores: Dict[str, score_type]) -> None:
    with open(f'{MODELSELECTION_DIR}/scores.json', 'w') as f:
        json.dump(scores, f)
def retreieve_scores() -> Dict[str, score_type]:
    with open(f'{MODELSELECTION_DIR}/scores.json', 'r') as f:
        # load the JSON data from the file
        scores = json.load(f)
    return scores
def get_baseline_scores(data_type:str) -> Tuple[List[float], List[List[float]]]:
    """Retreive baseline (random) scores of ep-mean"""
    ep_scores = np.loadtxt(Path(RANDOM_MODELS_DIR) / f'ep_scores_{data_type}.csv', delimiter=',')
    ep_scores_series = np.loadtxt(Path(RANDOM_MODELS_DIR) / f'ep_scores_series_{data_type}.csv', delimiter=',')

    return list(ep_scores), list(ep_scores_series)
def save_baseline_scores(ep_scores:List[float], ep_scores_series:List[List[float]], DE_type:str) -> None:
    np.savetxt(f'{Path(RANDOM_MODELS_DIR)}/ep_scores_{DE_type}.csv', np.asarray(ep_scores),
               delimiter=',')
    np.savetxt(f'{Path(RANDOM_MODELS_DIR)}/ep_scores_series_{DE_type}.csv', np.asarray(ep_scores_series),
               delimiter=',')


def lineplot(ax, x_data, data_stack, line_names, title, yticks):
    """ Line plot for the epr scores versus top quantiles

    """
    serif_font()
    colors = ['grey', 'lightpink', 'lightblue', 'lightgreen', 'blue', 'orange', 'cyan', 'cyan', 'cyan','cyan']
    linestyles = np.repeat(['-', '--', '-.', ':'], 2)
    for i, line_name in enumerate(line_names):
        ax.plot(x_data,
                       data_stack[i],
                       label=line_names[i],
                       color=colors[i],
                       alpha=1,
                       linewidth=2,
                       linestyle=linestyles[i],
                       marker='o')
    ax.set_ylabel('Early precision ratio')
    ax.set_ymargin(.2)
    ax.set_xlabel('Top quantile')
    ax.set_title(title)
    if yticks is not None:
        ax.set_yticks(yticks)

def create_random_links(links_stack: List[pd.DataFrame], n:int) -> List[pd.DataFrame]:
    """
    Creates n number of random links using weights of given set of links
    """
    links_stack = flatten(links_stack)
    links_stack = [normalize_links(links) for links in links_stack]
    weights = [links['Weight'].values.tolist() for links in links_stack]
    weights = [i for j in weights for i in j] #flatten
    sample_links = links_stack[0]
    random_links = []
    for i in range(n):
        rep = sample_links.copy()
        rep.loc[:,'Weight'] = random.sample(weights, len(rep))
        random_links.append(rep)

    return random_links

def calculate_early_precision(links, golden_links, top_quantiles) -> List[float]:
    '''  calculate normalized EP scores for the top_quantiles
        Compare extracted links by GRN to those suggested by model_selection.
        Normalized with length of golden link to make it comparible between cases with different lengths.
    '''
    ns = []
    for tpq in top_quantiles:
        links_top = choose_top_quantile(links, quantile=tpq)

        l_regs = links_top['Regulator'].to_numpy(str)
        l_targs = links_top['Target'].to_numpy(str)

        s_regs = golden_links['Regulator'].to_numpy(str)
        s_targs = golden_links['Target'].to_numpy(str)

        n = 0
        for reg, target in zip(s_regs, s_targs):
            if ((l_regs == reg) & (l_targs == target)).any():
                n+=1
        ns.append(n)
    return list(np.asarray(ns)/len(golden_links)*100)
