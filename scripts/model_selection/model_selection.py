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
from utils import make_title_pretty
from imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, CALIBRATION_DIR, RANDOM_MODELS_DIR, top_quantiles
from utils import calibration
from utils.model_selection import lineplot, violinplot, create_random_links, calculate_early_precision

score_type: TypeAlias = Dict[str, Any]  # template to store score. score_name, e.g. ep: score value




def extract_tag_based_data_from_modeldata(data: Dict[str, pd.DataFrame], tag: str) -> List[pd.DataFrame]:
    """
    Since data is given for each model name, e.g. early_KNN_Portia, this function extracts those items that has specific name in them.
    """
    selected = []
    for key, value in data.items():
        if tag in key:
            selected.append(value)
    return selected



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


def lineplot_all_models(scores_series, top_quantiles, line_names):
    ncols = 2
    nrows = 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True,
                                           figsize=(3 * ncols, 2.25 * nrows))

    for idx, DE_type in enumerate(DE_types):
        i = int(idx / ncols)
        j = idx % ncols
        ax_series = axes[i][j]
        eprs = extract_tag_based_data_from_modeldata(scores_series, DE_type)
        epr_scores_all = np.vstack([[1 for _ in top_quantiles], eprs])
        lineplot(ax=ax_series, idx=idx, x_data=top_quantiles, data_stack=epr_scores_all,
                 line_names=line_names, title=make_title_pretty(DE_type), yticks=[0, 1, 2])
        if idx == 1:
            ax_series.legend(loc='upper center', bbox_to_anchor=(1.35, 1), fancybox=False, frameon=False)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_all_models.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_all_models.pdf'))
def get_baseline_scores(data_type):
    """Retreive baseline (random) scores of ep-mean"""
    ep_mean_scores_random = np.loadtxt(Path(RANDOM_MODELS_DIR) / f'ep_scores_{data_type}.csv', delimiter=',')
    return ep_mean_scores_random
def calculate_scores() -> Dict[str, score_type]:
    all_scores: Dict[str, score_type] = {} # model_name: dict of different scores
    for DE_type in F_DE_data().keys():
        # - get the baseline scores
        ep_scores_random = get_baseline_scores(DE_type) #1000
        ep_score_random = np.mean(ep_scores_random)
        all_scores['_'.join([DE_type, 'baseline'])] = {'ep':ep_scores_random,
                                                       'epr':ep_scores_random/np.mean(ep_scores_random),
                                                       'percentage':None, 'sig_flag':None}
        # - get the golden links
        golden_links = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
        #- scores of all DE_type
        for method in GRN_methods:
            # - model name
            model_name = '_'.join([DE_type, method])
            # - get the links
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{target_study}.csv'), index_col=False)
            # - calculate epr for each model
            ep_scores = calculate_early_precision(links, golden_links, top_quantiles)
            ep_score = np.mean(ep_scores)
            # - calculate what percentage of random values are bigger than random score
            count = sum(1 for x in ep_scores_random if
                        x > ep_score)  # Count how many elements in 'a' are bigger than 'b'
            percentage = int((count / len(ep_scores_random)) * 100)  # Calculate the percentage
            # - check if AUC values are sig larger than random values
            s, p = scipy.stats.ttest_1samp(ep_scores_random, ep_score)
            if (p < 0.05) & (ep_score > np.mean(ep_scores_random)):
                sig_flag = True
            else:
                sig_flag = False
            # - get the test score
            # test_scores_method = retrieve_test_scores(DE_type, method, studies)  # test score is calculated for both ctr and sample
            all_scores[model_name] = {'ep':ep_score, 'epr':ep_score/ep_score_random, 'sig_flag':sig_flag, 'percentage':percentage}

    return all_scores
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--GRN_methods', nargs='+', default=['RF', 'ridge', 'portia'])
    parser.add_argument('--target_study', type=str, default='ctr')
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    args, remaining_args = parser.parse_known_args()

    GRN_methods = args.GRN_methods
    target_study = args.target_study
    studies = args.studies
    top_quantiles = np.linspace(.75, .9, 10)
    #- to be displayed on the graph
    methods_preferred_names = {'RF':'RF', 'ridge':'Ridge', 'portia':'Portia', 'baseline':'Baseline'}

    #- calculate all scores for all models
    scores = calculate_scores() # model_name: scores[dict]

    #- plot epr: line and violin
    # titles = [make_title_pretty(DE_type) for DE_type in DE_types]
    # lineplot_all_models(ep_scores_all_models, top_quantiles, line_names= methods_preferred_names)
    violinplot_all_models(scores, 1000 ,methods_preferred_names)

    # #- filter those that have non sig AUC vs random
    # shortlisted_model_names = list(filter(lambda x: sig_AUC_vs_random[x], sig_AUC_vs_random))
    #
    # #- filter based on test score: it should be bigger than 0 across both ctr and sample
    # shortlisted_model_names = list(filter(lambda key:(test_scores[key] is None) or all(study_score>0 for study_score in test_scores[key]), shortlisted_model_names))
    # np.savetxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), shortlisted_model_names, delimiter=",", fmt="%s")
    #
    # #- plot the filtered models. line plot
    # shortlisted_scores = [epr_scores_series[name] for name in shortlisted_model_names]
    # shortlisted_scores_random_included = np.vstack([[1 for _ in top_quantiles], shortlisted_scores])
    #
    # line_names = [f'{make_title_pretty(name)} ({round(np.sum(AUC_scores[name]), 2)})' for name in
    #                                   shortlisted_model_names]
    # line_names_random_included = [f'Arbitrary (1)'] + line_names
    #
    # fig_main, ax_main = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.5))
    # lineplot(ax=ax_main, idx=0, x_data=top_quantiles, data_stack=shortlisted_scores_random_included,
    #           line_names=line_names_random_included, title='', yticks= None)
    # ax_main.legend(loc='upper center', bbox_to_anchor=(1.54, 1), fancybox=False, frameon=False, title='Model (AUC)', title_fontproperties={'weight':'bold'} )
    #
    # fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_shortlisted_models.png'), bbox_inches="tight", dpi=300, transparent=True)
    # fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_shortlisted_models.pdf'), bbox_inches="tight")
    #
    # #- select top performing model for early and late
    # def find_best_model(phase):
    #     phase_model_names = list(filter(lambda name: name.split('_')[0]==phase, shortlisted_model_names))
    #     phase_scores = [AUC_scores[name] for name in phase_model_names]
    #     best_model_name = phase_model_names[np.argmax(phase_scores)]
    #     return best_model_name
    # best_model_names = [find_best_model('early'), find_best_model('late')]
    # np.savetxt(os.path.join(MODELSELECTION_DIR, f'selected_models.txt'), best_model_names, delimiter=",", fmt="%s")


