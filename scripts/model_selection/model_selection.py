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
from imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, CALIBRATION_DIR, RANDOM_MODELS_DIR
from utils import calibration
from utils.model_selection import lineplot, violinplot, create_random_links, calculate_early_precision

score_type: TypeAlias = Dict[str, Any]  # template to store score. score_name, e.g. ep: score value
def retrieve_test_scores(DE_type: str, method: str, studies: List[str]) -> Optional[Tuple[float, float]]:
    """
    Retrieves calibration test scores for ridge and RF and sets portia's score to None
    """
    if method == 'portia':
        scores_studies = None
    else:
        scores_studies = []
        for study in studies:
            best_scores, _ = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
            mean_score = np.mean(best_scores)
            scores_studies.append(mean_score)
    return scores_studies



def extract_tag_based_data_from_modeldata(data: Dict[str, pd.DataFrame], tag: str) -> List[pd.DataFrame]:
    """
    Since data is given for each model name, e.g. early_KNN_Portia, this function extracts those items that has specific name in them.
    """
    selected = []
    for key, value in data.items():
        if tag in key:
            selected.append(value)
    return selected




def violinplot_all_models(scores, methods_preferred_names):
    methods_names = ['Baseline'] + [methods_preferred_names[method] for method in GRN_methods]
    ncols = 2
    nrows = 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.3 * ncols, 2.1 * nrows))

    ep_mean_scores = {model_name: item['epr'] for model_name, item in scores.items()}
    sig_flags = {model_name: item['sig_flag'] for model_name, item in scores.items()}
    percentages = {model_name: item['percentage'] for model_name, item in scores.items()}
    for idx, DE_type in enumerate(F_DE_data().keys()):
        random_scores = get_baseline_scores(DE_type)
        AUC_scores_methods = extract_tag_based_data_from_modeldata(ep_mean_scores, DE_type)
        sig_flags_methods = extract_tag_based_data_from_modeldata(sig_flags, DE_type)
        percentages_methods = extract_tag_based_data_from_modeldata(percentages, DE_type)
        i = int(idx / ncols)
        j = idx % ncols
        ax_dist = axes[i][j]
        sig_signs = [r'$*$' if flag else '' for flag in sig_flags_methods]
        percentages_methods = [f'{item}%' for item in percentages_methods]
        sig_signs.insert(0, '')  # for random
        scores_dist = np.vstack([random_scores/np.mean(random_scores), [np.repeat(score, len(random_scores)) for score in AUC_scores_methods]])
        violinplot(ax=ax_dist, idx=idx, data_stack=scores_dist, x_labels=methods_names, sig_signs=sig_signs, percentages=percentages_methods, title=make_title_pretty(DE_type))

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
    methods_preferred_names = {'RF':'RF', 'ridge':'Ridge', 'portia':'Portia'}

    #- calculate all scores for all models
    scores = calculate_scores() # model_name: scores[dict]

    #- plot epr: line and violin
    # titles = [make_title_pretty(DE_type) for DE_type in DE_types]
    # lineplot_all_models(ep_scores_all_models, top_quantiles, line_names= methods_preferred_names)
    violinplot_all_models(scores, methods_preferred_names)

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


