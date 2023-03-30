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


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from typing import Dict, List, Tuple, Callable
from utils import make_title_pretty


from imports import ENRICH_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, CALIBRATION_DIR
from utils import calibration
from utils.model_selection import lineplot, violinplot, create_random_links, calculate_early_precision
def retrieve_test_scores(DE_type: str, method: str, studies: List[str]) -> List[float]:
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

def retrieve_links(DE_types: List[str], methods: List[str], model_name_F: callable) -> Tuple[Dict[str, pd.DataFrame], Dict[str, List[pd.DataFrame]]]:
    """
    Retrieves regulatory links for ridge, RF, and Portia. Using these links, it creates random links
    """
    method_links_stack: Dict[str, pd.DataFrame] = {} # {model:links} # only for ctr -> we need it for epr calculation
    random_links_stack: Dict[str, List[pd.DataFrame]] = {} # {DE_type: [n*links]}
    for DE_type in DE_types:
        links_methods = []
        for method in methods:
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_ctr.csv'), index_col=False)
            method_links_stack[model_name_F(DE_type, method)] = links
            links_methods.append(links)
        random_links = create_random_links(links_methods, n_repeat)
        random_links_stack[DE_type] = random_links
    return method_links_stack, random_links_stack

def extract_tag_based_data_from_modeldata(data: Dict[str, pd.DataFrame], tag: str) -> List[pd.DataFrame]:
    """
    Since data is given for each model name, e.g. early_KNN_Portia, this function extracts those items that has specific name in them.
    """
    selected = []
    for key, value in data.items():
        if tag in key:
            selected.append(value)
    return selected
def calculate_early_precision_ratio_random(golden_links: pd.DataFrame, links_random: List[pd.DataFrame],
                                            name: str, top_quantiles: List[float]) -> Tuple[List[float], List[List[float]]]:
    """
    Calculates EP ratio for randomly generated links.

    Outputs:
        ep_scores_random_series: EP scores calculated for each random links for each top quantile value
        AUC_scores_random: AUC calculated for each random links, summing on epr scores across all top quantiles
    """
    AUC_scores_random: List[float] = []
    ep_scores_random_series: List[List[float]] = []  # 100*10
    # - epr scores for random
    for i in tqdm.tqdm(range(len(links_random)), desc=f'Progress for {name}'):
        scores = calculate_early_precision(links_random[i], golden_links, top_quantiles) #10 scores one for each top quantile
        ep_scores_random_series.append(scores)
        AUC_scores_random.append(np.sum(scores))
    return AUC_scores_random, ep_scores_random_series


def calculate_scores(DE_types: List[str], methods: List[str], studies: List[str], top_quantiles: List[int],
                     model_name_F: Callable[[str, str], str]) -> Tuple[
    Dict[str, Dict[str, float]], Dict[str, float], Dict[str, np.ndarray], Dict[str, Dict[str, np.ndarray]], Dict[
        str, bool]]:
    """ The main function to calculate scores for model selection.
    Outputs:
        test_scores: test scores for regression models
        epr_scores_series: EP ratio for each GRN method, for each top quantile
        sig_AUC_vs_random: whether AUC values for GRN methods are significantly larger than AUC values obtained for random links
        AUC_scores_n:  AUC scores for GRN methods normalized with mean random AUC
        AUC_scores_random_n:  AUC scores for random links normalzied with mean random AUC

    """
    # - retrieve links and calculate random links
    links_methods_stack, links_random_stack = retrieve_links(DE_types, methods, model_name_F)

    # - calculate test score and epr for each model
    test_scores: Dict[str, Dict[str, float]] = {}  # test scores for all models
    AUC_scores_n: Dict[str, float] = {}  # AUC values for all models. normalized to AUC_scores_random
    epr_scores_series: Dict[
        str, Dict[str, np.ndarray]] = {}  # for line plot of epr vs top quantiles. DE_type: scores_series_methods
    AUC_scores_random_n: Dict[
        str, np.ndarray] = {}  # random AUC values for each de_type. normalized to AUC_scores_random
    sig_AUC_vs_random: Dict[str, bool] = {}  # whether AUC values are significantly larger than random values

    for DE_type in DE_types:
        # - get the golden links
        golden_links = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)

        # - calculate ep and AUC for random -> for each DE_type and not for each model
        AUC_scores_random_detype, ep_scores_random_series_detype = calculate_early_precision_ratio_random(golden_links, links_random_stack[DE_type], DE_type, top_quantiles)

        # - calculate and store AUC for de_type
        AUC_scores_random_detype_n = AUC_scores_random_detype / np.mean(AUC_scores_random_detype)
        AUC_scores_random_n[DE_type] = AUC_scores_random_detype_n

        for method_i, method in enumerate(methods):
            model_name = model_name_F(DE_type, method)

            # - get the test score
            test_scores[model_name] = retrieve_test_scores(DE_type, method, studies)  # test score is calculated for both ctr and sample

            # - calculate AUC and epr for each model
            ep_scores_series_method = calculate_early_precision(links_methods_stack[model_name], golden_links,
                                                                top_quantiles)
            AUC_scores_method = np.sum(ep_scores_series_method)
            AUC_scores_method_n = AUC_scores_method / np.mean(AUC_scores_random_detype)
            epr_scores_series_method = ep_scores_series_method / np.mean(ep_scores_random_series_detype,
                                                                         axis=0)  # this is divided by the mean value of line

            # - store AUC and epr values
            AUC_scores_n[model_name] = AUC_scores_method_n  # only calculated for ctr
            epr_scores_series[model_name] = epr_scores_series_method

            # - check if AUC values are sig larger than random values
            s, p = scipy.stats.ttest_1samp(AUC_scores_random_detype_n, AUC_scores_method_n)
            if (p < 0.05) & (AUC_scores_method_n > np.mean(AUC_scores_random_detype_n)):
                sig_AUC_vs_random[model_name] = True
            else:
                sig_AUC_vs_random[model_name] = False

    return test_scores, AUC_scores_n, AUC_scores_random_n, epr_scores_series, sig_AUC_vs_random


def violinplot_all_models(AUC_scores, AUC_scores_random, sig_flags, methods_names):
    ncols = 2
    nrows = 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.3 * ncols, 2.1 * nrows))

    for idx, DE_type in enumerate(DE_types):
        random_values = AUC_scores_random[DE_type]
        AUC_scores_methods = extract_tag_based_data_from_modeldata(AUC_scores, DE_type)
        sig_flags_methods = extract_tag_based_data_from_modeldata(sig_flags, DE_type)
        i = int(idx / ncols)
        j = idx % ncols
        ax_dist = axes[i][j]
        sig_signs = [r'$*$' if flag else '' for flag in sig_flags_methods]
        sig_signs.insert(0, '')  # for random
        scores_dist = np.vstack([random_values, [np.repeat(score, n_repeat) for score in AUC_scores_methods]])
        violinplot(ax=ax_dist, idx=idx, data_stack=scores_dist, x_labels=methods_names, sig_signs=sig_signs, title=make_title_pretty(DE_type))

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--GRN_methods', nargs='+', default=['RF', 'ridge', 'portia'])
    parser.add_argument('--n_random_links', default=1000, type=int, help='Number of randomly generated noisy links')
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    methods = args.GRN_methods
    n_repeat = args.n_random_links
    top_quantiles = np.linspace(.75, .9, 10)
    #- to be displayed on the graph
    methods_preferred_names = ['Arbitrary', 'RF','Ridge', 'Portia']

    DE_types = F_DE_data().keys()

    model_name_F = lambda DE_type, method: '_'.join([DE_type, method])

    # - create directories
    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)

    #- calculate test score, epr, AUC, and sig flags
    test_scores, AUC_scores, AUC_scores_random, epr_scores_series, sig_AUC_vs_random = \
        calculate_scores(DE_types, methods, studies, top_quantiles, model_name_F)

    #- plot epr: line and violin
    titles = [make_title_pretty(DE_type) for DE_type in DE_types]
    lineplot_all_models(epr_scores_series, top_quantiles, line_names= methods_preferred_names)
    violinplot_all_models(AUC_scores, AUC_scores_random, sig_flags=sig_AUC_vs_random, methods_names=methods_preferred_names)

    #- filter those that have non sig AUC vs random
    shortlisted_model_names = list(filter(lambda x: sig_AUC_vs_random[x], sig_AUC_vs_random))

    #- filter based on test score: it should be bigger than 0 across both ctr and sample
    shortlisted_model_names = list(filter(lambda key:(test_scores[key] is None) or all(study_score>0 for study_score in test_scores[key]), shortlisted_model_names))
    np.savetxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), shortlisted_model_names, delimiter=",", fmt="%s")

    #- plot the filtered models. line plot
    shortlisted_scores = [epr_scores_series[name] for name in shortlisted_model_names]
    shortlisted_scores_random_included = np.vstack([[1 for _ in top_quantiles], shortlisted_scores])

    line_names = [f'{make_title_pretty(name)} ({round(np.sum(AUC_scores[name]), 2)})' for name in
                                      shortlisted_model_names]
    line_names_random_included = [f'Arbitrary (1)'] + line_names

    fig_main, ax_main = plt.subplots(nrows=1, ncols=1, figsize=(3.5, 2.5))
    lineplot(ax=ax_main, idx=0, x_data=top_quantiles, data_stack=shortlisted_scores_random_included,
              line_names=line_names_random_included, title='', yticks= None)
    ax_main.legend(loc='upper center', bbox_to_anchor=(1.54, 1), fancybox=False, frameon=False, title='Model (AUC)', title_fontproperties={'weight':'bold'} )

    fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_shortlisted_models.png'), bbox_inches="tight", dpi=300, transparent=True)
    fig_main.savefig(os.path.join(MODELSELECTION_DIR, f'lineplots_shortlisted_models.pdf'), bbox_inches="tight")

    #- select top performing model for early and late
    def find_best_model(phase):
        phase_model_names = list(filter(lambda name: name.split('_')[0]==phase, shortlisted_model_names))
        phase_scores = [AUC_scores[name] for name in phase_model_names]
        best_model_name = phase_model_names[np.argmax(phase_scores)]
        return best_model_name
    best_model_names = [find_best_model('early'), find_best_model('late')]
    np.savetxt(os.path.join(MODELSELECTION_DIR, f'selected_models.txt'), best_model_names, delimiter=",", fmt="%s")


