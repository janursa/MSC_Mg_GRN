"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import typing
import random


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
sys.path.append(os.path.join(os.path.dirname(SCRIPT_DIR), '..'))
sys.path.append(os.path.join(os.path.dirname(SCRIPT_DIR), '..', '..'))

from typing import Dict

import tqdm

from scripts.imports import ENRICH_DIR, VS_STRING_DIR, F_DE_data, GRN_DIR, F_DE_protiens, CALIBRATION_DIR
from scripts.utils.links import plot_match_counts_series, normalize_links, format_links_string
from scripts.utils import calibration, serif_font, make_title_pretty

from geneRNI.evaluation import precision_recall_curve,calculate_auc_roc, calculate_PR


def create_random_links(links_assembly, n=1000):
    #- TODO: create matrix (using protname)
    links_assembly = [normalize_links(links) for links in links_assembly]
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j] #flatten
    sample_links = links_assembly[0]
    random_links = pd.DataFrame({key:sample_links[key] for key in ['Regulator', 'Target']})
    weightpoolvector = []
    for i in range(len(links_assembly[0])):
        weightpoolvector.append(random.sample(weights, n))
    random_links['Weight']= np.mean(weightpoolvector, axis=1)
    random_links_pool = random_links.copy()
    random_links_pool['WeightPool'] = weightpoolvector
    return random_links, random_links_pool


class PermutationTest:

    def __init__(self, golden_links: pd.DataFrame):
        self.golden_links: pd.DataFrame = golden_links
        self.background_aurocs: List[float] = []
        self.background_auprs: List[float] = []

    def add_random_predictions(self, links: pd.DataFrame):
        self.background_aurocs.append(calculate_auc_roc(links, self.golden_links))
        self.background_auprs.append(calculate_PR(links, self.golden_links))

    @staticmethod
    def compute_one_sided_p_value(background: np.ndarray, value: float) -> float:
        return np.mean(value < background)

    def compute_overall_score(self, links: pd.DataFrame) -> Dict[str, float]:
        auroc = calculate_auc_roc(links, self.golden_links)
        p_auroc = PermutationTest.compute_one_sided_p_value(
            np.asarray(self.background_aurocs),
            auroc
        )
        aupr = calculate_PR(links, self.golden_links)
        p_aupr = PermutationTest.compute_one_sided_p_value(
            np.asarray(self.background_auprs),
            aupr
        )
        overall_score = -0.5 * (np.log10(p_auroc) + np.log10(p_aupr))
        return {
            'auroc': auroc,
            'aupr': aupr,
            'p-value-auroc': p_auroc,
            'p-value-aupr': p_aupr,
            'overall-score': overall_score
        }


def remove_mg(links):
    return links.loc[(links['Regulator']!='mg') & (links['Target']!='mg'),:]


if __name__ == '__main__':
    if not os.path.isdir(VS_STRING_DIR):
        os.makedirs(VS_STRING_DIR)
    n_repeat = 10

    arbitrary_dir = os.path.join(GRN_DIR, 'arbitrary')
    if not os.path.isdir(arbitrary_dir):
        os.makedirs(arbitrary_dir)

    DE_types = F_DE_data().keys()
    study='all-in'
    methods = ['ridge', 'portia', 'arbitrary']
    # methods = ['portia', 'arbitrary']

    # ncols = 2
    # nrows = 3
    # fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(5 * ncols, 4 * nrows))

    links_string_dict = {}
    for DE_type in DE_types:
        links_string = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
        links_string_dict[DE_type] = links_string

    # serif_font()
    # colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey', 'lightcoral']
    # linestyles = ['-', '--', '-.', ':', '-', '--']

    for idx, DE_type in enumerate(DE_types):

        # Load golden links
        golden_links = format_links_string(links_string_dict[DE_type], F_DE_protiens()[DE_type])

        permutation_test = PermutationTest(golden_links)

        # i = int(idx/ncols)
        # j = idx%ncols
        # ax = axes[i][j]
        # - retreive or create the links
        links_stack = []

        for method in methods:
            if method in ['portia', 'ridge']:
                # - single weight
                links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
                links_stack.append(links)
            elif method == 'RF':
                links_pool = pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{DE_type}_{study}.csv'))
                links_stack.append(links_pool)
            elif method == 'arbitrary':
                for _ in tqdm.tqdm(range(200), desc='Making permutations'):
                    links, links_pool = create_random_links(links_stack, n=n_repeat)
                    permutation_test.add_random_predictions(remove_mg(links))
                links.to_csv(os.path.join(arbitrary_dir, f'links_{DE_type}_{study}.csv'), index=False)
                links_pool.to_pickle(os.path.join(arbitrary_dir, f'links_pool_{DE_type}_{study}.csv'))
                links_stack.append(links_pool)
        # assert (len(links_stack)==2)
        #-- compare to golden links
        protnames = F_DE_protiens()[DE_type]
        # - get the links, filter mask based on scores, and calculate recall
        for method_i, method in enumerate(methods):
            if method == 'arbitrary':
                continue

            links = links_stack[method_i]
            links = remove_mg(links)
            # if method in ['ridge', 'RF']: #skipped
            #     # best_scores, best_params = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
            #     # protnames_left = np.asarray(protnames)[np.asarray(best_scores > 0)]
            #     # filter_mask = links['Target'].isin(protnames_left)
            #     filter_mask = None
            # else:
            #     filter_mask = None

            _, recall, _ = precision_recall_curve(links, golden_links)
            # calculate_auc_roc(links, golden_links)

            # label = method + f' (score = {str(round(np.mean(recall[:-1]), 2))})'
            # def normalize(aa):
            #     return (aa
            # ax.plot(, recall, label=label, color=colors[i], alpha=1, linewidth=2,
            #             linestyle=linestyles[i])
            # print(f'{DE_type} -- {method}---{round(np.mean(recall[:-1]), 2)}')

            results = permutation_test.compute_overall_score(links)
            print(f'Results for {method}: {results}')

        # ax.set_title(DE_type)
        # ax.legend(frameon=False)
        # ax.set_ylabel('Recall')
        # ax.set_xlabel('Thresholds')

    # fig_pr.savefig(os.path.join(VS_STRING_DIR, f'precision_recall.png'), dpi=300, transparent=True)
    # fig.savefig(os.path.join(VS_STRING_DIR, f'recall.pdf'))




