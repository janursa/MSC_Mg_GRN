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
import time


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
sys.path.append(os.path.join(os.path.dirname(SCRIPT_DIR), '..'))
sys.path.append(os.path.join(os.path.dirname(SCRIPT_DIR), '..', '..'))

from typing import Dict

import tqdm

from scripts.imports import ENRICH_DIR, VS_STRING_DIR, F_DE_data, GRN_DIR, F_DE_protiens, CALIBRATION_DIR
from scripts.utils.links import choose_top_quantile, normalize_links, format_links_string
from scripts.utils import calibration, serif_font, make_title_pretty

from geneRNI.evaluation import precision_recall_curve,calculate_auc_roc, calculate_PR


def create_random_links(links_assembly):
    #- TODO: create matrix (using protname)
    links_assembly = [normalize_links(links) for links in links_assembly]
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j] #flatten
    sample_links = links_assembly[0]
    random_links = pd.DataFrame({key:sample_links[key] for key in ['Regulator', 'Target']})
    random_links['Weight'] = random.sample(weights, len(random_links))
    return random_links


class PermutationTest:

    def __init__(self, golden_links: pd.DataFrame):
        self.golden_links: pd.DataFrame = golden_links
        # self.background_aurocs: List[float] = []
        # self.background_auprs: List[float] = []
        self.background_ep: List[float] = []

    def add_random_predictions(self, links: pd.DataFrame):
        # self.background_aurocs.append(calculate_auc_roc(links, self.golden_links))
        # self.background_auprs.append(calculate_PR(links, self.golden_links))
        self.background_ep.append(calculate_EP(links, self.golden_links))

    @staticmethod
    def compute_one_sided_p_value(background: np.ndarray, value: float) -> float:
        return np.mean(value < background)

    def compute_overall_score(self, links: pd.DataFrame) -> Dict[str, float]:
        # auroc = calculate_auc_roc(links, self.golden_links)
        # p_auroc = PermutationTest.compute_one_sided_p_value(
        #     np.asarray(self.background_aurocs),
        #     auroc
        # )
        # aupr = calculate_PR(links, self.golden_links)
        # p_aupr = PermutationTest.compute_one_sided_p_value(
        #     np.asarray(self.background_auprs),
        #     aupr
        # )
        ep = calculate_EP(links, self.golden_links)
        p_ep = PermutationTest.compute_one_sided_p_value(
            np.asarray(self.background_ep),
            ep
        )
        # overall_score = -0.5 * (np.log10(p_auroc) + np.log10(p_aupr))
        return {
            # 'auroc': auroc,
            # 'p-value-auroc': p_auroc,
            # 'aupr': aupr,
            # 'p-value-aupr': p_aupr,
            'ep': ep,
            'p-value-ep': p_ep,
            # 'overall-score': overall_score
        }
def _calculate_EP(links, string_df) -> int:
    '''
        Compare extracted links by GRN to those suggested by vs_string. Returns number of match links.
    '''
    # ns = []
    # for tpq in np.linspace(.75, .9, 10):
    tpq = .85
    links_top = choose_top_quantile(links, quantile=tpq)
    n = 0
    for reg, target in zip(string_df['Regulator'].values, string_df['Target'].values):
        if ((links_top['Regulator'] == reg) & (links_top['Target'] == target)).any():
            n+=1
        # ns.append(n)
    # return np.mean(ns)
    return n
def calculate_EP(links, string_df) -> int:
    '''
        Compare extracted links by GRN to those suggested by vs_string. Returns number of match links.
    '''
    ns = []
    for tpq in np.linspace(.75, .9, 10):
        links_top = choose_top_quantile(links, quantile=tpq)

        l_regs = links_top['Regulator'].to_numpy(str)
        l_targs = links_top['Target'].to_numpy(str)

        s_regs = string_df['Regulator'].to_numpy(str)
        s_targs = string_df['Target'].to_numpy(str)

        n = 0
        for reg, target in zip(s_regs, s_targs):
            if ((l_regs == reg) & (l_targs == target)).any():
                n+=1
        ns.append(n)
    return np.mean(ns)

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
    study='ctr'
    methods = ['portia', 'arbitrary']
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
        # if DE_type not in ['late_50']: #TODO: remove
        #     continue
        # Load golden links
        golden_links = format_links_string(links_string_dict[DE_type], F_DE_protiens()[DE_type])

        permutation_test = PermutationTest(links_string_dict[DE_type]) #TODO: golden links is different now

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
                start = time.time()
                for _ in tqdm.tqdm(range(200), desc='Making permutations'):
                    links = create_random_links(links_stack)
                    permutation_test.add_random_predictions(links)
                end = time.time()
                print('lapsed time ', (end-start))
                # links_stack.append(links)
        # assert (len(links_stack)==2)
        #-- compare to golden links
        protnames = F_DE_protiens()[DE_type]
        # - get the links, filter mask based on scores, and calculate recall
        for method_i, method in enumerate(methods):
            if method == 'arbitrary':
                continue

            links = links_stack[method_i]
            # if method in ['ridge', 'RF']: #skipped
            #     # best_scores, best_params = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
            #     # protnames_left = np.asarray(protnames)[np.asarray(best_scores > 0)]
            #     # filter_mask = links['Target'].isin(protnames_left)
            #     filter_mask = None
            # else:
            #     filter_mask = None

            # _, recall, _ = precision_recall_curve(links, golden_links)
            # calculate_auc_roc(links, golden_links)

            # label = method + f' (score = {str(round(np.mean(recall[:-1]), 2))})'
            # def normalize(aa):
            #     return (aa
            # ax.plot(, recall, label=label, color=colors[i], alpha=1, linewidth=2,
            #             linestyle=linestyles[i])
            # print(f'{DE_type} -- {method}---{round(np.mean(recall[:-1]), 2)}')

            results = permutation_test.compute_overall_score(links)
            print(f'Results for {DE_type}, {method}: {results}')

        # ax.set_title(DE_type)
        # ax.legend(frameon=False)
        # ax.set_ylabel('Recall')
        # ax.set_xlabel('Thresholds')

    # fig_pr.savefig(os.path.join(VS_STRING_DIR, f'precision_recall.png'), dpi=300, transparent=True)
    # fig.savefig(os.path.join(VS_STRING_DIR, f'recall.pdf'))



