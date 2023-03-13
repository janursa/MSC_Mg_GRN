"""
    Conduct vester's sensitivity analysis for non noised links, and plot VSA results for ctr and mg and the role change
"""
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import F_DE_protiens, VSA_DIR
from scripts.utils.VSA import role_analysis, RolePlot
from analyse import retreive_links_shortlisted
from scripts.utils import MG_noise_F, AG_noise_F, make_title_pretty
def filter_and_sort(batch_results_1, batch_results_2, target_genes):
    """
        Sorts batch VSA results based on target proteins
    """
    sorted_VSA = {}
    def extract_F(batch_results, prot, tag='Q'):
        batch_results_filtered = [oo.loc[oo['Entry']==prot,:] for oo in batch_results]
        batch_results_ = [oo[tag].values.tolist()[0] for oo in batch_results_filtered]
        return batch_results_

    for prot in target_genes:
        Qs_ctr = extract_F(batch_results_1, prot=prot, tag='Q')
        Ps_ctr = extract_F(batch_results_1, prot=prot, tag='P')
        Roles_ctr = extract_F(batch_results_1, prot=prot, tag='Role')

        Qs_sample = extract_F(batch_results_2, prot=prot, tag='Q')
        Ps_sample = extract_F(batch_results_2, prot=prot, tag='P')
        Roles_sample = extract_F(batch_results_2, prot=prot, tag='Role')

        sorted_VSA[prot] = {'ctr': {'Q': Qs_ctr, 'P': Ps_ctr, 'Role': Roles_ctr},
                            'Mg': {'Q': Qs_sample, 'P': Ps_sample, 'Role': Roles_sample}}
    return sorted_VSA

def batch_VSA(links_batch, gene_names, target_genes):
    """
        Runs VSA for the entire noised dfs, shortlists the target genes
    """
    oo_batch = []
    for links in links_batch:
        rr = role_analysis(links, gene_names)
        rr_short = rr.loc[rr['Entry'].isin(target_genes), :]
        oo_batch.append(rr_short)
    return oo_batch

if __name__ == '__main__':
    # - create dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)
    std_MF = .3
    rel_std_AF = .3
    methods = ['RF','ridge','portia']
    studies = ['ctr', 'mg']
    selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    DIR_NOISE = os.path.join(VSA_DIR, 'noise')
    if not os.path.isdir(DIR_NOISE):
        os.makedirs(DIR_NOISE)

    for DE_type, DE_proteins in F_DE_protiens().items():
        for method in methods:
            model_name = '_'.join([DE_type, method])
            if model_name not in selected_models:  # only selected models
                continue
            """
                Multi step function to conduct multiplicative gaussian noise analysis
            """
            # - retreive data
            links_stack = retreive_links_shortlisted(method, DE_type)
            role_change = pd.read_csv(os.path.join(VSA_DIR, f'role_change_{model_name}.csv'), index_col=False)
            target_genes = role_change['Entry'].values.tolist()
            print('Number of target genes ', len(target_genes))
            # - create noised dfs
            links_ctr_noised = MG_noise_F(links_stack[0], std=std_MF) # 100 noised links
            links_sample_noised = MG_noise_F(links_stack[1],std=std_MF)
            # - run VSA for noised dfs
            batch_results_ctr = batch_VSA(links_ctr_noised, gene_names=DE_proteins, target_genes=target_genes)
            batch_results_sample = batch_VSA(links_sample_noised, gene_names=DE_proteins, target_genes=target_genes)
            # - we plot noisy role change for each protein: filter batch results for target proteins
            rr_sorted_M = filter_and_sort(batch_results_ctr, batch_results_sample, target_genes=target_genes)

            """
                   Multi step function to conduct additive gaussian noise analysis
            """
            # - create noised dfs
            links_ctr_noised = AG_noise_F(links_stack[0], rel_std=rel_std_AF)
            links_sample_noised = AG_noise_F(links_stack[1], rel_std=rel_std_AF)
            # - run VSA for noised dfs
            batch_results_ctr = batch_VSA(links_ctr_noised, gene_names=DE_proteins, target_genes=target_genes)
            batch_results_sample = batch_VSA(links_sample_noised, gene_names=DE_proteins, target_genes=target_genes)
            # - filter series for target proteins and sort them as a dict
            rr_sorted_A = filter_and_sort(batch_results_ctr, batch_results_sample, target_genes=target_genes)


            # - plot
            ncols = len(target_genes)
            nrows = 2
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2 * ncols, 2 * nrows))
            RolePlot.plot_noise_analysis(axes[0][:], rr_sorted_M, study_1='ctr', study_2='Mg', show_title=True)
            axes[0][0].set_ylabel('Multiplicative noise', fontweight='bold')
            RolePlot.plot_noise_analysis(axes[1][:], rr_sorted_A, study_1='ctr', study_2='Mg', show_title=False)
            axes[1][0].set_ylabel('Additive noise', fontweight='bold')


            # fig.suptitle(make_title_pretty(model_name), fontweight='bold')

            fig.savefig(os.path.join(DIR_NOISE, f'noiseAnalysis_{model_name}.pdf'))
            fig.savefig(os.path.join(DIR_NOISE, f'noiseAnalysis_{model_name}.png'), dpi=300, transparent=True)


