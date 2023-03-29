"""
    Conduct vester's sensitivity analysis for non noised and noised links
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import logging as lg
import numpy as np
from typing import List
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, VSA_DIR, F_DE_protnames, MODELSELECTION_DIR
from scripts.utils.VSA import role_analysis, RolePlot, role_change
from scripts.utils import make_title_pretty, comic_font
from scripts.utils.links import choose_top_quantile
from scripts.utils import MG_noise_F, AG_noise_F
from scripts.utils import serif_font


def batch_VSA(links_batch: List[pd.DataFrame], gene_names: List[str], target_genes: List[str]) -> List[pd.DataFrame]:
    results_batch = [role_analysis(links, gene_names) for links in links_batch]
    results_batch = [oo.loc[oo['Entry'].isin(target_genes), :] for oo in results_batch]
    return results_batch
def filter_and_sort(batch_results_1, batch_results_2, target_genes):
    """
        Sorts batch VSA results based on target proteins
    """
    sorted_VSA = {}
    def extract_F(batch_results, prot, tag='Q'):
        batch_results_filtered = [oo.query(f"Entry == '{prot}'") for oo in batch_results]
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
if __name__ == '__main__':
    lg.basicConfig(filename=os.path.join(VSA_DIR, 'log.log'), filemode='w', format='%(name)s - %(levelname)s - %(message)s')
    #- create dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)
    studies = ['ctr', 'mg']
    methods = ['ensemble']
    std_MF = .3 #multiplicative noise std
    rel_std_AF = .3
    n_relica_noise =100

    model_i = 0
    for idx, (DE_type, DE_proteins) in enumerate(F_DE_protnames().items()):
        for method in methods:
            links_stack = []
            for study in studies:
                links_path = Path(GRN_DIR) / method / f'links_{study}.csv'
                links = pd.read_csv(links_path, index_col=False)
                links = choose_top_quantile(links, 0.5)
                links_stack.append(links)
            oo_stack = [role_analysis(links, DE_proteins) for links in links_stack]
            df_role_change = role_change(oo_stack[0], oo_stack[1], target_role=3)
            df_role_change.to_csv(os.path.join(VSA_DIR, f'role_change.csv'), index=False)
            # -------- plot --------------
            # ctr vs sample
            ncols, nrows = 3, 1
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(3 * ncols, 3 * nrows))
            serif_font()
            # fig.suptitle(t='(A) Study of protein roles using VSA', x=.5, y=.93, fontweight='bold', fontsize=12)
            axes_ctr_mg = axes[0:2]
            RolePlot.plot_ctr_vs_sample(axes_ctr_mg, oo_stack, studies, DE_proteins)
            axes_ctr_mg[1].set_ylabel('')
            # role change from ctr to sample
            ax_role_change = axes[2]
            RolePlot.plot_role_change(df_role_change, ax=ax_role_change)
            ax_role_change.set_ymargin(.4)
            ax_role_change.set_xmargin(.4)
            ax_role_change.set_ylabel('')
            # legends
            handles = []
            for i, color in enumerate(RolePlot.ctr_sample_colors):
                handles.append(
                    ax_role_change.scatter([], [], marker='o', label=studies[i], s=100, edgecolor='black', color=color, linewidth=.2))
            handles.append(ax_role_change.scatter([], [], marker='', label='--', s=1,
                                                  edgecolor='white', color='white',
                                                  linewidth=.5))  # create a blank place

            for i, color in enumerate(RolePlot.roles_colors):
                handles.append(ax_role_change.scatter([], [], marker='o', label=RolePlot.roles[i], s=100,
                                          edgecolor='black', color=color, linewidth=.5))
            ax_role_change.legend(handles=handles, loc='upper center',
                            bbox_to_anchor=(1.3, 1), prop={'size': 10}, frameon=False
                            )
            fig.savefig(os.path.join(VSA_DIR, f'role_changes.png'), dpi=300, transparent=True)
            fig.savefig(os.path.join(VSA_DIR, f'role_changes.pdf'))
            # ------------------- noise analysis
            target_genes = df_role_change['Entry'].values.tolist()
            lg.info(f'Number of target genes {len(target_genes)}')

            # - multiplicative noise
            links_ctr_noised = MG_noise_F(links_stack[0], n_relica=n_relica_noise, std=std_MF)  # 100 noised links
            links_sample_noised = MG_noise_F(links_stack[1], n_relica=n_relica_noise, std=std_MF)
            # run VSA for noised dfs
            batch_results_ctr = batch_VSA(links_ctr_noised, gene_names=DE_proteins, target_genes=target_genes)
            batch_results_sample = batch_VSA(links_sample_noised, gene_names=DE_proteins, target_genes=target_genes)
            # keep only target proteins and sort based on P, Q, and role
            rr_sorted_M = filter_and_sort(batch_results_ctr, batch_results_sample, target_genes=target_genes)
            lg.info('multipliciative noise is calculated')
            # - additive noise
            links_ctr_noised = AG_noise_F(links_stack[0], n_relica=n_relica_noise, rel_std=rel_std_AF)
            links_sample_noised = AG_noise_F(links_stack[1], n_relica=n_relica_noise, rel_std=rel_std_AF)
            # run VSA for noised dfs
            batch_results_ctr = batch_VSA(links_ctr_noised, gene_names=DE_proteins, target_genes=target_genes)
            batch_results_sample = batch_VSA(links_sample_noised, gene_names=DE_proteins, target_genes=target_genes)
            # keep only target proteins and sort based on P, Q, and role
            rr_sorted_A = filter_and_sort(batch_results_ctr, batch_results_sample, target_genes=target_genes)
            lg.info('additive noise is calculated')
            #- plot
            ncols, nrows = len(target_genes), 2
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2 * ncols, 2 * nrows))
            # fig.suptitle(t='(B) Robustness of protein role change to noise', x=.5, y=.95, fontweight='bold', fontsize=11)
            RolePlot.plot_noise_analysis(axes[0][:], rr_sorted_M, study_1='ctr', study_2='Mg', show_title=True)
            axes[0][0].set_ylabel('Multiplicative noise', fontweight='bold')
            RolePlot.plot_noise_analysis(axes[1][:], rr_sorted_A, study_1='ctr', study_2='Mg', show_title=False)
            axes[1][0].set_ylabel('Additive noise', fontweight='bold')

            fig.savefig(os.path.join(VSA_DIR, 'noise', f'noiseAnalysis.pdf'))
            fig.savefig(os.path.join(VSA_DIR, 'noise', f'noiseAnalysis.png'), dpi=300, transparent=True)

