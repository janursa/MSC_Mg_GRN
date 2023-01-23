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

from imports import GRN_DIR, protnames, VSA_DIR
from utils.links import read_write_links
from utils.VSA import VestersSA, VSA_plot, role_change
from step1 import retreive_grn
from utils import MG_noise_F, AG_noise_F
def filter_and_sort(batch_results_1, batch_results_2, target_genes):
    """
        Sorts batch VSA results based on target proteins
    """
    sorted_VSA = {}
    def extract_F(batch_results, prot, tag='AS'):
        batch_results_filtered = [oo.loc[oo['Entry']==prot,:] for oo in batch_results]
        batch_results_ = [oo[tag].values.tolist()[0] for oo in batch_results_filtered]
        return batch_results_

    for prot in target_genes:
        ASs_ctr = extract_F(batch_results_1, prot=prot, tag='AS')
        PSs_ctr = extract_F(batch_results_1, prot=prot, tag='PS')
        Roles_ctr = extract_F(batch_results_1, prot=prot, tag='Role')

        ASs_mg = extract_F(batch_results_2, prot=prot, tag='AS')
        PSs_mg = extract_F(batch_results_2, prot=prot, tag='PS')
        Roles_mg = extract_F(batch_results_2, prot=prot, tag='Role')

        sorted_VSA[prot] = {'ctr': {'AS': ASs_ctr, 'PS': PSs_ctr, 'Role': Roles_ctr},
                            'Mg': {'AS': ASs_mg, 'PS': PSs_mg, 'Role': Roles_mg}}
    return sorted_VSA

def batch_VSA(links_batch, gene_names, target_genes):
    """
        Runs VSA for the entire noised dfs, shortlists the target genes
    """
    oo_batch = []
    for links in links_batch:
        rr = VestersSA(links, gene_names)
        rr_short = rr.loc[rr['Entry'].isin(target_genes), :]
        oo_batch.append(rr_short)
    return oo_batch

if __name__ == '__main__':
    method = 'portia'
    """
        Multi step function to conduct multiplicative gaussian noise analysis
    """
    # - retreive data
    links_ctr, links_sample = retreive_grn(method)
    role_change = pd.read_csv(os.path.join(VSA_DIR, f'role_change_{method}.csv'), index_col=False)
    target_genes = role_change['Entry'].values.tolist()
    # - create noised dfs
    links_ctr_noised = MG_noise_F(links_ctr) # 100 noised links
    links_sample_noised = MG_noise_F(links_sample)
    # - run VSA for noised dfs
    batch_results_ctr = batch_VSA(links_ctr_noised, gene_names=protnames, target_genes=target_genes)
    batch_results_sample = batch_VSA(links_sample_noised, gene_names=protnames, target_genes=target_genes)
    # - we plot noisy role change for each protein: filter batch results for target proteins
    rr_sorted = filter_and_sort(batch_results_ctr, batch_results_sample, target_genes=target_genes)
    # - plot
    fig = VSA_plot.plot_noise_analysis(rr_sorted, study_1='ctr', study_2='Mg')
    fig.savefig(os.path.join(VSA_DIR, f'multiplicative_{method}.pdf'))

    """
           Multi step function to conduct additive gaussian noise analysis
    """
    # - create noised dfs
    links_ctr_noised = AG_noise_F(links_ctr)
    links_sample_noised = AG_noise_F(links_sample)
    # - run VSA for noised dfs
    batch_results_ctr = batch_VSA(links_ctr_noised, gene_names=protnames, target_genes=target_genes)
    batch_results_sample = batch_VSA(links_sample_noised, gene_names=protnames, target_genes=target_genes)
    # - filter series for target proteins and sort them as a dict
    rr_sorted = filter_and_sort(batch_results_ctr, batch_results_sample, target_genes=target_genes)
    # - plot
    fig = VSA_plot.plot_noise_analysis(rr_sorted, study_1='ctr', study_2='Mg')
    fig.savefig(os.path.join(VSA_DIR, f'addGaussNoise{method}.pdf'))


