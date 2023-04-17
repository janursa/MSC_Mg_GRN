"""
    Sets of functions useful for uncertainity_analysis analysis against noise
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple, Callable
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from common_tools import serif_font
from common_tools.VSA import role_analysis, RolePlot


class NoiseAnalysis:
    """
    This class handles the functions for the uncertainity_analysis analysis of VSA results for role change from ctr to
    sample. It introduces different types and levels of noise to input datasets and runs the GRN methods to infer
    noise outputs. Then, it uses VSA to determine the roles for the noisy datasets. Finally, it plots the noisy
    roles to examine if the roles are statistically different from ctr to sample.
    """

    def __init__(self, data_ctr: np.ndarray, data_sample: np.ndarray, target_genes: List[str], run_grn_func: Callable,
                 n_noisy_datasets, std_noise, noise_type,
                 kwargs_grn_func, kwargs_role_analysis):
        assert noise_type in ['mp', 'ad']
        self.data_ctr = data_ctr
        self.data_sample = data_sample
        self.target_genes = target_genes
        self.run_grn_func = run_grn_func
        self.kwargs_grn_func = kwargs_grn_func
        self.kwargs_role_analysis = kwargs_role_analysis
        self.studies = ['ctr', 'sample']
        self.n_noisy_datasets = n_noisy_datasets
        self.std_noise = std_noise
        self.noise_type = noise_type

    def analyse_noise(self):
        """ Runs the model for the stack of noisy datasets
        Add noise to the data and create a stack. Run GRN for each data in the stack.
        """
        vsa_results_stack_studies = [] # results for both control and sample
        for i_data, data in enumerate([self.data_ctr, self.data_sample]):
            # - create noisy data stack
            if self.noise_type == 'mp':
                noisy_data_stack = NoiseAnalysis._add_mpnoise(data, n_noisy_datasets=self.n_noisy_datasets,
                                                              std=self.std_noise)
            else:
                noisy_data_stack = NoiseAnalysis._add_adnoise(data, n_noisy_datasets=self.n_noisy_datasets,
                                                              std=self.std_noise)
            # - run GRN for noisy data
            noisy_links_stack = []
            with tqdm(total=self.n_noisy_datasets,
                      desc=f'Run GRN for noisy data: {self.noise_type} {self.std_noise} study {i_data}') as pbar:
                for noisy_data in noisy_data_stack:
                    noisy_links_stack.append(self.run_grn_func(data=noisy_data, i_data=i_data, **self.kwargs_grn_func))
                    pbar.update(1)
            # - VSA analysis on the noisy data
            vsa_results_stack = [role_analysis(links, **self.kwargs_role_analysis) for links in noisy_links_stack]
            # - filter VSA analysis for the target genes
            vsa_results_stack_filtered = [oo.loc[oo['Entry'].isin(self.target_genes), :] for oo in vsa_results_stack]
            vsa_results_stack_studies.append(vsa_results_stack_filtered)
        # - re structred the data to make ctr and sample data belong to the same target proteins
        results = self._organized(vsa_results_stack_studies[0], vsa_results_stack_studies[1])
        return results

    @staticmethod
    def _plot(axes, data, genenames, preferred_titles=None):
        """
            Plots Q/P for ctr and mg conditions, for n different runs of sensitivity analysis, one window for each prot

            oo_prots: Q/P data for ctr and sample, for each gene
        """
        assert len(data) == len(genenames)
        serif_font()

        for idx, (gene, data) in enumerate(zip(genenames, data)):  # - plot for each prot on a seperate window
            ax = axes[idx]
            for idxx, (study, study_data) in enumerate(data.items()):
                Ps = [item[1] for item in study_data]
                Qs = [item[0] for item in study_data]
                ax.scatter(Ps, Qs,
                           color=RolePlot.ctr_sample_colors[idxx],
                           alpha=.6,
                           linewidths=.2,
                           edgecolors='black',
                           s=35
                           )
            title = preferred_titles[idx] if (preferred_titles is not None) else gene
            RolePlot.postprocess(ax, title=title, show_axis_names=False)
            RolePlot.mark_x(ax, 0)
            RolePlot.mark_y(ax, 0)
            RolePlot.plot_roles(ax)

    @staticmethod
    def _role_change_test(ctr, sample):
        Qs_ctr = [point[0] for point in ctr]
        Ps_ctr = [point[1] for point in ctr]
        Qs_sample = [point[0] for point in sample]
        Ps_sample = [point[1] for point in sample]

        abs_diff_Qs = np.abs(np.mean(Qs_ctr) - np.mean(Qs_sample))
        abs_diff_Ps = np.abs(np.mean(Ps_ctr) - np.mean(Ps_sample))
        Q_std = np.std(Qs_ctr)
        P_std = np.std(Ps_ctr)

        if (abs_diff_Qs > Q_std) | (abs_diff_Ps > P_std):
            if (abs_diff_Qs > 2 * Q_std) | (abs_diff_Ps > 2 * P_std):
                if (abs_diff_Qs > 3 * Q_std) | (abs_diff_Ps > 3 * P_std):
                    return '***'
                else:
                    return '**'
            else:
                return '*'
        else:
            return ''

    @staticmethod
    def determine_sig_change_in_role(results, target_genes):
        """Determines the degree of sig change from ctr to sample by assigning
        one of *, **, ***
        """
        assert len(results) == len(target_genes)
        sig_signs = []
        for gene_i, gene in enumerate(target_genes):
            ctr = results[gene_i]['ctr']
            sample = results[gene_i]['sample']
            sig_signs.append(NoiseAnalysis._role_change_test(ctr, sample))
        return sig_signs

    @staticmethod
    def plot_results(results, target_genes, noise_type, std_noise, sig_flags):
        """ Plot role change from control to sample for n noisy results
        Scatter plots where each node shows the results of one noisy analysis.
        The nodes for ctr and sample are plotted together with different colors.
        Sig. changes are marked with *, **, ***
        """
        preferred_titles = [f'{gene} {sign}' for gene, sign in zip(target_genes, sig_flags)]

        ncols, nrows = len(target_genes), 1
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2 * ncols, 2 * nrows))

        ax = axes
        NoiseAnalysis._plot(ax, results, target_genes, preferred_titles=preferred_titles)
        if noise_type == 'mp':
            ax[0].set_ylabel(f'Multiplicative noise ({int(100 * std_noise)} %)', fontweight='bold')
        if noise_type == 'ad':
            ax[0].set_ylabel(f'Additive noise ({int(100 * std_noise)} %)', fontweight='bold')
        return fig

    @staticmethod
    def _add_mpnoise(data, n_noisy_datasets=100, std=.3) -> List[np.array]:
        """ Multiplicative noise
             Creates n_noisy_datasets data
        """
        data_std = np.std(data)
        applied_std = std * data_std
        noisy_data_stack = []
        for i in range(n_noisy_datasets):
            frac_noise = np.random.rand()
            noise_mask = np.random.choice([0, 1], size=data.shape, p=[1 - frac_noise, frac_noise])
            rand_values = np.random.normal(loc=1, scale=applied_std, size=data.shape)
            noisy_data = data * (1 + (rand_values - 1) * noise_mask)
            noisy_data_stack.append(noisy_data)

        return noisy_data_stack

    @staticmethod
    def _add_adnoise(data, n_noisy_datasets=100, std=.05) -> List[np.array]:
        """ Additive noise
                 Creates n_relica noised links
        """
        data_std = np.std(data)
        applied_std = data_std * std
        noisy_data_stack = []
        for i in range(n_noisy_datasets):
            rand_values = np.random.normal(loc=1, scale=applied_std, size=data.shape)
            noisy_data = rand_values + data
            noisy_data_stack.append(noisy_data)

        return noisy_data_stack

    def _organized(self, batch_results_ctr, batch_results_sample):
        """
            Sorts batch VSA results based on target proteins
        """
        batch_results_studies = [batch_results_ctr, batch_results_sample]

        def extract_F(batch_results, prot_name, tag='Q'):
            # - extract results for prot
            batch_results_filtered = [results.query(f"Entry == '{prot_name}'") for results in batch_results]
            batch_results_ = np.asarray([results[tag].to_numpy()[0] for results in batch_results_filtered])
            return batch_results_

        results_stack = []
        for prot in self.target_genes:
            prot_results = {}
            for i_study, study in enumerate(self.studies):
                Qs = extract_F(batch_results_studies[i_study], prot_name=prot, tag='Q')
                Ps = extract_F(batch_results_studies[i_study], prot_name=prot, tag='P')
                prot_results[study] = [(Q, P) for Q, P in zip(Qs, Ps)]
            results_stack.append(prot_results)
        return results_stack
