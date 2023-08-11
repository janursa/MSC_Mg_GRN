"""
    Conduct vester's sensitivity analysis for non noised and noised links
"""
import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import json
import matplotlib.pyplot as plt

from typing import List, Dict, Tuple, TypeAlias, Any

from proteomics_MSC.imports import VSA_ROBUSTNESS_DIR, VSA_DIR, VSA_VERBOSE_DIR
from proteomics_MSC.common_tools import F_model_name_2_method_and_DE_type, F_DE_genenames, F_selected_models
from proteomics_MSC.common_tools import serif_font, flatten
from proteomics_MSC.common_tools.role_analysis import role_analysis, RolePlot
from proteomics_MSC.common_tools.noise_analysis import NoiseAnalysis

ASPS_series_type: TypeAlias = List[Tuple[float, float]]  # data type for list of (AS,PS)

def axis_limits(model_name):
    if model_name == 'early_MinProb_Portia':
        active_sum_range = [-2, 12]
        passive_sum_range = [0, 10]
    else:
        active_sum_range = [-3, 13]
        passive_sum_range = [-1, 9]
    return active_sum_range, passive_sum_range
class VSA_noise_analysis(NoiseAnalysis):
    """Main class to conduct noise analysis for the protein role change using VSA"""
    def __init__(self, selected_genes, plot_flag=False, **kwargs):
        super().__init__(**kwargs)
        self.sig_signs = {} # keeps record of sig signs for all noise types and noise intensities
        self.selected_genes = selected_genes
        self.plot_flag = plot_flag
    def _save_file(self, noise_type, noise_intensity):
        save_dir = self._save_dir()
        save_file = f'{save_dir}/results_{noise_type}_{noise_intensity}.json'
        return save_file
    @staticmethod
    def __re_organize_vsa_results(batch_results_ctr: List[pd.DataFrame],
                                batch_results_sample: List[pd.DataFrame],
                                genes_to_extract: List[str]) -> Dict[str, Tuple[ASPS_series_type, ASPS_series_type]]:
        """
            Sorts batch VSA results based on aimed proteins. gene: [[(Q1, P1), ...], [(Q1, P1), ..]]}, which contains repeated
            results for ctr and sample.
        """

        # - extract the data
        def extract_function(batch_results, targ_gene) -> ASPS_series_type:
            """Gets a list of vsa results, extracts tag values for given gene"""
            ASs_PSs = []
            for results in batch_results:
                mask = results['Entry'] == targ_gene
                ASs_PSs.append(results.loc[mask, ['AS', 'PS']].values.tolist()[0]) # index 0 is AS
            return ASs_PSs

        results_stack: Dict[
            str, Tuple[ASPS_series_type, ASPS_series_type]] = {}  # gene <- (ctr, sample) <- list of ASs_PSs for each
        for gene in genes_to_extract:
            ASs_PSs_ctr = extract_function(batch_results_ctr, targ_gene=gene)
            ASs_PSs_sample = extract_function(batch_results_sample, targ_gene=gene)
            results_stack[gene] = (ASs_PSs_ctr, ASs_PSs_sample) # index 0 is ctr
        return results_stack

    def _run_analysis(self, noise_type:str, noise_intensity:float) -> None:
        save_file = self._save_file(noise_type, noise_intensity)
        studies = self.studies
        # - run the analysis if the output file doesnt exist
        if not os.path.exists(save_file) or self.force:
            print(f'Noise analysis for {self.model_name} {noise_type} {noise_intensity}')
            # - run grn to get the links
            noisy_links_stack_ctr, noisy_links_stack_sample = self._create_noisy_links_studies(model_name=self.model_name,
                                                                                         studies=studies,
                                                                                         noise_type=noise_type,
                                                                                         noise_intensity=noise_intensity,
                                                                                         n_repeat=self.n_repeat)
            # - run vsa
            method, DE_type = F_model_name_2_method_and_DE_type(self.model_name)
            genenames = F_DE_genenames()[DE_type]
            vsa_results_stack_ctr = [role_analysis(links, genenames) for links in noisy_links_stack_ctr]
            vsa_results_stack_sample = [role_analysis(links, genenames) for links in noisy_links_stack_sample]
            # - re structred the data to have: gene: {ctr:[(Q1, P1), ...], sample:[(Q1, P1), ..]}
            gene_based_vsa_results = self.__re_organize_vsa_results(vsa_results_stack_ctr, vsa_results_stack_sample,
                                                                    genes_to_extract=self.selected_genes)
            with open(save_file, 'w') as ff:
                print(f'output -> {save_file}')
                json.dump(gene_based_vsa_results, ff)
    def _run_post_analysis(self,noise_type:str, noise_intensity:float) -> None:
        """The analysis of the results and plots"""
        # - retrieve the results
        save_file = self._save_file(noise_type, noise_intensity)
        with open(save_file, 'r') as ff:
            results = json.load(ff)
        # - determine sig level from ctr to sample: *, **, ***
        sig_signs_noisetype_intensity = self.__determine_sig_change(results)
        # - store the sig level for robustness check
        if noise_type not in self.sig_signs.keys():
            self.sig_signs[noise_type] = {}
        self.sig_signs[noise_type][noise_intensity] = sig_signs_noisetype_intensity
        # - plot scatters showing ctr to sample change
        if self.plot_flag:
            save_plots_dir = self._save_plot_dir()
            self.__plot_results(results, noise_type, noise_intensity, sig_signs_noisetype_intensity, save_plots_dir)
    def determine_robust_targs(self) -> List[str]:
        """Determines robust changes by taking into account sig changes across different noise types"""
        robust_targs = []
        for noise_type, noise_type_data in self.sig_signs.items():
            threshold = self.robust_thresholds[noise_type]  # the noise intensity that we consider for robustness
            target_results = noise_type_data[threshold]
            robust_genes_noisetype = [gene for gene, sig_sig in target_results.items() if (sig_sig != '')]
            robust_targs.append(robust_genes_noisetype)
        list_of_sets = [set(sublist) for sublist in robust_targs]
        # Find the intersection of all sets
        robust_targs = list(set.intersection(*list_of_sets))
        return robust_targs
    @staticmethod
    def __plot_results(results, noise_type, std_noise, sig_flags, save_plots_dir):
        """ Plot role change from control to sample for n noisy results
        Scatter plots where each node shows the results of one noisy analysis.
        """

        preferred_titles = [f'{gene} {sign}' for gene, sign in sig_flags.items()]

        ncols, nrows = len(results), 1
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2 * ncols, 2 * nrows))
        serif_font()
        for gene_i, (gene, gene_results) in enumerate(results.items()):  # - plot for each gene on a seperate window
            ax = axes[gene_i]
            for study_i, study_data in enumerate(gene_results):
                PSs = [item[1] for item in study_data]
                ASs = [item[0] for item in study_data]
                ax.scatter(PSs, ASs,
                           color=RolePlot.ctr_sample_colors[study_i],
                           alpha=.6,
                           linewidths=.2,
                           edgecolors='black',
                           s=35
                           )
            title = preferred_titles[gene_i] if (preferred_titles is not None) else gene
            xlim, ylim = axis_limits(model_name)
            RolePlot.postprocess(ax, title=title, show_axis_names=True, xlim=xlim, ylim=ylim)
            # RolePlot.mark_x(ax, 0)
            # RolePlot.mark_y(ax, 0)
            RolePlot.plot_roles(ax)

        fig.savefig(Path(save_plots_dir) / f'{noise_type}_{std_noise}.pdf')
        fig.savefig(Path(save_plots_dir) / f'{noise_type}_{std_noise}.png', dpi=150, transparent=True)
    @staticmethod
    def __determine_sig_change(results: Dict[str, Tuple[List[float]]]) -> Dict[str, str]:
        """Determines the degree of sig change from ctr to sample by assigning
        one of *, **, ***.
        Results: ctr and sample data for each gene
        """
        sig_signs = {}
        for gene, gene_results in results.items():
            ctr, sample = gene_results
            # - runs the test
            Qs_ctr = [point[0] for point in ctr]
            Ps_ctr = [point[1] for point in ctr]
            Qs_sample = [point[0] for point in sample]
            Ps_sample = [point[1] for point in sample]

            abs_diff_Qs = np.abs(np.mean(Qs_ctr) - np.mean(Qs_sample))
            abs_diff_Ps = np.abs(np.mean(Ps_ctr) - np.mean(Ps_sample))
            Q_std = np.std(Qs_ctr)
            P_std = np.std(Ps_ctr)
            sign = ''
            if (abs_diff_Qs > Q_std) | (abs_diff_Ps > P_std):
                if (abs_diff_Qs > 2 * Q_std) | (abs_diff_Ps > 2 * P_std):
                    if (abs_diff_Qs > 3 * Q_std) | (abs_diff_Ps > 3 * P_std):
                        sign = '***'
                    else:
                        sign = '**'
                else:
                    sign = '*'
            sig_signs[gene] = sign
        return sig_signs

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    # - parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--force', type=bool, default=False,
                        help='To overwrite if the result files exist')
    args, remaining_args = parser.parse_known_args()

    studies = config['studies']
    force = args.force
    n_repeat = config['n_repeat_noise_analysis']

    # - run for each model
    for model_name in F_selected_models():
        selected_genes = np.loadtxt(f'{VSA_VERBOSE_DIR}/top_genes_{model_name}.csv', dtype=str)
        VSA_noise_analysis_obj = VSA_noise_analysis(model_name=model_name, studies=studies, force=force, base_dir=VSA_ROBUSTNESS_DIR, n_repeat=n_repeat, selected_genes=selected_genes, plot_flag=False)
        VSA_noise_analysis_obj.run_decorator()
        # - robust targs
        robust_targs = VSA_noise_analysis_obj.determine_robust_targs()
        to_save = f'{VSA_DIR}/robust_genes_{model_name}.csv'
        np.savetxt(to_save, robust_targs, fmt='%s')
        print(f'output -> {to_save}')

