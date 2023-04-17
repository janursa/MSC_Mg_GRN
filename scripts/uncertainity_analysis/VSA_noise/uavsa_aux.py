import pandas as pd
from typing import *
import numpy as np
import matplotlib.pyplot as plt
from common_tools import serif_font
from common_tools.VSA import RolePlot
QP_series_type: TypeAlias = List[Tuple[float, float]]  # data type for list of (Q,P)

def re_organize_vsa_results(batch_results_ctr: List[pd.DataFrame],
                            batch_results_sample: List[pd.DataFrame],
                            genes_to_extract: List[str]) -> Dict[str, Tuple[QP_series_type, QP_series_type]]:
    """
        Sorts batch VSA results based on target proteins. gene: [[(Q1, P1), ...], [(Q1, P1), ..]]}, which contains repeated
        results for ctr and sample.
    """
    # - extract the data
    def extract_function(batch_results, targ_gene) -> QP_series_type:
        """Gets a list of vsa results, extracts tag values for given gene"""
        Qs_Ps = []
        for results in batch_results:
            mask = results['Entry'] == targ_gene
            Qs_Ps.append(results.loc[mask, ['Q', 'P']].values.tolist()[0])
        return Qs_Ps

    results_stack: Dict[str, Tuple[QP_series_type, QP_series_type]] = {} # gene <- (ctr, sample) <- list of Qs_Ps for each
    for gene in genes_to_extract:
        Qs_Ps_ctr = extract_function(batch_results_ctr, targ_gene=gene)
        Qs_Ps_sample = extract_function(batch_results_sample, targ_gene=gene)
        results_stack[gene] = (Qs_Ps_ctr, Qs_Ps_sample)
    return results_stack


def determine_sig_change_in_role(results: Dict[str, Tuple[QP_series_type]]) -> Dict[str, str]:
    """Determines the degree of sig change from ctr to sample by assigning
    one of *, **, ***
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

def plot_results(results, noise_type, std_noise, sig_flags):
    """ Plot role change from control to sample for n noisy results
    Scatter plots where each node shows the results of one noisy analysis.
    """
    preferred_titles = [f'{gene} {sign}' for gene, sign in sig_flags.items()]

    ncols, nrows = len(results), 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2 * ncols, 2 * nrows))
    serif_font()
    for gene_i, (gene, gen_results) in enumerate(results.items()):  # - plot for each gene on a seperate window
        ax = axes[gene_i]
        for study_i, study_data in enumerate(gen_results):
            Ps = [item[1] for item in study_data]
            Qs = [item[0] for item in study_data]
            ax.scatter(Ps, Qs,
                       color=RolePlot.ctr_sample_colors[study_i],
                       alpha=.6,
                       linewidths=.2,
                       edgecolors='black',
                       s=35
                       )
        title = preferred_titles[gene_i] if (preferred_titles is not None) else gene
        RolePlot.postprocess(ax, title=title, show_axis_names=False, xlim=[-1.2, 1.2], ylim=[-1.2, 1.2])
        RolePlot.mark_x(ax, 0)
        RolePlot.mark_y(ax, 0)
        RolePlot.plot_roles(ax)
    if noise_type == 'mp':
        axes[0].set_ylabel(f'Multiplicative noise ({int(100 * std_noise)} %)', fontweight='bold')
    if noise_type == 'ad':
        axes[0].set_ylabel(f'Additive noise ({int(100 * std_noise)} %)', fontweight='bold')
    return fig
