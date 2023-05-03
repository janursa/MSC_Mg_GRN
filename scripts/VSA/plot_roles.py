import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List
import numpy as np

import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_model_name_2_method_and_DE_type, VSA_DIR, F_selected_models, \
    F_DE_genenames, F_protnames_to_genenames, VSA_NOISE_DIR
from common_tools.role_analysis import  RolePlot
from uncertainity_analysis.ua_aux import NoiseAnalysis

def custom_annotation_ctr_vs_sample(gene_name):
    offset_x, offset_y = None, None
    if model_name == 'early_MinProb_portia':
        pass
    else:
        if gene_name in ['APOA2']:
            offset_y = .1
            offset_x = -0.05
        elif gene_name in ['TRIM28']:
            offset_y = 0.1
            offset_x = -0.4
    return offset_x, offset_y

def custom_annotation_role_change(prot):
    offset_x, offset_y, rad_arrow = None, None, None
    if model_name == 'early_MinProb_portia':
        # - short-term\

        if prot == 'MYL1':
            offset_x = .3
            offset_y = -.2
            rad_arrow = 0.1
        elif prot == 'SF3A1':
            offset_x = -.1
            offset_y = .2
            rad_arrow = -0.3
        elif prot == 'RPL17':
            offset_x = -.3
            offset_y = .2
            rad_arrow = -.3
        elif prot == 'EIF3E':
            offset_x = 0
            offset_y = -.3
            rad_arrow = -.2

    else: # - long-term
        if prot == 'GLS':
            offset_x = -.2
            offset_y = .1
            rad_arrow = -0.3
        elif prot == 'MYL1':
            offset_x = .3
            offset_y = -.2
            rad_arrow = 0.1
        elif prot == 'XRCC6':
            offset_x = +.4
            offset_y = -.1
            rad_arrow = -.3
        elif prot == 'TRIM28':
            offset_x = -.2
            offset_y = -.15
            rad_arrow = -.5
        elif prot == 'MAP4':
            offset_x = +.2
            offset_y = -0.2
            rad_arrow = .3

    return offset_x, offset_y, rad_arrow

def axis_limits(model_name):
    if model_name == 'early_MinProb_Portia':
        active_sum_range = [-2, 12]
        passive_sum_range = [0, 10]
    else:
        active_sum_range = [-3, 13]
        passive_sum_range = [-1, 9]
    return active_sum_range, passive_sum_range


def plot_roles(studies_names, data_ctr, data_treatment, top_role_changes, gene_names):
    # - define titles. by default, we have 3 columns. it can be 4 if critical_role_change is not empty
    titles = ['(1) Control', '(2) Treatment', '(3) Significant role change']
    active_sum_range, passive_sum_range = axis_limits(model_name)
    # - rows and columns
    ncols, nrows = len(titles), 1
    fig_roles, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.5 * ncols, 2.5 * nrows))
    # - plot the first 2 columns, role changes for ctr and sample
    for ii, data in enumerate([data_ctr, data_treatment]):
        RolePlot.plot_ctr_vs_sample(axes[ii], data, gene_names, custom_annotation=custom_annotation_ctr_vs_sample, active_sum_range=active_sum_range, passive_sum_range=passive_sum_range)
        if ii == 1:
            axes[ii].set_ylabel('')
    # - 3rd col, top role change
    n_top_role_change = len(top_role_changes)
    # - top role change a
    def sub_plot(ax_i, i_start, i_end):
        top_role_changes_a = top_role_changes[i_start: i_end]
        RolePlot.plot_role_change(top_role_changes_a, ax=axes[ax_i], custom_annotation=custom_annotation_role_change,
                                  active_sum_range=active_sum_range, passive_sum_range=passive_sum_range)
        axes[ax_i].set_ylabel('')
    sub_plot(2, 0, n_top_role_change)

    # - set the titles
    for ii, ax in enumerate(axes):
        ax.set_title(titles[ii], fontweight='bold', fontsize=10)
    # legends. ctr and mg as well as role names
    handles = []
    for i, color in enumerate(RolePlot.ctr_sample_colors):
        handles.append(
            axes[-1].scatter([], [], marker='o', label=studies_names[i], s=100, edgecolor='black', color=color,
                             linewidth=.2))
    handles.append(axes[-1].scatter([], [], marker='', label='--', s=1,
                                    edgecolor='white', color='white',
                                    linewidth=.5))  # create a blank place

    handles.extend(RolePlot.create_role_legends(axes[-1]))
    axes[-1].legend(handles=handles, loc='upper center',
                    bbox_to_anchor=(1.3, 1), prop={'size': 10}, frameon=False
                    )
    return fig_roles


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    # - create VSA dir
    VSA_PLOT_DIR = Path(VSA_DIR) / 'plots'
    if not os.path.isdir(VSA_PLOT_DIR):
        os.makedirs(VSA_PLOT_DIR)
    # - the main function for VSA analysis.
    for model_name in F_selected_models():
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        genenames = F_DE_genenames()[DE_type]
        # - get the data
        vsa_results_studies_ctr = pd.read_csv(Path(VSA_DIR) / f'vsa_{model_name}_{studies[0]}.csv', index_col=False)
        vsa_results_studies_sample = pd.read_csv(Path(VSA_DIR) / f'vsa_{model_name}_{studies[1]}.csv', index_col=False)
        top_role_change = pd.read_csv(Path(VSA_DIR) / f'top_role_change_{model_name}.csv', index_col=False)
        # - robust role change genes
        with open(f"{VSA_NOISE_DIR}/target_genes_{model_name}.txt", 'r') as file:
            target_genes = list(np.loadtxt(file, dtype=str, delimiter=','))
        top_role_change_shortlist = top_role_change.loc[top_role_change['Entry'].isin(target_genes), :].reset_index()
        # - plot roles for ctr and sample, top role changes and critical role changes in one plot
        fig = plot_roles(studies, vsa_results_studies_ctr, vsa_results_studies_sample, top_role_change_shortlist,
                         genenames)
        fig.savefig(os.path.join(VSA_PLOT_DIR, f'role_changes_{model_name}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(VSA_PLOT_DIR, f'role_changes_{model_name}.pdf'))

