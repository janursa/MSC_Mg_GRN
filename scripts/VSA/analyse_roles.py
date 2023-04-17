"""
    Conduct vester's sensitivity analysis for non noised and noised links
"""
import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_model_name_2_method_and_DE_type, VSA_DIR, F_selected_models, \
    F_DE_genenames, F_protnames_to_genenames
from common_tools.VSA import role_analysis, RolePlot, determine_top_role_change, determine_critical_role_change
from VSA.vsa_aux import retreive_links_with_genenames

def custom_annotation_ctr_vs_sample(gene_name):
    offset_x, offset_y, rad = None, None, None
    if gene_name in ['AK1']:
        offset_y = 0.15
        offset_x = 0.15
        rad = -.2
    if gene_name in ['APOA2']:
        offset_x = 0.15
    if gene_name in ['OGDH']:
        offset_y = 0.2
        offset_x = -0.2
        rad = -.2

    return offset_x, offset_y, rad

def custom_annotation_role_change(prot):
    offset_x, offset_y, rad_arrow = None, None, None
    offset_text_location = False

    if prot == 'HDGFL2':
        offset_x = .25
        offset_y = +.15
        rad_arrow = +.5
        offset_text_location = True
    elif prot == 'KTN1':
        offset_x = 0.2
        offset_y = .8
        rad_arrow = -.5
        offset_text_location = True
    elif prot == 'SF3A1':
        offset_x = -.3
        offset_y = .9
        rad_arrow = +.5
        offset_text_location = True
    elif prot == 'MYL1':
        offset_y = .7
        rad_arrow = 0.5
        offset_text_location = True
    elif prot == 'KTN1':
        offset_x = .2
        offset_y = .3
        rad_arrow = -1
        offset_text_location = True
    elif prot == 'UBAP2L':
        offset_x = -.5
        offset_y = +.7
        rad_arrow = +.6
        offset_text_location = True
    # - long-term
    elif prot == 'MARS1':
        offset_x = -.55
        offset_y = -.3
        rad_arrow = +.2
    elif prot == 'AK1':
        offset_x = -.8
        offset_y = -.4
        rad_arrow = +.2
        offset_text_location = True
    elif prot == 'OGDH':
        offset_x = .2
        offset_y = -.4
        rad_arrow = +.05
        offset_text_location = True
    elif prot == 'TRIM28':
        offset_x = -.8
        # offset_y = -.15
        rad_arrow = -.05
    elif prot == 'NUCB1':
        # offset_x = -.5
        offset_y = .15
        rad_arrow = -.05
    elif prot == 'LRP1':
        offset_x = .1
        offset_y = 0.05
    elif prot == 'NUFIP2':
        offset_x = .05
        # offset_y = 0.1
        # rad = -.1

    return offset_x, offset_y, rad_arrow, offset_text_location


def plot_roles(studies_names, data_ctr, data_treatment, role_changes, critical_role_changes, gene_names):
    # - determine whether there is 3 or 4 columns depending on if critical_role_change is empty
    if len(critical_role_changes) == 0:
        ncols, nrows = 3, 1
    else:
        ncols, nrows = 4, 1
    fig_roles, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.5 * ncols, 2.5 * nrows))
    # - plot the first 2 columns, role changes for ctr and sample
    for ii, data in enumerate([data_ctr, data_treatment]):
        RolePlot.plot_ctr_vs_sample(axes[ii], data, gene_names, custom_annotation=custom_annotation_ctr_vs_sample)
        if ii == 1:
            axes[ii].set_ylabel('')
    # - define titles. by default, we have 3 columns. it can be 4 if critical_role_change is not empty
    titles = ['(A1) Control', '(A2) Treatment', '(A3) Top role change']

    # role change from ctr to sample for top role changes
    def plot_role_change(df_role_change, axes_i):
        ax_role_change = axes[axes_i]
        RolePlot.plot_role_change(df_role_change, ax=ax_role_change, custom_annotation=custom_annotation_role_change)
        ax_role_change.set_xlim([-1.4, 1.4])
        ax_role_change.set_ylim([-1.4, 1.4])
        ax_role_change.set_ylabel('')

    plot_role_change(role_changes, axes_i=2)
    # role change from ctr to sample for critical role changes
    if len(critical_role_changes) > 0:
        plot_role_change(critical_role_changes, axes_i=3)
        titles.append('(A4) Critical role change')
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
    parser.add_argument('--top_quantile_role_change', type=float, default=.9,
                        help="Top quantile of biggest distance in the role change from ctr to sample")
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    top_quantile_role_change = args.top_quantile_role_change
    # - create VSA dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)
    # - load names of selected models
    selected_models = F_selected_models()

    # - read the links and store it for selected models
    links_all: Dict[str, List[pd.DataFrame]] = retreive_links_with_genenames(selected_models,
                                                                             F_protnames_to_genenames(), studies)

    # - the main function for VSA analysis.
    for model_name in selected_models:
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        genenames = F_DE_genenames()[DE_type]
        # - VSA analysis for both studies
        vsa_results_studies = []
        links_studies = links_all[model_name]
        for study_i, _ in enumerate(studies):
            vsa_results = role_analysis(links_studies[study_i], genenames)
            vsa_results_studies.append(vsa_results)
        # - determine top role change from ctr to sample
        top_role_change = determine_top_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                    top_quantile=top_quantile_role_change)
        # - determine those proteins with a critical role change from ctr to sample
        critical_role_change = determine_critical_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                              target_role=3)
        # - write VSA analysis for ctr and sample to files
        for study_i, study in enumerate(studies):
            vsa_results_studies[study_i].to_csv(Path(VSA_DIR) / f'vsa_{model_name}_{study}.csv', index=False)
        # - write top role changes and critical role changes to file
        top_role_change.to_csv(Path(VSA_DIR) / f'top_role_change_{model_name}.csv', index=False)
        critical_role_change.to_csv(Path(VSA_DIR) / f'critical_role_change_{model_name}.csv', index=False)

        # - plot roles for ctr and sample, top role changes and critical role changes in one plot
        fig = plot_roles(studies, vsa_results_studies[0], vsa_results_studies[1], top_role_change, critical_role_change,
                         genenames)
        fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.pdf'))
