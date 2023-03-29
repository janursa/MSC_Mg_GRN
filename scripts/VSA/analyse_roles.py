"""
    Conduct vester's sensitivity analysis for non noised and noised links
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import numpy as np
import json
from typing import Dict, List

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import GRN_DIR, VSA_DIR, F_DE_protnames, MODELSELECTION_DIR, F_DE_genenames, F_protnames_to_genenames
from utils.VSA import role_analysis, RolePlot, determine_top_role_change, determine_critical_role_change
from utils import serif_font



def write_to_file(df, file_name):
    with open(file_name, 'w') as f:
        df.to_csv(f, index=False)
def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data
def custom_annotation_role_change(prot):
    offset_x, offset_y, arrow_t, rad = None, None, None, None
    if prot == 'Q02790':
        offset_x = .1
        offset_y = 0.1
        # rad = -.1
        # arrow_t = None
    elif prot == 'Q07866':
        offset_x = -0.45
        offset_y = 0.1
        rad = -.1
    elif prot == 'Q99613':
        offset_x = -0.3
        offset_y = -.2
        rad = +.2
    elif prot == 'Q9UKY7':
        offset_x = -.1
        offset_y = -.2
        rad = +.2
    elif prot == 'P67936':
        offset_x = -.43
        offset_y = .07
        rad = +.1
    elif prot == 'Q9HBL0':
        offset_x = 0.08
        offset_y = -.03
        rad = +.2
    elif (prot == 'P56192') & (model_name=='day1_1_KNN_RF'):
        rad = +.1
    elif (prot == 'P56192') & (model_name=='day1_21_KNN_portia'):
        offset_x = -.25
        offset_y = -.2
        rad = +.2
    elif prot == 'P00568':
        offset_x = -.5
        offset_y = -.15
        rad = +.2
    elif prot == 'Q02218':
        offset_x = .1
        offset_y = -.15
        rad = +.05
    elif prot == 'Q13263':
        offset_x = -.5
        # offset_y = -.15
        rad = -.05
    elif prot == 'Q02818':
        # offset_x = -.5
        # offset_y = -.15
        rad = -.05
    elif prot == 'Q07954':
        offset_x = .1
        offset_y = 0
    elif prot == 'Q07954':
        offset_x = .1
        offset_y = 0
    return offset_x, offset_y, arrow_t, rad

def custom_annotation_ctr_vs_sample(gene_name):
    offset_x, offset_y, arrow_t, rad = None, None, None, None
    if gene_name in ['Q02790']:
        offset_y = 0.2
        offset_x = 0.2
        rad = 0
    if gene_name in ['Q99613']:
        offset_y = 0.1
        offset_x = 0.15
        rad = -.2
    if gene_name in ['Q07866']:
        offset_y = 0.2
        offset_x = 0.1
        rad = -.2
    if gene_name in ['P00568']:
        offset_y = 0.15
        offset_x = 0.15
        rad = -.2
    if gene_name in ['P02652']:
        offset_x = 0.15
    if gene_name in ['Q02218']:
        offset_y = 0.2
        offset_x = -0.2
        rad = -.2

    return offset_x, offset_y, arrow_t, rad

def plot_roles(data_ctr, data_treatment, top_role_change, critical_role_change, genenames):
    #- determine wether there is gonna be 3 or 4 collumns depending on if critical_role_change is empty
    if len(critical_role_change)==0:
        ncols, nrows = 3, 1
    else:
        ncols, nrows = 4, 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(3 * ncols, 3 * nrows))
    serif_font()
    #- plot the first 2 columns, role changes for ctr and sample
    for ii, data in enumerate([data_ctr, data_treatment]):
        RolePlot.plot_ctr_vs_sample(axes[ii], data, genenames, custom_annotation=custom_annotation_ctr_vs_sample)
        if ii == 1:
            axes[ii].set_ylabel('')
    #- define titles. by default, we have 3 columns. it can be 4 if critical_role_change is not empty
    titles = ['(A1) Control', '(A2) Treatment-mg', '(A3) Top role change']
    # role change from ctr to sample for top role changes
    def plot_role_change(df_role_change, axes_i):
        ax_role_change = axes[axes_i]
        RolePlot.plot_role_change(df_role_change, ax=ax_role_change, custom_annotation= custom_annotation_role_change)
        ax_role_change.set_ymargin(.4)
        ax_role_change.set_xmargin(.4)
        ax_role_change.set_ylabel('')

    plot_role_change(top_role_change, axes_i=2)
    # role change from ctr to sample for critical role changes
    if len(critical_role_change) > 0:
        plot_role_change(critical_role_change, axes_i=3)
        titles.append('(A4) Critical role change')
    #- set the titles
    for ii, ax in enumerate(axes):
        ax.set_title(titles[ii], fontweight='bold')
    # legends. ctr and mg as well as role names
    handles = []
    for i, color in enumerate(RolePlot.ctr_sample_colors):
        handles.append(
            axes[-1].scatter([], [], marker='o', label=studies[i], s=100, edgecolor='black', color=color,
                                   linewidth=.2))
    handles.append(axes[-1].scatter([], [], marker='', label='--', s=1,
                                          edgecolor='white', color='white',
                                          linewidth=.5))  # create a blank place

    handles.extend(RolePlot.create_role_legends(axes[-1]))
    axes[-1].legend(handles=handles, loc='upper center',
                          bbox_to_anchor=(1.3, 1), prop={'size': 10}, frameon=False
                          )
    return fig
def change_protnames_to_genenames_in_links(links, map_protnames_genenames):
    """Changes regulator and target names from protnames to genenames"""
    regulators = links['Regulator']
    targets = links['Target']
    regulators = [map_protnames_genenames[protname] for protname in regulators]
    targets = [map_protnames_genenames[protname] for protname in targets]
    links['Regulator'] = regulators
    links['Target'] = targets

    return links
def retreive_links_with_genenames(selected_models, map_protnames_genenames, studies):
    """Reads links for selected models and both studies, convert the links annotation from protname
    to genename and returns links"""
    links_all = {}
    for model_name in selected_models:
        #- break the model's name into DE_type format
        model_name_parts = model_name.split('_')
        period, imput, method = model_name_parts

        DE_type = '_'.join([period, imput])
        #- VSA analysis for both studies
        links_studies = []
        for study in studies:
            links = read_from_file(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv')
            #- change the names from protnames to genenames
            links = change_protnames_to_genenames_in_links(links, map_protnames_genenames)
            links_studies.append(links)
        links_all[model_name] = links_studies
    return links_all

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--top_quantile_role_change', type=float, default=.9, help="Top quantile of biggest distance in the role change from ctr to sample")
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    top_quantile_role_change = args.top_quantile_role_change
    #- create VSA dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)
    # - load names of selected models
    selected_models = np.loadtxt(os.path.join(MODELSELECTION_DIR, f'selected_models.txt'), dtype=str, delimiter=",")

    #- read the links and store it for selected models
    links_all: Dict[str, List[pd.DataFrame]] = retreive_links_with_genenames(selected_models, F_protnames_to_genenames(), studies)

    #- the main function for VSA analysis.
    for model_name in selected_models:
        DE_type = '_'.join(model_name.split('_')[0:2])
        genenames = F_DE_genenames()[DE_type]
        #- VSA analysis for both studies
        vsa_results_studies = []
        links_studies = links_all[model_name]
        for study_i, _ in enumerate(studies):
            vsa_results = role_analysis(links_studies[study_i], genenames)
            vsa_results_studies.append(vsa_results)
        #- determine top role change from ctr to sample
        top_role_change = determine_top_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                       top_quantile=top_quantile_role_change)
        #- determine those proteins with a critical role change from ctr to sample
        critical_role_change = determine_critical_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                                     target_role=3)
        #- write VSA analysis for ctr and sample to files
        for study_i, study in enumerate(studies):
            write_to_file(vsa_results_studies[study_i], Path(VSA_DIR) / f'vsa_{model_name}_{study}.csv')
        #- write top role changes and critical role changes to file
        write_to_file(top_role_change, Path(VSA_DIR)/ f'top_role_change_{model_name}.csv')
        write_to_file(critical_role_change, Path(VSA_DIR) / f'critical_role_change_{model_name}.csv')

        #- plot roles for ctr and sample, top role changes and critical role changes in one plot
        fig= plot_roles(vsa_results_studies[0], vsa_results_studies[1], top_role_change, critical_role_change, genenames)
        fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.pdf'))

