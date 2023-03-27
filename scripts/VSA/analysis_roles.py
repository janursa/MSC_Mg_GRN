"""
    Conduct vester's sensitivity analysis for non noised and noised links
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import GRN_DIR, VSA_DIR, F_DE_proteins
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


def plot_roles(data_ctr, data_treatment, top_role_change, critical_role_change):

    # ctr vs treatment
    ncols, nrows = 4, 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(3 * ncols, 3 * nrows))
    serif_font()
    for ii, data in enumerate([data_ctr, data_treatment]):
        RolePlot.plot_ctr_vs_sample(axes[ii], data, DE_proteins, custom_annotation=custom_annotation_ctr_vs_sample)
        if ii == 1:
            axes[ii].set_ylabel('')
    # role change from ctr to treatment
    for ii, df_role_change in enumerate([ critical_role_change, top_role_change]):
        ax_role_change = axes[2 + ii]
        RolePlot.plot_role_change(df_role_change, ax=ax_role_change, custom_annotation= custom_annotation_role_change)
        ax_role_change.set_ymargin(.4)
        ax_role_change.set_xmargin(.4)
        ax_role_change.set_ylabel('')
    titles = ['(A1) Control', '(A2) Treatment-mg', '(A3) Critical role change', '(A4) Top role change', ]
    for ii, ax in enumerate(axes):
        ax.set_title(titles[ii], fontweight='bold')
    # legends for ctr and mg as well as vsa roles
    handles = []
    for i, color in enumerate(RolePlot.ctr_sample_colors):
        handles.append(
            ax_role_change.scatter([], [], marker='o', label=studies[i], s=100, edgecolor='black', color=color,
                                   linewidth=.2))
    handles.append(ax_role_change.scatter([], [], marker='', label='--', s=1,
                                          edgecolor='white', color='white',
                                          linewidth=.5))  # create a blank place

    handles.extend(RolePlot.create_role_legends(ax_role_change))
    ax_role_change.legend(handles=handles, loc='upper center',
                          bbox_to_anchor=(1.3, 1), prop={'size': 10}, frameon=False
                          )
    return fig
if __name__ == '__main__':
    #- create dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)

    studies = ['ctr', 'mg']
    methods = ['RF','ridge','portia']

    # selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    selected_models = ['day1_21_KNN_portia']

    model_i = 0
    for idx, (DE_type, DE_proteins) in enumerate(F_DE_proteins().items()):
        for method in methods:
            model_name = '_'.join([DE_type,method])
            if model_name not in selected_models: #only selected models
                continue
            # -------- analyse --------------
            vsa_results_studies = []
            for study in studies:
                links = read_from_file(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv')
                vsa_results = role_analysis(links, DE_proteins)
                vsa_results_studies.append(vsa_results)

            top_role_change = determine_top_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                           top_quantile=.9)
            critical_role_change = determine_critical_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                                         target_role=3)
            #- write to file
            write_to_file(vsa_results_studies[0], Path(VSA_DIR) / f'vsa_{model_name}_ctr.csv')
            write_to_file(vsa_results_studies[1], Path(VSA_DIR) / f'vsa_{model_name}_mg.csv')
            write_to_file(top_role_change, Path(VSA_DIR)/f'top_role_change_{model_name}.csv')
            write_to_file(critical_role_change, Path(VSA_DIR) / f'critical_role_change_{model_name}.csv')

            # -------- plot --------------
            fig= plot_roles(vsa_results_studies[0], vsa_results_studies[1], top_role_change, critical_role_change)
            fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.png'), dpi=300, transparent=True)
            fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.pdf'))

