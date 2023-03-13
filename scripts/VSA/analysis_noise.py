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
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, VSA_DIR, F_DE_protiens, F_DE_data,CALIBRATION_DIR, time_points
from scripts.utils.VSA import NoiseAnalysis,role_analysis, RolePlot, determine_top_role_change, determine_critical_role_change
from scripts.utils import make_title_pretty, comic_font
from scripts.utils.links import choose_top_quantile, run_generni, run_portia
from scripts.utils import serif_font
from scripts.utils import read_write_data
from scripts.utils.calibration import retrieve_data



def write_to_file(df, file_name):
    with open(file_name, 'w') as f:
        df.to_csv(f, index=False)
def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data

def plot_roles(data_ctr, data_treatment, top_role_change, critical_role_change):

    # ctr vs treatment
    ncols, nrows = 4, 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(3 * ncols, 3 * nrows))
    serif_font()
    for ii, data in enumerate([data_ctr, data_treatment]):
        RolePlot.plot_ctr_vs_sample(axes[ii], data, DE_proteins)
        if ii == 1:
            axes[ii].set_ylabel('')
    # role change from ctr to treatment
    for ii, df_role_change in enumerate([top_role_change, critical_role_change]):
        ax_role_change = axes[2 + ii]
        RolePlot.plot_role_change(df_role_change, ax=ax_role_change)
        ax_role_change.set_ymargin(.4)
        ax_role_change.set_xmargin(.4)
        ax_role_change.set_ylabel('')
    titles = ['(1) Control', '(2) Treatment', '(3) Top role changes', '(4) Critical role change']
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
    VSA_NOISE_DIR = Path(VSA_DIR)/'noise'
    if not os.path.isdir(VSA_NOISE_DIR):
        os.makedirs(VSA_NOISE_DIR)
    studies = ['ctr', 'mg']
    methods = ['RF','ridge','portia']

    # selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    selected_models = ['day1_21_KNN_portia']

    run_role_analysis = True # if false, use the previous analysis results.
    run_noise_analysis = True

    def run_RF(data, i_data, DE_type, method, gene_names, **kwargs):
        param = dict(estimator_t = 'RF')
        study = ['ctr','mg'][i_data]
        n_timepoints = data.shape[0]
        days = time_points()[0:n_timepoints]
        _, param_unique = retrieve_data(study=study, DE_type=DE_type, method=method,
                                        output_dir=CALIBRATION_DIR)
        ests, train_scores, links_df, oob_scores, test_scores = run_generni(
            data=data, time_points=days, gene_names=gene_names, param=param, param_unique=param_unique)
        return links_df

    grn_functions = {'portia':run_portia,
                     'RF': run_RF
                     }
    def run_grn_func(**kwargs):
        return grn_functions[method](**kwargs)


    model_i = 0
    for idx, (DE_type, DE_proteins) in enumerate(F_DE_protiens().items()):
        for method in methods:
            model_name = '_'.join([DE_type,method])
            if model_name not in selected_models: #only selected models
                continue

            if run_role_analysis:
                # -------- analyse --------------
                vsa_results_studies = []
                for study in studies:
                    links = read_from_file(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv')
                    vsa_results = role_analysis(links, DE_proteins)
                    vsa_results_studies.append(vsa_results)
                    # if study == 'mg':
                    #     print(vsa_results.iloc[0:10,:])
                    # role change
                top_role_change = determine_top_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                               top_quantile=.9)
                critical_role_change = determine_critical_role_change(vsa_results_studies[0], vsa_results_studies[1],
                                                                             target_role=3)

                write_to_file(top_role_change, Path(VSA_DIR)/f'top_role_change_{model_name}.csv')
                write_to_file(critical_role_change, Path(VSA_DIR) / f'critical_role_change_{model_name}.csv')

                # -------- plot --------------
                fig= plot_roles(vsa_results_studies[0], vsa_results_studies[1], top_role_change, critical_role_change)
                fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.png'), dpi=300, transparent=True)
                fig.savefig(os.path.join(VSA_DIR, f'role_changes_{model_name}.pdf'))

            if run_noise_analysis:
                n_rep = 5
                std_mnoise = 0.000001 #multiplicative noise std
                std_anoise = 0.1

                top_role_change = read_from_file(Path(VSA_DIR) / f'top_role_change_{model_name}.csv')
                critical_role_change = read_from_file(Path(VSA_DIR) / f'critical_role_change_{model_name}.csv')

                # ------------------- noise analysis--------------------------

                target_genes = np.concatenate((top_role_change['Entry'].to_numpy(), critical_role_change['Entry'].to_numpy()))
                data_studies = []
                for study in studies:
                    data_org = read_write_data(mode='read', tag=f'{DE_type}_{study}')
                    data_studies.append(data_org)
                noise_analysis_obj = NoiseAnalysis(data_ctr=data_studies[0], data_sample=data_studies[1], target_genes=target_genes,
                                                   run_grn_func=run_grn_func, kwargs_grn_func = dict(method = method, DE_type=DE_type, gene_names=DE_proteins),
                                                   kwargs_role_analysis={'gene_names':DE_proteins})
                results_mp = noise_analysis_obj.analyse_mp_noise(n_rep=n_rep, std_noise=std_mnoise)
                # results_ad = noise_analysis_obj.analyse_ad_noise(n_rep=n_rep, std_noise=std_anoise)
                #- save files
                with open(f'{VSA_NOISE_DIR}/results_mp_{model_name}_{std_mnoise}.json', "w") as outfile:
                    json.dump(results_mp, outfile)
                # with open(f'{VSA_NOISE_DIR}/results_ad_{model_name}_{std_anoise}.json', "w") as outfile:
                #     json.dump(results_ad, outfile)
                #- plot
                fig = noise_analysis_obj.plot_results()

                fig.savefig(os.path.join(VSA_NOISE_DIR, f'noiseAnalysis_{model_name}.pdf'))
                fig.savefig(os.path.join(VSA_NOISE_DIR, f'noiseAnalysis_{model_name}.png'), dpi=300, transparent=True)

