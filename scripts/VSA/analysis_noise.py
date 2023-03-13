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


def merge_results(new, previous):
    """ merges new and old results for warm analysis
    results are organized as each item in the list is for a gene. each item has ctr and sample, where each
    has n elements of (Q,P)
    """
    merged = []
    for i_gene, _ in enumerate(target_genes):
        merged_i = {study: new[i_gene][study] + previous[i_gene][study] for study in new[0].keys()}
        merged.append(merged_i)
    return merged

def write_to_file(df, file_name):
    with open(file_name, 'w') as f:
        df.to_csv(f, index=False)
def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data
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


def run_grn_func(**kwargs):
    return grn_functions[method](**kwargs)
def output_file_name(model_name, std):
    return f'{VSA_NOISE_DIR}/results_mp_{model_name}_{std}.json'

if __name__ == '__main__':
    #- create dir
    VSA_NOISE_DIR = Path(VSA_DIR)/'noise'
    if not os.path.isdir(VSA_NOISE_DIR):
        os.makedirs(VSA_NOISE_DIR)
    studies = ['ctr', 'mg']
    methods = ['RF','ridge','portia']

    # selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    selected_models = ['day1_11_KNN_RF']

    grn_functions = {'portia': run_portia,
                     'RF': run_RF
                     }

    n_rep = 5
    std_mnoise = 0.02  # multiplicative noise std
    std_anoise = 0.1

    run_analysis = True
    run_plot = True
    warm_start = True # to use the previous run's results

    model_i = 0
    for idx, (DE_type, DE_proteins) in enumerate(F_DE_protiens().items()):
        for method in methods:
            model_name = '_'.join([DE_type,method])
            if model_name not in selected_models: #only selected models
                continue

            top_role_change = read_from_file(Path(VSA_DIR) / f'top_role_change_{model_name}.csv')
            critical_role_change = read_from_file(Path(VSA_DIR) / f'critical_role_change_{model_name}.csv')
            target_genes = np.concatenate((top_role_change['Entry'].to_numpy(), critical_role_change['Entry'].to_numpy()))
            data_studies = []
            for study in studies:
                data_org = read_write_data(mode='read', tag=f'{DE_type}_{study}')
                data_studies.append(data_org)
            # ------------------- noise analysis--------------------------
            if warm_start:
                file_exist = True
                if not os.path.isdir(output_file_name(model_name, str(std_mnoise))):
                    print('No preivous file for warm start')
                    file_exist = False
                    remaining_n_rep = n_rep
                if file_exist:
                    with open(output_file_name(model_name, str(std_mnoise)), "r") as outfile:
                        previous_results_mp = json.load(outfile)
                    with open(output_file_name(model_name, str(std_anoise)), "r") as outfile:
                        previous_results_ad = json.load(outfile)
                    previous_n_rep = len(previous_results_mp[0]['ctr'])
                    remaining_n_rep = n_rep - previous_n_rep
                    print(f'warm start with {previous_n_rep} from previous run')
            else:
                remaining_n_rep = n_rep

            noise_analysis_obj = NoiseAnalysis(data_ctr=data_studies[0], data_sample=data_studies[1],
                                               target_genes=target_genes,
                                               run_grn_func=run_grn_func,
                                               n_rep=remaining_n_rep, std_mpnoise=std_mnoise, std_adnoise=std_anoise,
                                               kwargs_grn_func=dict(method=method, DE_type=DE_type,
                                                                    gene_names=DE_proteins),
                                               kwargs_role_analysis={'gene_names': DE_proteins})
            if run_analysis & (remaining_n_rep>0):
                results_mp = noise_analysis_obj.analyse_mp_noise()
                results_ad = noise_analysis_obj.analyse_ad_noise()

                if warm_start and file_exist:
                    results_mp = merge_results(results_mp, previous_results_mp)
                    results_ad = merge_results(results_ad, previous_results_ad)
                #- save files
                with open(output_file_name(model_name, str(std_mnoise)), "w") as outfile:
                    json.dump(results_mp, outfile)
                with open(output_file_name(model_name, str(std_anoise)), "w") as outfile:
                    json.dump(results_ad, outfile)
            if run_plot:
                #- plot
                with open(output_file_name(model_name, str(std_mnoise)), "r") as outfile:
                    results_mp = json.load(outfile)
                with open(output_file_name(model_name, str(std_anoise)), "r") as outfile:
                    results_ad = json.load(outfile)
                fig = noise_analysis_obj.plot_results(results_mp, results_ad)

                fig.savefig(os.path.join(VSA_NOISE_DIR, f'noiseAnalysis_{model_name}.pdf'))
                fig.savefig(os.path.join(VSA_NOISE_DIR, f'noiseAnalysis_{model_name}.png'), dpi=300, transparent=True)

