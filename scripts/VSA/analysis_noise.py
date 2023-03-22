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

# from scripts.imports import VSA_DIR, F_DE_protiens,CALIBRATION_DIR, time_points
# from scripts.utils.VSA import NoiseAnalysis
# from scripts.utils.links import run_generni, run_portia
# from scripts.utils import read_write_data
# from scripts.utils.calibration import retrieve_data

from imports import VSA_DIR, F_DE_protiens,CALIBRATION_DIR, time_points
from utils.VSA import NoiseAnalysis
from utils.links import run_generni, run_portia
from utils import read_write_data
from utils.calibration import retrieve_data


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
def output_file_name(noise_type, model_name, std):
    return f'{VSA_NOISE_DIR}/results_{noise_type}_{model_name}_{std}.json'

def main_function(data_studies, noise_type, n_rep, std_noise, warm_start, run_plot, run_analysis, model_name, gene_names, target_genes):
    print(f'Analysis for noise {noise_type}')
    if warm_start:
        file_exist = True
        if not os.path.exists(output_file_name(noise_type, model_name, str(std_noise))):
            print(f'No preivous file for warm start')
            file_exist = False
            remaining_n_rep = n_rep
        else:
            with open(output_file_name(noise_type, model_name, str(std_noise)), "r") as outfile:
                previous_results = json.load(outfile)
            previous_n_rep = len(previous_results[0]['ctr'])
            remaining_n_rep = n_rep - previous_n_rep
            print(f'warm start with {previous_n_rep} from previous run')
    else:
        remaining_n_rep = n_rep

    noise_analysis_obj = NoiseAnalysis(data_ctr=data_studies[0], data_sample=data_studies[1],
                                       target_genes=target_genes,
                                       run_grn_func=run_grn_func,
                                       n_rep=remaining_n_rep, std_noise=std_noise, noise_type=noise_type,
                                       kwargs_grn_func=dict(method=method, DE_type=DE_type,
                                                            gene_names=gene_names),
                                       kwargs_role_analysis={'gene_names': gene_names})
    if run_analysis & (remaining_n_rep > 0):
        print('Running analysis')
        results = noise_analysis_obj.analyse_noise()

        if warm_start and file_exist:
            results = merge_results(results, previous_results)
        # - save files
        with open(output_file_name(noise_type, model_name, str(std_noise)), "w") as outfile:
            json.dump(results, outfile)
    if run_plot:
        print('Visualization')
        # - plot
        with open(output_file_name(noise_type, model_name, str(std_noise)), "r") as outfile:
            results = json.load(outfile)

        fig = noise_analysis_obj.plot_results(results)
        # fig2 = noise_analysis_obj.plot_results_significant(results)

        TO_SAVE_pdfs = os.path.join(VSA_NOISE_DIR, model_name, 'pdfs')
        TO_SAVE_pngs = os.path.join(VSA_NOISE_DIR, model_name, 'pngs')
        if not os.path.isdir(TO_SAVE_pdfs):
            os.makedirs(TO_SAVE_pdfs)
        if not os.path.isdir(TO_SAVE_pngs):
            os.makedirs(TO_SAVE_pngs)

        fig.savefig(Path(TO_SAVE_pdfs)/f'{noise_type}_{std_noise}.pdf')
        fig.savefig(Path(TO_SAVE_pngs)/f'{noise_type}_{std_noise}.png', dpi=290, transparent=True)

        # fig2.savefig(Path(TO_SAVE_pdfs) / f'{noise_type}_{std_noise}_sig.pdf')
        # fig2.savefig(Path(TO_SAVE_pngs) / f'{noise_type}_{std_noise}_sig.png', dpi=300, transparent=True)


if __name__ == '__main__':
    #- create dir
    VSA_NOISE_DIR = Path(VSA_DIR)/'noise'
    if not os.path.isdir(VSA_NOISE_DIR):
        os.makedirs(VSA_NOISE_DIR)


    studies = ['ctr', 'mg']
    methods = ['RF','ridge','portia']

    selected_models = ['day1_11_KNN_RF', 'day1_21_KNN_portia']
    # selected_models = ['day1_11_KNN_RF']

    grn_functions = {'portia': run_portia,
                     'RF': run_RF
                     }

    n_rep = 100
    std_mnoises = [0, 0.02, 0.05, .1]  # multiplicative noise std
    std_anoises = [0, 0.02, 0.05, 0.1, 0.2]
    # std_mnoises = [0.05]  # multiplicative noise std
    # std_anoises = [0.1]

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
            # target_genes = list(set(target_genes))
            data_studies = []
            for study in studies:
                data_org = read_write_data(mode='read', tag=f'{DE_type}_{study}')
                data_studies.append(data_org)

            # ------------------- noise analysis--------------------------
            for std_mnoise in std_mnoises:
                #- multiplicative
                noise_type = 'mp'
                std_noise = std_mnoise
                gene_names = DE_proteins
                main_function(data_studies, noise_type, n_rep, std_noise, warm_start, run_plot,
                              run_analysis, model_name, gene_names, target_genes)
            for std_anoise in std_anoises:
                #- additive
                noise_type = 'ad'
                std_noise = std_anoise
                gene_names = DE_proteins
                main_function(data_studies, noise_type, n_rep, std_noise, warm_start, run_plot,
                              run_analysis, model_name, gene_names, target_genes)