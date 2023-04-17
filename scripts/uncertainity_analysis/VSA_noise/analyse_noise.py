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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from imports import F_model_name_2_method_and_DE_type, F_DE_genenames, F_selected_models, VSA_NOISE_DIR
from common_tools import read_write_data
from VSA.vsa_aux import retrieve_VSA_results, run_noise_analysis, postprocess_noiseanalysis_results



if __name__ == '__main__':
    #- parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--n_noisy_datasets', type=int, default=100, help="Number of noisy data to produce")
    parser.add_argument('--stds_mnoise', nargs='+', type=float, default=[0.0, 0.02, 0.05, .1], help="Standard deviation of multiplicative gaussian noise")
    parser.add_argument('--stds_anoise', nargs='+', type=float, default=[0.0, 0.02, 0.05, 0.1, 0.2], help="Standard deviation of additive gaussian noise")
    parser.add_argument('--astd_sig_check', type=float, default=.1, help="Standard deviation of multiplicative gaussian noise to evaluate significant change in protein roles")
    parser.add_argument('--mstd_sig_check', type=float, default=.05, help="Standard deviation of additive gaussian noise to evaluate significant change in protein roles")

    parser.add_argument('--warm_start', type=bool, default=True, help="To use previous results.")

    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    n_noisy_datasets = args.n_noisy_datasets
    stds_mnoise = args.stds_mnoise
    stds_anoise = args.stds_anoise
    warm_start = args.warm_start
    mstd_sig_check = args.mstd_sig_check
    astd_sig_check = args.astd_sig_check

    # - load names of selected models
    selected_models = F_selected_models()

    #- create noise dir
    if not os.path.isdir(VSA_NOISE_DIR):
        os.makedirs(VSA_NOISE_DIR)
    VSA_NOISE_RUNS_DIR = Path(VSA_NOISE_DIR)/'runs'
    if not os.path.isdir(VSA_NOISE_RUNS_DIR):
        os.makedirs(VSA_NOISE_RUNS_DIR)
    VSA_NOISE_PLOTS_DIR = Path(VSA_NOISE_DIR)/'plots'
    if not os.path.isdir(VSA_NOISE_PLOTS_DIR):
        os.makedirs(VSA_NOISE_PLOTS_DIR)
    #- define functions that run GRN inference. We need them to introdue noise and run the models

    sig_target_genes_all = {}
    for model_name in selected_models:
        print(f'Noise analysis for {model_name}')
        #- determine method and DE_type from model name
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        genenames = F_DE_genenames()[DE_type]
        # - retreive model VSA results
        target_genes = retrieve_VSA_results(model_name)
        #- get proteomics data
        data_studies = []
        for study in studies:
            data_org = read_write_data(mode='read', tag=f'{DE_type}_{study}')
            data_studies.append(data_org)
        #- run the function for multiplicative noise
        run_noise_analysis(noise_type='mp', stds_noise=stds_mnoise,
                           data_ctr=data_studies[0],
                           data_sample=data_studies[1],
                           n_noisy_datasets=n_noisy_datasets,
                           warm_start=warm_start,
                           model_name=model_name, gene_names=genenames,
                           target_genes=target_genes, method=method,
                           DE_type=DE_type,
                           noise_runs_dir=VSA_NOISE_RUNS_DIR)
        #- run the function for additive noise
        run_noise_analysis(noise_type='ad', stds_noise=stds_anoise,
                           data_ctr=data_studies[0],
                           data_sample=data_studies[1],
                           n_noisy_datasets=n_noisy_datasets,
                           warm_start=warm_start,
                           model_name=model_name, gene_names=genenames,
                           target_genes=target_genes, method=method,
                           DE_type=DE_type,
                           noise_runs_dir=VSA_NOISE_RUNS_DIR)
        # - plot results and determine sig signs
        sig_signs_allstds_mp = postprocess_noiseanalysis_results(noise_type='mp', stds_noise=stds_mnoise, target_genes=target_genes, model_name=model_name, noise_plots_dir=VSA_NOISE_PLOTS_DIR, noise_runs_dir=VSA_NOISE_RUNS_DIR)
        sig_signs_allstds_ad = postprocess_noiseanalysis_results(noise_type='ad', stds_noise=stds_anoise, target_genes=target_genes, model_name=model_name, noise_plots_dir=VSA_NOISE_PLOTS_DIR, noise_runs_dir=VSA_NOISE_RUNS_DIR)
        #- determine sig changes. those sig signs with any of *, **, ***
        def find_sig_targetgenes(sig_signs_allstds, std_sig_check):
            sig_signs = sig_signs_allstds[std_sig_check]
            selected_indices = []
            for index, sig_sig in enumerate(sig_signs):
                if sig_sig != '':
                    selected_indices.append(index)
            return np.asarray(target_genes)[selected_indices]

        sig_genenames_mp = find_sig_targetgenes(sig_signs_allstds_mp, mstd_sig_check)
        sig_genenames_ad = find_sig_targetgenes(sig_signs_allstds_ad, astd_sig_check)

        sig_target_genes = list(set(sig_genenames_mp) & set(sig_genenames_ad))

        sig_target_genes_all[model_name] = sig_target_genes
    with open(Path(VSA_NOISE_DIR)/ "target_genes.json", 'w') as file:
        json.dump(sig_target_genes_all, file)
