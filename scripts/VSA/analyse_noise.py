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
from typing import List, Dict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import VSA_DIR, F_DE_genenames, CALIBRATION_DIR, time_points, MODELSELECTION_DIR, VSA_NOISE_DIR, F_protnames_to_genenames
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
def run_RF(data, i_data, DE_type, gene_names, **kwargs):
    """ Runs random forest GRN inference
    Receives proteomics data and rund GRN inference using geneRNI
    """
    param = dict(estimator_t = 'RF')
    method = 'RF'
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
    return f'{VSA_NOISE_RUNS_DIR}/results_{noise_type}_{model_name}_{std}.json'

def run_noise_analysis(data_ctr, data_sample, noise_type, n_noisy_datasets, stds_noise, warm_start, model_name, gene_names, target_genes, method, DE_type):
    # for each std
    for std_noise in stds_noise:
        #- if warm start, get the status of previous runs
        if warm_start:
            file_exist = True
            if not os.path.exists(output_file_name(noise_type, model_name, str(std_noise))):
                print(f'No preivous file for warm start')
                file_exist = False
                remaining_n_noisy_datasets = n_noisy_datasets
            else:
                with open(output_file_name(noise_type, model_name, str(std_noise)), "r") as outfile:
                    previous_results = json.load(outfile)
                previous_n_noisy_datasets = len(previous_results[0]['ctr'])
                remaining_n_noisy_datasets = n_noisy_datasets - previous_n_noisy_datasets
                print(f'warm start with {previous_n_noisy_datasets} from previous run')
        else:
            remaining_n_noisy_datasets = n_noisy_datasets
        #- contruct noise object
        noise_analysis_obj = NoiseAnalysis(data_ctr=data_ctr, data_sample=data_sample,
                                           target_genes=target_genes,
                                           run_grn_func=run_grn_func,
                                           n_noisy_datasets=remaining_n_noisy_datasets, std_noise=std_noise, noise_type=noise_type,
                                           kwargs_grn_func=dict(method=method, DE_type=DE_type,
                                                                gene_names=gene_names),
                                           kwargs_role_analysis={'gene_names': gene_names})
        if remaining_n_noisy_datasets > 0:
            results = noise_analysis_obj.analyse_noise()

            if warm_start and file_exist:
                results = merge_results(results, previous_results)
            # - save files
            with open(output_file_name(noise_type, model_name, str(std_noise)), "w") as outfile:
                json.dump(results, outfile)
def determine_( sig_std=.1):
    """sig_std: at this std, we determine if role change is sig"""
def postprocess_noiseanalysis_results(noise_type, stds_noise, target_genes, model_name):
    """ Postprocess the results of noise analysis.
    Determines sig level in role change from ctr to mg by assigning *, **, ***.
    Scatter plots the roles for ctr and sampls on the sample graph with different colors for
    all noisy data
    """
    #- create plot dir
    # - for each model_name, we store pdfs an pngs in the same-name folder
    TO_SAVE_pdfs = os.path.join(VSA_NOISE_PLOTS_DIR, model_name, 'pdfs')
    TO_SAVE_pngs = os.path.join(VSA_NOISE_PLOTS_DIR, model_name, 'pngs')
    if not os.path.isdir(TO_SAVE_pdfs):
        os.makedirs(TO_SAVE_pdfs)
    if not os.path.isdir(TO_SAVE_pngs):
        os.makedirs(TO_SAVE_pngs)
    #- post process for each std
    sig_signs_allstds = {}
    for std_noise in stds_noise:
        # - retreive files from noise analysis
        with open(output_file_name(noise_type, model_name, str(std_noise)), "r") as outfile:
            results = json.load(outfile)
        #- determine sig level for role change from ctr to sample. *, **, ***
        sig_signs =  NoiseAnalysis.determine_sig_change_in_role(results, target_genes)
        sig_signs_allstds[std_noise] = sig_signs
        #- plot the results
        print(f'Plot {model_name} {noise_type}')
        fig = NoiseAnalysis.plot_results(results, target_genes, noise_type, std_noise, sig_signs)
        fig.savefig(Path(TO_SAVE_pdfs) / f'{noise_type}_{std_noise}.pdf')
        fig.savefig(Path(TO_SAVE_pngs) / f'{noise_type}_{std_noise}.png', dpi=290, transparent=True)
        plt.close(fig)
    return sig_signs_allstds

def retrieve_VSA_results(model_name):
    """ Get the results of VSA analysis
    """
    top_role_change = read_from_file(Path(VSA_DIR) / f'top_role_change_{model_name}.csv')
    critical_role_change = read_from_file(Path(VSA_DIR) / f'critical_role_change_{model_name}.csv')
    target_genes = np.concatenate((top_role_change['Entry'].to_numpy(), critical_role_change['Entry'].to_numpy()))
    return target_genes

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
    selected_models = np.loadtxt(os.path.join(MODELSELECTION_DIR, f'selected_models.txt'), dtype=str, delimiter=",")

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
    grn_functions = {'portia': run_portia, 'RF': run_RF}

    for model_name in selected_models:
        print(f'Noise analysis for {model_name}')
        #- determine method and DE_type from model name
        method = model_name.split('_')[-1]
        DE_type = '_'.join(model_name.split('_')[0:2])
        genenames = F_DE_genenames()[DE_type]
        #- get proteomics data
        data_studies = []
        for study in studies:
            data_org = read_write_data(mode='read', tag=f'{DE_type}_{study}')
            data_studies.append(data_org)
        #- retreive model VSA results
        target_genes = retrieve_VSA_results(model_name)
        #- construct the generic noise anlysis function
        run_noise_analysis_F = lambda noise_type, noise_stds: run_noise_analysis (noise_type=noise_type, stds_noise=noise_stds, data_ctr=data_studies[0], data_sample=data_studies[1], n_noisy_datasets=n_noisy_datasets,
                      warm_start=warm_start,
                       model_name=model_name, gene_names=genenames, target_genes=target_genes, method=method, DE_type=DE_type)
        #- run the function for multiplicative noise
        run_noise_analysis_F(noise_type='mp', noise_stds=stds_mnoise)
        #- run the function for additive noise
        run_noise_analysis_F(noise_type='ad', noise_stds=stds_anoise)
        # - plot results and determine sig signs
        sig_signs_allstds_mp = postprocess_noiseanalysis_results(noise_type='mp', stds_noise=stds_mnoise, target_genes=target_genes, model_name=model_name)
        sig_signs_allstds_ad = postprocess_noiseanalysis_results(noise_type='ad', stds_noise=stds_anoise, target_genes=target_genes, model_name=model_name)
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

        sig_genenames = list(set(sig_genenames_mp) & set(sig_genenames_ad))


        np.savetxt(os.path.join(VSA_NOISE_DIR, f'sig_genenames_{model_name}.txt'), sig_genenames, delimiter=",", fmt="%s")
