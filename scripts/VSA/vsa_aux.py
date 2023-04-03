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

from imports import F_model_name_2_method_and_DE_type, VSA_DIR, GRN_DIR, CALIBRATION_DIR, time_points, F_selected_models, VSA_NOISE_DIR
from utils.VSA import NoiseAnalysis, RolePlot
from utils.links import run_generni, run_portia
from utils import read_write_data,serif_font
from utils.calibration import retrieve_data


def merge_results(new, previous, target_genes):
    """ merges new and old results for warm analysis
    results are organized as each item in the list is for a gene. each item has ctr and sample, where each
    has n elements of (Q,P)
    """
    merged = []
    for i_gene, _ in enumerate(target_genes):
        merged_i = {study: new[i_gene][study] + previous[i_gene][study] for study in new[0].keys()}
        merged.append(merged_i)
    return merged
def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data



def read_from_file(file_name):
    with open(file_name, 'r') as f:
        data = pd.read_csv(f, index_col=False)
    return data

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
        method, DE_type = F_model_name_2_method_and_DE_type(model_name)
        period, imput = DE_type.split('_')

        #- VSA analysis for both studies
        links_studies = []
        for study in studies:
            links = read_from_file(Path(GRN_DIR) / method / f'links_{DE_type}_{study}.csv')
            #- change the names from protnames to genenames
            links = change_protnames_to_genenames_in_links(links, map_protnames_genenames)
            links_studies.append(links)
        links_all[model_name] = links_studies
    return links_all
def run_grn_func(method, **kwargs):
    grn_functions = {'portia': run_portia}
    return grn_functions[method](**kwargs)
def output_file_name(noise_runs_dir,noise_type, model_name, std):
    return f'{noise_runs_dir}/results_{noise_type}_{model_name}_{std}.json'

def run_noise_analysis(data_ctr, data_sample, noise_type, n_noisy_datasets, stds_noise, warm_start, model_name, gene_names, target_genes, method, DE_type, noise_runs_dir):
    # for each std
    for std_noise in stds_noise:
        #- if warm start, get the status of previous runs
        if warm_start:
            file_exist = True
            if not os.path.exists(output_file_name(noise_runs_dir, noise_type, model_name, str(std_noise))):
                print(f'No preivous file for warm start')
                file_exist = False
                remaining_n_noisy_datasets = n_noisy_datasets
            else:
                with open(output_file_name(noise_runs_dir, noise_type, model_name, str(std_noise)), "r") as outfile:
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
                results = merge_results(results, previous_results, target_genes)
            # - save files
            with open(output_file_name(noise_runs_dir, noise_type, model_name, str(std_noise)), "w") as outfile:
                json.dump(results, outfile)
def determine_( sig_std=.1):
    """sig_std: at this std, we determine if role change is sig"""
def postprocess_noiseanalysis_results(noise_type, stds_noise, target_genes, model_name, noise_plots_dir, noise_runs_dir):
    """ Postprocess the results of noise analysis.
    Determines sig level in role change from ctr to mg by assigning *, **, ***.
    Scatter plots the roles for ctr and sampls on the sample graph with different colors for
    all noisy data
    """
    #- create plot dir
    # - for each model_name, we store pdfs an pngs in the same-name folder
    TO_SAVE_pdfs = os.path.join(noise_plots_dir, model_name, 'pdfs')
    TO_SAVE_pngs = os.path.join(noise_plots_dir, model_name, 'pngs')
    if not os.path.isdir(TO_SAVE_pdfs):
        os.makedirs(TO_SAVE_pdfs)
    if not os.path.isdir(TO_SAVE_pngs):
        os.makedirs(TO_SAVE_pngs)
    #- post process for each std
    sig_signs_allstds = {}
    for std_noise in stds_noise:
        # - retreive files from noise analysis
        with open(output_file_name(noise_runs_dir, noise_type, model_name, str(std_noise)), "r") as outfile:
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
    # target_genes = list(set(target_genes))
    return target_genes
