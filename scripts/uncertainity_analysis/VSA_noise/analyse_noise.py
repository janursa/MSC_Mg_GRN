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
from typing import List, Dict, Tuple, TypeAlias
from common_tools.VSA import role_analysis

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from imports import F_model_name_2_method_and_DE_type, F_DE_genenames, F_selected_models, VSA_NOISE_DIR
from common_tools import read_write_data
from VSA.vsa_aux import retrieve_VSA_results
from common_tools.uncertainity import create_noisy_links
from uncertainity_analysis.VSA_noise.uavsa_aux import re_organize_vsa_results, determine_sig_change_in_role, plot_results
from common_tools.links import run_portia

def parse_args():
    class Args:
        noise_types = ['mp', 'ad']
        noise_intensities = {'mp':[0.05], 'ad':[.1]}
        robust_thresholds = {'mp':0.05, 'ad':.1}
        n_repeat = 500
        studies = ['ctr', 'mg']
    return Args()
def create_dir(model_name):
    results_save_dir = Path(VSA_NOISE_DIR) / model_name
    if not os.path.isdir(results_save_dir):
        os.makedirs(results_save_dir)
    plots_save_dir = Path(results_save_dir) / 'plots'
    if not os.path.isdir(plots_save_dir):
        os.makedirs(plots_save_dir)
    return results_save_dir, plots_save_dir
def run_noise_analysis(model_name: str, noise_type:str,
                       noise_intensity:float, save_file:str, genes_under_study: List[str]) -> None:
    """Main function to run analysis for a given model, noise type, and noise intensity"""
    print(f'Noise analysis for {model_name} {noise_type} {noise_intensity}')
    method, DE_type = F_model_name_2_method_and_DE_type(model_name)
    genenames = F_DE_genenames()[DE_type]
    # - get data
    data_ctr, data_sample = [read_write_data(mode='read', tag=f'{DE_type}_{study}') for study in args.studies]
    # - create noisy links for ctr and sample
    kwargs = dict(noise_type=noise_type, n_repeat=args.n_repeat,
                  noise_intensity=noise_intensity, gene_names=genenames, grn_function=run_portia)
    noisy_links_stack_ctr = create_noisy_links(data_org=data_ctr, **kwargs)
    noisy_links_stack_sample = create_noisy_links(data_org=data_sample, **kwargs)
    # - VSA
    vsa_results_stack_ctr = [role_analysis(links, genenames) for links in noisy_links_stack_ctr]
    vsa_results_stack_sample = [role_analysis(links, genenames) for links in noisy_links_stack_sample]
    # - re structred the data to have: gene: {ctr:[(Q1, P1), ...], sample:[(Q1, P1), ..]}
    gene_based_vsa_results = re_organize_vsa_results(vsa_results_stack_ctr, vsa_results_stack_sample,
                                                     genes_to_extract=genes_under_study)
    # - save to file
    with open(save_file, 'w') as ff:
        json.dump(gene_based_vsa_results, ff)

def determine_robust_genes(_sig_signs: List[str], _genes_under_study: List[str]) -> List[str]:
    """To determine which genes are robust
    For any given sig signs, *, **, **, a gene is robust
    """
    robust_indices = [gene_i for gene_i, sig_sig in enumerate(_sig_signs) if (sig_sig != '')]
    return np.asarray(_genes_under_study)[robust_indices]
if __name__ == '__main__':
    # - parse the arguments
    args = parse_args()
    force = False # force run
    # - load names of selected models
    selected_models = F_selected_models()
    # - run for each model
    targets_genes = {}
    for model_name in selected_models:
        robust_genes_model = []
        save_dir, save_plots_dir, = create_dir(model_name)
        genes_under_study = retrieve_VSA_results(model_name)
        # - run for each noise type
        for noise_type in args.noise_types:
            sig_signs_noisetype = {intensity:[] for intensity in args.noise_intensities[noise_type]}  # for final evaluation of robust genes
            # - run for each noise intensity
            for noise_intensity in args.noise_intensities[noise_type]:
                save_file = f'{save_dir}/results_{noise_type}_{noise_intensity}.json' # file to save the results
                # - if it doesn't exist, run the analysis
                if not os.path.exists(save_file) or force:
                    run_noise_analysis(model_name, noise_type, noise_intensity, save_file, genes_under_study)
                # - retrieve the results
                with open(save_file, 'r') as ff:
                        gene_based_vsa_results = json.load(ff)
                # - determine sig level from ctr to sample: *, **, ***
                sig_signs = determine_sig_change_in_role(gene_based_vsa_results)
                # - plot scatters showing ctr to sample change
                fig = plot_results(gene_based_vsa_results, noise_type, noise_intensity, sig_signs)
                fig.savefig(Path(save_plots_dir) / f'{noise_type}_{noise_intensity}.pdf')
                fig.savefig(Path(save_plots_dir) / f'{noise_type}_{noise_intensity}.png', dpi=150, transparent=True)
                # - save sig signs for rubustness analysis
                sig_signs_noisetype[noise_intensity] = sig_signs
            #- determine robust genes for noise type. those sig signs with any of *, **, ***
            threshold = args.robust_thresholds[noise_type]  # the noise intensity that we consider for robustness
            robust_genes = determine_robust_genes(sig_signs_noisetype[threshold], genes_under_study)
            robust_genes_model.append(robust_genes)
        targets_genes[model_name] = list(set(np.asarray(robust_genes_model).flatten()))
    with open(f"{VSA_NOISE_DIR}/target_genes.json", 'w') as file:
        json.dump(targets_genes, file)
