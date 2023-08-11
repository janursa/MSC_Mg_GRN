import json
import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path


from proteomics_MSC.imports import TARGETS_DIR
from proteomics_MSC.common_tools import F_model_name_2_method_and_DE_type, F_selected_models, change_protnames_to_genenames_in_links
from proteomics_MSC.common_tools.noise_analysis import NetworkNoiseAnalysis
from proteomics_MSC.VSA.vsa_aux import retreive_links_with_genenames
from proteomics_MSC.common_tools.links import choose_top_quantile

def extract_target_genes_connections(df, target_genes, top_n=2):
    """Gets two most significantly altered connections of the target genes"""
    # Filter the DataFrame for target_genes in 'Target'
    target_df = df[df['Target'].isin(target_genes)].sort_values(by='divergence_score', ascending=False)
    top_target_df = target_df.groupby('Target').head(top_n)

    # Filter the DataFrame for target_genes in 'Regulator'
    regulator_df = df[df['Regulator'].isin(target_genes)].sort_values(by='divergence_score', ascending=False)
    top_regulator_df = regulator_df.groupby('Regulator').head(top_n)

    # Concatenate both top_target_df and top_regulator_df
    final_df = pd.concat([top_target_df, top_regulator_df], ignore_index=True)
    return final_df


import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':

    force = False #To overwrite if the result files exist'
    n_repeat = config['n_repeat_noise_analysis']


    studies = config['studies']
    top_quantile = config['top_quantile_links_change']

    # get target genes
    with open(f'{TARGETS_DIR}/target_genes.json', 'r') as file:
        targets_all = json.load(file)
    # - run noise analysis to determine distance scores from ctr to sample
    for model_name in F_selected_models():
        targets = targets_all[model_name]
        network_noise_analysis_obj = NetworkNoiseAnalysis(model_name=model_name, studies=studies,
                                                          base_dir=TARGETS_DIR, force=force, n_repeat=n_repeat)
        # - run the analysis and post analysis
        network_noise_analysis_obj.run_decorator()
        # - calculate overall divergence scores
        divergence_df = network_noise_analysis_obj.determine_overall_divergence_scores()
        # - get those that are connected to target genes
        divergence_short_df = extract_target_genes_connections(divergence_df, targets)

        # - save
        to_save = Path(TARGETS_DIR) / model_name / 'links_divergence_scores.csv'
        divergence_short_df.to_csv(to_save, index=False)
        print(f'\n output -> {to_save}')