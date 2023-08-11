import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path


from proteomics_MSC.imports import DIVERGENCE_ANALYSIS_DIR, DIVERGENCE_ANALYSIS_ROBUSTNESS_DIR
from proteomics_MSC.common_tools import F_model_name_2_method_and_DE_type, F_selected_models, change_protnames_to_genenames_in_links
from proteomics_MSC.common_tools.noise_analysis import NetworkNoiseAnalysis
from proteomics_MSC.VSA.vsa_aux import retreive_links_with_genenames
from proteomics_MSC.common_tools.links import choose_top_quantile
# def shortlist(model_name, studies, divergence_scores_df, top_percentile):
#     """Select the top percentile from the links"""
#     links_studies = retreive_links_with_genenames(model_name, studies)
#     # - assert
#     assert (len(links_studies[0]) == len(divergence_scores_df))
#     # - calculate abs difference between ctr and sample
#     divergence_scores_df['difference']  = np.abs(links_studies[0]['Weight']-links_studies[1]['Weight'])
#     # - choose the top percentile based on difference
#     divergence_scores_short.drop(columns=['difference'], inplace=True)
#     return divergence_scores_short


import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':

    force = True #To overwrite if the result files exist'
    n_repeat = config['n_repeat_noise_analysis']

    studies = config['studies']
    # - run noise analysis to determine distance scores from ctr to sample
    for model_name in F_selected_models():
        network_noise_analysis_obj = NetworkNoiseAnalysis(model_name=model_name, studies=studies,
                                                          base_dir=DIVERGENCE_ANALYSIS_ROBUSTNESS_DIR, force=force, n_repeat=n_repeat)
        # - run the analysis and post analysis
        network_noise_analysis_obj.run_decorator()
        # - calculate overall distance scores for this model by accounting for different noise types
        divergence_scores_df = network_noise_analysis_obj.determine_overall_divergence_scores()
        # - shortlist to those identified as top role change previously
        top_links = pd.read_csv(Path(DIVERGENCE_ANALYSIS_DIR) / f'top_links_{model_name}.csv')
        divergence_scores_short = pd.merge(divergence_scores_df, top_links, on=['Regulator', 'Target'])
        # - save
        to_save = Path(DIVERGENCE_ANALYSIS_ROBUSTNESS_DIR) / model_name / 'links_divergence_scores.csv'
        divergence_scores_short.to_csv(to_save, index=False)
        print(f'\n output -> {to_save}')