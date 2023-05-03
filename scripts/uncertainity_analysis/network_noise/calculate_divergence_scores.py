import sys
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))

from imports import F_model_name_2_method_and_DE_type, F_selected_models, TOP_LINKS_DIR, TOP_LINKS_WITH_DISTANCE_SCORES_DIR, NETWORK_NOISE_DIR, NETWORK_ANALYSIS_DIR, \
    change_protnames_to_genenames_in_links
from uncertainity_analysis.ua_aux import NetworkNoiseAnalysis
from VSA.vsa_aux import retreive_links_with_genenames
from common_tools.links import choose_top_quantile

# def target_links_around_target_genes_with_divergence_scores(model_name, divergence_scores_df) -> None:
#     """ Extract top links around the target genes
#     """
#     method, DE_type = F_model_name_2_method_and_DE_type(model_name)
#     # - get the top links
#     top_links = pd.read_csv(Path(NETWORK_ANALYSIS_DIR) / f'top_links_around_target_proteins_{model_name}.csv', index_col=False)
#     # - shortlist
#     top_links_with_scores = divergence_scores_df.merge(top_links, on=['Regulator','Target'])
#     # - save
#     file_name = Path(NETWORK_ANALYSIS_DIR) / f'top_links_around_target_genes_with_divergence_scores_{model_name}.csv'
#     top_links_with_scores.to_csv(file_name, index=False)

# def divergence_score(model_name, top_links, divergence_scores_df, tag: str) -> None:
#     """ Add distance score to the top links and save
#     """
#     if not os.path.isdir(TOP_LINKS_WITH_DISTANCE_SCORES_DIR):
#         os.makedirs(TOP_LINKS_WITH_DISTANCE_SCORES_DIR)
#     # - add distance scores
#     top_links_with_scores = divergence_scores_df.merge(top_links, on=['Regulator', 'Target'])
#     file_name = Path(TOP_LINKS_WITH_DISTANCE_SCORES_DIR) / f'{tag}_{model_name}.csv'
#     top_links_with_scores.to_csv(file_name, index=False)

def shortlist(model_name, studies, divergence_scores_df, top_percentile):
    links_studies = retreive_links_with_genenames(model_name, studies)
    # - assert
    assert (len(links_studies[0]) == len(divergence_scores_df))
    # - calculate abs difference between ctr and sample
    divergence_scores_df['difference']  = np.abs(links_studies[0]['Weight']-links_studies[1]['Weight'])
    # - choose the top percentile based on difference
    divergence_scores_short = choose_top_quantile(divergence_scores_df, top_percentile, col_name='difference')
    divergence_scores_short.drop(columns=['difference'], inplace=True)
    return divergence_scores_short


if __name__ == '__main__':
    # - parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--force', type=bool, default=False,
                        help='To overwrite if the result files exist')
    parser.add_argument('--top_quantile', type=float, default=0.9,
                        help='Top percentile of the links with the largest divergence')
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    force = args.force
    top_quantile = args.top_quantile
    # - run noise analysis to determine distance scores from ctr to sample
    for model_name in F_selected_models():
        network_noise_analysis_obj = NetworkNoiseAnalysis(model_name=model_name, studies=studies,
                                                          base_dir=NETWORK_NOISE_DIR, force=force)
        # - run the analysis and post analysis
        network_noise_analysis_obj.run_decorator()
        # - calculate overall distance scores for this model by accounting for different noise types
        divergence_scores_df = network_noise_analysis_obj.determine_overall_divergence_scores()
        # - shortlist it to only top percentile change
        divergence_scores_short = shortlist(model_name, studies, divergence_scores_df, top_quantile)
        # - save
        network_noise_analysis_obj.save_divergence_scores_df(model_name, divergence_scores_short)