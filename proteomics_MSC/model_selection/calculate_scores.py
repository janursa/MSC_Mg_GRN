"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import pandas as pd
import argparse
from typing import Dict, List, Tuple, Optional


from proteomics_MSC.imports import CALIBRATION_DIR,GRN_DIR, STRING_DIR
from proteomics_MSC.common_tools import F_DE_data, calibration
from proteomics_MSC.common_tools.model_selection import score_type, get_baseline_scores, save_scores, determine_sig_flag, calculate_early_precision

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

def retrieve_prediction_scores(DE_type: str, method: str, studies: List[str]) -> Optional[Tuple[float, float]]:
    """
    Retrieves calibration test scores for ridge and RF and sets portia's score to None
    """
    if method == 'portia':
        scores_studies = None
    else:
        scores_studies = []
        for study in studies:
            best_scores, _ = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
            mean_score = np.mean(best_scores)
            scores_studies.append(mean_score)
    return scores_studies

def calculate_scores() -> Dict[str, score_type]:
    """Calculates scores for all models including random models
    ep, epr
    """
    all_scores: Dict[str, score_type] = {} # model_name: dict of different scores
    for DE_type in F_DE_data().keys():
        # - get the baseline scores
        ep_scores_random, ep_scores_series_random = get_baseline_scores(DE_type) #1000
        print('mean: ', np.mean(ep_scores_random), 'max: ', np.max(ep_scores_random))
        all_scores['_'.join([DE_type, 'baseline'])] = {'ep':ep_scores_random,
                                                       'ep_series': [list(item) for item in ep_scores_series_random],
                                                       'epr':list(ep_scores_random/np.mean(ep_scores_random)),
                                                       'percentile_rank':None, 'sig_flag':None}
        # - get the golden links
        golden_links = pd.read_csv(os.path.join(STRING_DIR, f'links_{DE_type}.csv'), index_col=False)
        #- scores of all DE_type
        for method in GRN_methods:
            # - model name
            model_name = '_'.join([DE_type, method])
            # - get the links
            links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{target_study}.csv'), index_col=False)
            # - calculate epr for each model
            ep_scores = calculate_early_precision(links, golden_links, top_quantiles)
            ep_score = np.mean(ep_scores)
            # - calculate what percentage of random values are bigger than random score
            count = sum(1 for x in ep_scores_random if
                        x > ep_score)  # Count how many elements in 'a' are bigger than 'b'
            percentile_rank = int((count / len(ep_scores_random)) * 100)  # Calculate the percentile_rank
            # - check if AUC values are sig larger than random values
            sig_flag = determine_sig_flag(ep_scores_random, ep_score)
            # - get the test score
            test_scores_method = retrieve_prediction_scores(DE_type, method, studies)  # test score is calculated for both ctr and sample
            all_scores[model_name] = {'R2':test_scores_method,
                                      'ep_series':ep_scores,
                                      'epr':ep_score/np.mean(ep_scores_random),
                                      'sig_flag':sig_flag,
                                      'percentile_rank':percentile_rank}

    return all_scores
if __name__ == '__main__':
    GRN_methods =config['GRN_methods']
    studies = config['studies']
    target_study = config['model_selection_study']
    top_quantiles = config['top_quantiles']

    scores = calculate_scores() # model_name: scores[dict]
    save_scores(scores)
