"""
    Evaluate different models using conventional approaches such as precision recall
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Callable, Optional, TypeAlias, Any
import json


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from imports import ENRICH_DIR, CALIBRATION_DIR, MODELSELECTION_DIR, F_DE_data, GRN_DIR, RANDOM_MODELS_DIR, top_quantiles
from utils import calibration
from md_aux import calculate_early_precision
from md_aux import score_type, get_baseline_scores, save_scores


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
        all_scores['_'.join([DE_type, 'baseline'])] = {'ep':ep_scores_random,
                                                       'ep_series': [list(item) for item in ep_scores_series_random],
                                                       'epr':list(ep_scores_random/np.mean(ep_scores_random)),
                                                       'percentile_rank':None, 'sig_flag':None}
        # - get the golden links
        golden_links = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
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
            s, p = scipy.stats.ttest_1samp(ep_scores_random, ep_score)
            if (p < 0.05) & (ep_score > np.mean(ep_scores_random)):
                sig_flag = True
            else:
                sig_flag = False
            # - get the test score
            test_scores_method = retrieve_prediction_scores(DE_type, method, studies)  # test score is calculated for both ctr and sample
            all_scores[model_name] = {'R2':test_scores_method,
                                      'ep_series':ep_scores,
                                      'epr':ep_score/np.mean(ep_scores_random),
                                      'sig_flag':sig_flag,
                                      'percentile_rank':percentile_rank}

    return all_scores
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--GRN_methods', nargs='+', default=['RF', 'ridge', 'portia'])
    parser.add_argument('--target_study', type=str, default='ctr')
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    args, remaining_args = parser.parse_known_args()

    GRN_methods = args.GRN_methods
    target_study = args.target_study
    studies = args.studies

    scores = calculate_scores() # model_name: scores[dict]
    save_scores(scores)
