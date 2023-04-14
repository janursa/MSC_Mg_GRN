"""
    Filter models based on R2 score and sig flags. Then, use epr score to choose the best performing models for early and late
"""
import sys
import os
import numpy as np
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from md_aux import retreieve_scores
from imports import MODELSELECTION_DIR

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args, remaining_args = parser.parse_known_args()

    # - retreive scores
    scores = retreieve_scores()
    # - non random models
    no_random_scores = {key:value for key, value in scores.items() if 'baseline' not in key}
    # - filter based on R2 score: it should be bigger than 0 across both ctr and sample
    filtered_scores_R2 = {key:value for key, value in no_random_scores.items() if ((value['R2'] is None) or all(score > 0 for score in value['R2']))}
    # - filter based on sig epr
    filtered_scores_R2_sig = {key:value for key, value in filtered_scores_R2.items() if value['sig_flag']}
    # - output shortlisted model
    shortlisted_modelnames = list(filtered_scores_R2_sig.keys())
    np.savetxt(os.path.join(MODELSELECTION_DIR, f'shortlisted_models.txt'), shortlisted_modelnames, delimiter=",",
               fmt="%s")
    # - select top performing model for early and late
    shortlisted_modelnames_early = list(filter(lambda name: name.split('_')[0]=='early', shortlisted_modelnames))
    shortlisted_modelnames_late = list(filter(lambda name: name.split('_')[0]=='late', shortlisted_modelnames))

    shortlisted_epr_scores_early = [value['epr'] for key,value in filtered_scores_R2_sig.items() if key in shortlisted_modelnames_early]
    shortlisted_epr_scores_late = [value['epr'] for key,value in filtered_scores_R2_sig.items() if key in shortlisted_modelnames_late]

    best_model_names = [shortlisted_modelnames_early[np.argmax(shortlisted_epr_scores_early)],
                        shortlisted_modelnames_late[np.argmax(shortlisted_epr_scores_late)]]
    for model_name in best_model_names:
        print(f"{model_name} -> ep scores: {np.mean(scores[model_name]['ep_series'])}")
    np.savetxt(os.path.join(MODELSELECTION_DIR, f'selected_models.txt'), best_model_names, delimiter=",", fmt="%s")
