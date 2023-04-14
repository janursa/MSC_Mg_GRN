import argparse
import os
from typing import List, Tuple, Optional
from pathlib import Path

import numpy as np

from post_statistical_analysis.psa_aux import load_imputated_df
from utils.links import run_portia
from utils import process_data
from utils.enrich_analysis import string_network_analysis
from imports import top_quantiles, F_selected_models, RANDOM_FEATURE_DIR
from model_selection.md_aux import calculate_early_precision
def random_feature_selection(imput_method:str, period:str, study:str, n_features:int, n_repeat:int) -> Tuple[List[np.ndarray], List[List[str]]]:
    # - load the imputed df
    df_imput = load_imputated_df(imput_method, period)
    # - shortlist df by choosing n_features randomly
    data_stack = []
    genenames_stack = []
    for i in range(n_repeat):
        df_shortlist = df_imput.sample(n=n_features, ignore_index=True)
        data = process_data(df_shortlist, study)
        data_stack.append(data)
        genenames_stack.append(df_shortlist['Protein'].tolist())

    return data_stack, genenames_stack


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--periods', nargs='+', default=['early', 'late'])
    parser.add_argument('--n_features', type=int, default=50, help="Number of features (DE proteins) to select randomly")
    parser.add_argument('--n_repeat', type=int, default=100, help="How many random selection")
    args, remaining_args = parser.parse_known_args()

    periods = args.periods
    n_features = args.n_features
    n_repeat = args.n_repeat

    correct_names = {'early':'day1_11', 'late':'day1_21'} # imputated data are in day* format
    study = 'ctr' # calculate early precision on ctr

    # - create directories
    if not os.path.isdir(RANDOM_FEATURE_DIR):
        os.makedirs(RANDOM_FEATURE_DIR)

    # - load names of selected models
    selected_models = F_selected_models()
    for selected_model in selected_models:
        period, imput_method = selected_model.split('_')[0:2]
        period = correct_names[period]

        # - get shortlisted data and genenames
        data_stack, genenames_stack = random_feature_selection(imput_method, period, study, n_features, n_repeat)
        assert (len(genenames_stack) == n_repeat)
        # - run grn and get links
        links_stack = []
        for data, genenames in zip(data_stack, genenames_stack):
            links = run_portia(data, genenames)
            links_stack.append(links)
        # - run string enquiry and get
        golden_links = []
        ep_score_stack = []
        for data_i, genenames in enumerate(genenames_stack):
            string_links = string_network_analysis(genenames)
            golden_links.append(string_links)
            print(data_i)

            ep_scores = calculate_early_precision(links_stack[data_i], string_links, top_quantiles)
            ep_score_stack.append(np.mean(ep_scores))
        with open(Path(RANDOM_FEATURE_DIR)/f'ep_scores_random_{selected_model}.csv', 'w') as ff:
            np.savetxt(ff, ep_score_stack, delimiter=',')





