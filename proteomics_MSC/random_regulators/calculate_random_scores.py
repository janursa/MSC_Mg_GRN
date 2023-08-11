import argparse
import os
import sys
from typing import List, Tuple
from pathlib import Path

import numpy as np
from tqdm import tqdm


from proteomics_MSC.preprocess.process_DE_analysis import load_imputated_df
from proteomics_MSC.common_tools.links import run_portia
from proteomics_MSC.common_tools import process_data, F_selected_models
from proteomics_MSC.common_tools.model_selection import string_functions_analysis, calculate_early_precision
from proteomics_MSC.imports import RANDOM_REGULATORS_DIR

def random_regulators_selection(imput_method:str, period:str, study:str, n_features:int, n_repeat:int) -> Tuple[List[np.ndarray], List[List[str]]]:
    """Randomly select proteins as regulators and create corrosponding dataset"""
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

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':

    periods = config['periods']
    n_features = config['n_random_regulators']
    top_quantiles = config['top_quantiles']
    n_repeat = 100

    correct_names = {'early':'day1_11', 'late':'day1_21'} # imputated data are in day* format
    study = 'ctr' # calculate early precision on ctr

    os.makedirs(RANDOM_REGULATORS_DIR, exist_ok=True)
    # - load names of selected models
    selected_models = F_selected_models()
    for selected_model in selected_models:
        period, imput_method = selected_model.split('_')[0:2]
        period = correct_names[period]
        # - get shortlisted data and genenames
        data_stack, genenames_stack = random_regulators_selection(imput_method, period, study, n_features, n_repeat)
        assert (len(genenames_stack) == n_repeat)
        # - run grn and get links
        links_stack = []
        for data, genenames in tqdm(zip(data_stack, genenames_stack), desc='Run portia'):
            links = run_portia(data, genenames)
            links_stack.append(links)
        # - calculate early precision for each sample of random regulators
        golden_links = []
        ep_score_stack = []
        for data_i, genenames in tqdm(enumerate(genenames_stack), desc='Calculate the scores'):
            string_links = string_functions_analysis(genenames)
            golden_links.append(string_links)
            ep_scores = calculate_early_precision(links_stack[data_i], string_links, top_quantiles)
            ep_score_stack.append(np.mean(ep_scores))
        to_save = Path(RANDOM_REGULATORS_DIR) / f'ep_scores_random_{selected_model}.csv'
        with open(to_save, 'w') as ff:
            np.savetxt(ff, ep_score_stack, delimiter=',')
        print(f'output -> {to_save}')





