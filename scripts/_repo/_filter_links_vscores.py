"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, F_DE_protnames, CALIBRATION_DIR
from scripts.utils import calibration

if __name__ == '__main__':
    studies = ['ctr', 'mg']
    methods = ['RF', 'Ridge']

    for DE_type, _ in F_DE_protnames().items():
        protnames = F_DE_protnames()[DE_type]
        for method in methods:
            dir_method_filtered = os.path.join(GRN_DIR, f'{method}_filtered')
            if not os.path.isdir(dir_method_filtered):
                os.makedirs(dir_method_filtered)

            for study in studies:
                links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=[0])
                best_scores, best_params = calibration.retrieve_data(study, method, DE_type, CALIBRATION_DIR)
                protnames_left = np.asarray(protnames)[np.asarray(best_scores > 0)] #TODO: needs revision
                filter_mask = links['Target'].isin(protnames_left)
                links_f = links.loc[filter_mask,:]

                links_f.to_csv(os.path.join(dir_method_filtered, f'links_{DE_type}_{study}.csv'), index=False)
