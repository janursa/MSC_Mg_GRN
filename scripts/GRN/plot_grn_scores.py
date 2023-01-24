"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import OUTPUT_DIR
from utils.links import retrieve_scores
from utils.calibration import plot_scores


def main(method):
    # - retreive scores
    trainscores_ctr, testscores_ctr = retrieve_scores(method=method, study='ctr', output_dir=OUTPUT_DIR)
    trainscores_mg, testscores_mg = retrieve_scores(method=method, study='mg', output_dir=OUTPUT_DIR)

    # utils.links.plot_scores(data_ctr=[trainscores_ctr, trainscores_mg], data_sample=[testscores_ctr, testscores_mg])
    fig = plot_scores(data_ctr=testscores_ctr, data_sample=testscores_mg, ylabel='Train score')
    fig.savefig(os.path.join(OUTPUT_DIR, 'GRN', method, 'testscores.pdf'))
    fig = plot_scores(data_ctr=trainscores_ctr, data_sample=trainscores_mg, ylabel='Test score')
    fig.savefig(os.path.join(OUTPUT_DIR, 'GRN', method, 'trainscores.pdf'))

if __name__ == '__main__':
    main(method='RF')
    main(method='ridge')
