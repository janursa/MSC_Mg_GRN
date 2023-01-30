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

from scripts.imports import GRN_DIR
from scripts.utils.links import retrieve_scores
from scripts.utils.calibration import plot_scores


def main(method):
    for study in ['ctr', 'mg']:
        # - retreive scores
        trainscores, testscores = retrieve_scores(method=method, study=study, output_dir=GRN_DIR)

        fig = plot_scores(data_ctr=trainscores, data_sample=testscores, ylabel='Train score')
        fig.savefig(os.path.join(GRN_DIR, method, 'trainscores.pdf'))

        fig = plot_scores(data_ctr=testscores, data_sample=testscores, ylabel='Test score')
        fig.savefig(os.path.join(GRN_DIR, method, 'testscores.pdf'))


if __name__ == '__main__':
    main(method='RF')
    main(method='ridge')
