"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *

def add_score_to_links(links, scores, protnames):
    links['FitScore'] = None
    for target, score in zip(protnames, scores):
        links.loc[links['Target'] == target, 'FitScore'] = score
    return links

if __name__ == '__main__':
    method = 'RF'
    replica_n = 10
    # - pool data
    utils.links.pool_data(method=method, study='ctr', replica_n=replica_n, output_dir=OUTPUT_DIR)
    utils.links.pool_data(method=method, study='mg', replica_n=replica_n, output_dir=OUTPUT_DIR)
    #- retreive scores
    trainscores_ctr, testscores_ctr = utils.links.retreive_scores(method=method, study='ctr', output_dir=OUTPUT_DIR)
    trainscores_mg, testscores_mg = utils.links.retreive_scores(method=method, study='mg', output_dir=OUTPUT_DIR)

    # utils.links.plot_scores(data_ctr=[trainscores_ctr, trainscores_mg], data_sample=[testscores_ctr, testscores_mg])
    fig = utils.calibration.plot_scores(data_ctr=testscores_ctr, data_sample=testscores_mg, ylabel='Training score')
    fig.savefig(os.path.join(OUTPUT_DIR, 'GRN', method, 'oobscores.pdf'))
    fig = utils.calibration.plot_scores(data_ctr=trainscores_ctr, data_sample=trainscores_mg, ylabel='OOB score')
    fig.savefig(os.path.join(OUTPUT_DIR, 'GRN', method, 'trainscores.pdf'))
    #
    plt.show()
