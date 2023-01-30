"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import OUTPUT_DIR
from utils.links import pool_GRN_oo


if __name__ == '__main__':
    method = 'RF'
    replica_n = 10

    for study in ['ctr', 'mg', 'combined']:
        pool_GRN_oo(method=method, study=study, replica_n=replica_n, output_dir=OUTPUT_DIR)
    print('GRN results are successfully pooled for RF')