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
from utils.links import pool_GRN_oo, retreive_scores

if __name__ == '__main__':
    method = 'RF'
    replica_n = 10
    # - pool data
    pool_GRN_oo(method=method, study='ctr', replica_n=replica_n, output_dir=OUTPUT_DIR)
    pool_GRN_oo(method=method, study='mg', replica_n=replica_n, output_dir=OUTPUT_DIR)
    print('GRN results are successfully pooled for RF')