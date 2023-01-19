"""
Reads the results of hyper param tuning (RF), 
pool them (because there are multiple runs), average them, and output them to files.
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *
from utils import calibration

if __name__ == '__main__':
    # replica_n = 5
    # utils.calibration.pool_oo('ctr', 'ridge', n=replica_n,OUTPUT_DIR=OUTPUT_DIR)
    # utils.calibration.pool_oo('mg', 'ridge', n=replica_n,OUTPUT_DIR=OUTPUT_DIR)
    calibration.plot_oo('ridge', param_grid_ridge, protnames, OUTPUT_DIR)

