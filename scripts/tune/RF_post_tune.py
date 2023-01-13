"""
Reads the results of hyper param tuning (RF), 
pool them (because there are multiple runs), average them, and output them to files.
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *

if __name__ == '__main__':
    utils.calibration.pool_oo('ctr', 'RF', n=100, OUTPUT_DIR=OUTPUT_DIR)
    utils.calibration.pool_oo('mg', 'RF', n=100, OUTPUT_DIR=OUTPUT_DIR)

    utils.calibration.plot_oo('RF', param_grid_RF, protnames, OUTPUT_DIR)
    # plt.show()

