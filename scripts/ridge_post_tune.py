"""
Reads the results of hyper param tuning (RF), 
pool them (because there are multiple runs), average them, and output them to files.
"""
from imports import *

if __name__ == '__main__':
    utils.calibration.pool_oo('ctr', 'ridge', n=100,OUTPUT_DIR=OUTPUT_DIR)
    utils.calibration.pool_oo('mg', 'ridge', n=100,OUTPUT_DIR=OUTPUT_DIR)

    utils.calibration.plot_oo('ridge', param_grid_ridge, protnames, OUTPUT_DIR)

