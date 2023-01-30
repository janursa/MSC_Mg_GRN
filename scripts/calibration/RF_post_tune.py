"""
Reads the results of hyper param tuning (RF), 
pool them (because there are multiple runs), average them, and output them to files.
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import param_grid_RF, protnames, CALIBRATION_DIR
from scripts.utils import calibration

if __name__ == '__main__':
    calibration.plot_oo('RF', param_grid_RF(), protnames(), CALIBRATION_DIR)

