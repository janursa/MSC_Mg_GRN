"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *

if __name__ == '__main__':
    method = 'RF'
    param_grid = param_grid_RF
    test_size = 0

    istart = 1
    iend = 5
    # - read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
    data_mg = utils.process_data(df_target, study='mg', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    # - define settings and tune
    param = dict(estimator_t=method)
    specs = dict(
        random_state=0
        #     train_flag=True # maximizing training data
    )
    utils.calibration.batch_tune('ctr', method, data_ctr, param, protnames, test_size, time, param_grid, specs,
                                 istart=istart, iend=iend, OUTPUT_DIR=OUTPUT_DIR)
    utils.calibration.batch_tune('mg', method, data_mg, param, protnames, test_size, time, param_grid, specs,
                                 istart=istart, iend=iend, OUTPUT_DIR=OUTPUT_DIR)

