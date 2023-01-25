"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *
from utils import process_data, calibration
if __name__ == '__main__':
    method = 'ridge'
    param_grid = param_grid_ridge
    # - read the data
    data_ctr = process_data(df_target, study='ctr', time_points=time_points(), standardize=False)
    data_mg = process_data(df_target, study='mg', time_points=time_points(), standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #- create the dataaset
    # test_size = 1 / len(time)
    test_size = 0
    cv=5
    #- define settings and calibration
    param = dict(estimator_t = method)
    specs = dict(
        i_start=0,
        i_end=1,
        param=param,
        gene_names=protnames,
        time_points=time_points(),
        test_size=test_size,
        cv=cv,
        method=method,
        param_grid=param_grid,
        OUTPUT_DIR=OUTPUT_DIR,
        random_state=0,
        n_jobs=1,
        loo=True, #leave one out
    )

    calibration.calibrate(study='ctr', data=data_ctr, **specs)
    calibration.calibrate(study='mg', data=data_mg, **specs)
