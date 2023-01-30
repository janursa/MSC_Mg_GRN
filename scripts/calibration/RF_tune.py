"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import df_target, time_points, param_grid_RF, protnames, CALIBRATION_DIR
from scripts.utils import process_data, calibration

if __name__ == '__main__':
    method = 'RF'
    param_grid = param_grid_RF()
    test_size = 0
    cv=None
    # - read the data
    data_ctr = process_data(df_target(), study='ctr', time_points=time_points(), standardize=False)
    data_mg = process_data(df_target(), study='mg', time_points=time_points(), standardize=False)
    # data_combined = process_data(df_target, study='combined', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    # - define settings and calibration
    param = dict(estimator_t=method)
    specs = dict(
        i_start=0,
        i_end=10,
        param=param,
        gene_names=protnames(),
        time_points=time_points(),
        test_size=test_size,
        cv=cv,
        method=method,
        param_grid=param_grid,
        output_dir=CALIBRATION_DIR,
        random_state=0,
        n_jobs=10,
    )

    calibration.calibrate(study='ctr', data=data_ctr, **specs)
    calibration.calibrate(study='mg', data=data_mg, **specs)