"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *

if __name__ == '__main__':
    method = 'ridge'
    param_grid = param_grid_ridge
    # - read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=True)
    data_mg = utils.process_data(df_target, study='mg', standardize=True)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #- create the dataaset
    # test_size = 1 / len(time)
    test_size = 0
    cv=5
    #- define settings and calibration
    param = dict(estimator_t = method)
    specs = dict(
        i_start=0,
        i_end=10,
        param=param,
        gene_names=protnames,
        time_points=time,
        test_size=test_size,
        cv=cv,
        method=method,
        param_grid=param_grid,
        OUTPUT_DIR=OUTPUT_DIR,
        random_state=0,
        n_jobs=10,
    )

    utils.calibration.main(study='ctr', data=data_ctr, **specs)
    # utils.calibration.main(study='mg', data=data_mg, **specs)
