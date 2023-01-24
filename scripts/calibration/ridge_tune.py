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

    for study in ['ctr', 'mg', 'combined']:
        data = process_data(df_target, study=study, time_points=time_points(), standardize=False)
        print(f'Data shape: (n_samples_time_series, n_genes) = {data.shape}')
        param = dict(estimator_t=method, alpha=1)
        specs = dict(
            i_start=0,
            i_end=5,
            param=param,
            gene_names=protnames,
            time_points=time_points(),
            test_size=0,
            cv=5,
            method=method,
            param_grid=param_grid,
            OUTPUT_DIR=OUTPUT_DIR,
            random_state=0,
            n_jobs=1,
        )

        calibration.main(study=study, data=data, **specs)
