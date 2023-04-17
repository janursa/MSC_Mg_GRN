"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
import sys
import os
import numpy as np
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_DE_data, time_points, param_grid_ridge, CALIBRATION_DIR
from common_tools import read_write_data, calibration
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    args, remaining_args = parser.parse_known_args()

    studies = args.studies

    method = 'ridge'

    param_grid = param_grid_ridge()
    for DE_type, DE_data in F_DE_data().items():
        for study in studies:
            # - read the data
            data = data = read_write_data(mode='read', tag=f'{DE_type}_{study}')

            n_timepoints = data.shape[0]
            days = time_points()[0:n_timepoints]
            # print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
            #- create the dataaset
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values
            #- define settings and calibration
            param = dict(estimator_t = method)
            specs = dict(
                i_start=1,
                i_end=1,
                param=param,
                gene_names=gene_names,
                time_points=days,
                loo=True,
                method=method,
                param_grid=param_grid,
                output_dir= os.path.join(CALIBRATION_DIR, method, DE_type, study),
                random_state=0,

                n_jobs=10,
            )
            best_scores, best_params = calibration.calibrate(study=study, data=data, **specs)
            np.save(os.path.join(CALIBRATION_DIR, method, DE_type, f'best_scores_{study}.npy'), best_scores)
            np.save(os.path.join(CALIBRATION_DIR, method, DE_type, f'best_params_{study}.npy'), best_params)