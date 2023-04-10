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

from imports import F_DE_data, time_points, param_grid_RF, F_DE_protnames, CALIBRATION_DIR
from utils import read_write_data, calibration

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--n_replica', type=int, default=10)
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    n_replica = args.n_replica

    method = 'RF'

    # - define settings and calibration
    param = dict(estimator_t=method)
    for DE_type, DE_data in F_DE_data().items():
        param_grid = param_grid_RF(len(DE_data['Protein']))
        for study in studies:
            # - read the data
            data = read_write_data(mode='read', tag=f'{DE_type}_{study}')
            # print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
            n_timepoints = data.shape[0]
            days = time_points()[0:n_timepoints]
            #- create the dataaset
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist()+['mg']
            else:
                gene_names = DE_data['Protein'].values.tolist()
            #- define settings and calibration
            specs = dict(
                i_start=n_replica,
                i_end=n_replica,
                param=param,
                gene_names=gene_names,
                time_points=days,
                method=method,
                param_grid=param_grid,
                output_dir=os.path.join(CALIBRATION_DIR, method, DE_type, study),
                random_state=0,
                n_jobs=1,
            )
            best_scores, best_params = calibration.calibrate(study=study, data=data, **specs)
            np.save(os.path.join(CALIBRATION_DIR, method, DE_type, f'best_scores_{study}.npy'), best_scores)
            np.save(os.path.join(CALIBRATION_DIR, method, DE_type, f'best_params_{study}.npy'), best_params)