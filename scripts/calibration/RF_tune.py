"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import F_DE_data, time_points, param_grid_RF, F_DE_protiens, CALIBRATION_DIR
from scripts.utils import process_data, calibration

if __name__ == '__main__':
    method = 'RF'

    # - define settings and calibration
    param = dict(estimator_t=method)
    studies = ['ctr', 'mg', 'all-in']
    for DE_type, DE_data in F_DE_data().items():
        # if DE_type in ['early_30', 'early_50', 'late_30']: #completed ones
        #     continue
        param_grid = param_grid_RF(len(DE_data['Protein']))
        for study in studies:
            # - read the data
            data = process_data(DE_data, study=study, time_points=time_points(), standardize=False)
            print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
            #- create the dataaset
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist()+['mg']
            else:
                gene_names = DE_data['Protein'].values.tolist()
            #- define settings and calibration
            specs = dict(
                i_start=0,
                i_end=10,
                param=param,
                gene_names=gene_names,
                time_points=time_points(),
                method=method,
                param_grid=param_grid,
                output_dir=os.path.join(CALIBRATION_DIR, method, DE_type, study),
                random_state=0,
                n_jobs=10,
            )
            best_scores, best_params = calibration.calibrate(study=study, data=data, **specs)
            np.save(os.path.join(CALIBRATION_DIR, method, DE_type, f'best_scores_{study}.npy'), best_scores)
            np.save(os.path.join(CALIBRATION_DIR, method, DE_type, f'best_params_{study}.npy'), best_params)