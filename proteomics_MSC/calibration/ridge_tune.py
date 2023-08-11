import sys
import os
import numpy as np
import yaml
from importlib.resources import open_text

from proteomics_MSC.imports import CALIBRATION_DIR
from proteomics_MSC.common_tools import F_DE_data, F_DE_protnames
from proteomics_MSC.common_tools.calibration import param_grid_ridge
from proteomics_MSC.common_tools import read_write_data, calibration

with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    studies = config['studies']
    time_points = config['time_points']

    method = 'ridge'

    param_grid = param_grid_ridge()
    for DE_type, DE_data in F_DE_data().items():
        for study in studies:
            # - read the data
            data = data = read_write_data(mode='read', tag=f'{DE_type}_{study}')

            n_timepoints = data.shape[0]
            days = time_points[0:n_timepoints]
            # print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
            #- create the dataaset
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values
            #- define settings and calibration
            param = dict(estimator_t = method)
            specs = dict(
                i_start=0,
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
            to_save = os.path.join(CALIBRATION_DIR, method, DE_type, f'best_scores_{study}.npy')
            np.save(to_save, best_scores)
            print(f'output -> {to_save}')
            to_save = os.path.join(CALIBRATION_DIR, method, DE_type, f'best_params_{study}.npy')
            np.save(to_save, best_params)
            print(f'output -> {to_save}')