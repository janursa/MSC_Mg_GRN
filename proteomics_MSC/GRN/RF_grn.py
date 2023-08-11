import sys
import os
import numpy as np
import yaml
from importlib.resources import open_text

from proteomics_MSC.imports import GRN_DIR, CALIBRATION_DIR
from proteomics_MSC.common_tools import read_write_data, calibration, F_DE_data
from proteomics_MSC.common_tools.links import batch_run_generni


with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    studies = config['studies']
    n_replica = config['n_replica_grn_run_RF']
    time_points = config['time_points']

    method = 'RF'
    param = dict(estimator_t=method)
    if not os.path.isdir(os.path.join(GRN_DIR, method)):
        os.makedirs(os.path.join(GRN_DIR, method))

    for DE_type, DE_data in F_DE_data().items():
        DE_data = F_DE_data()[DE_type]
        for study in studies:
            print(f'----------------------{DE_type}, {study} ---------------')
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values.tolist()
            data = read_write_data(mode='read', tag=f'{DE_type}_{study}')
            n_timepoints = data.shape[0]
            days = time_points[0:n_timepoints]
            specs = dict(
                i_start=0,
                i_end=n_replica,
                param=param,
                gene_names=gene_names,
                time_points=days,
                method=method,
                output_dir=GRN_DIR,
                verbose=False,

            )

            data = np.asarray(data)
            print(f'Data shape: (n_samples_time_series, n_genes) = {data.shape}')
            _, param_unique = calibration.retrieve_data(study=study, method=method, DE_type=DE_type, output_dir=CALIBRATION_DIR)
            batch_run_generni(study=study, data=data, DE_type=DE_type, param_unique=param_unique, **specs)
