import sys
import os
import numpy as np
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import time_points, F_DE_data, GRN_DIR, CALIBRATION_DIR
from utils import read_write_data, calibration
from utils.links import batch_run_generni


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--n_replica', type=int, default=100)
    args, remaining_args = parser.parse_known_args()

    studies = args.studies
    n_replica = args.n_replica

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
            days = time_points()[0:n_timepoints]
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
