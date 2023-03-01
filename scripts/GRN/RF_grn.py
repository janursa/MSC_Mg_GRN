import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import time_points, F_DE_data, GRN_DIR, CALIBRATION_DIR
from utils import process_data, calibration
from utils.links import batch_GRN


if __name__ == '__main__':
    method = 'RF'
    param = dict(estimator_t=method)
    if not os.path.isdir(os.path.join(GRN_DIR, method)):
        os.makedirs(os.path.join(GRN_DIR, method))

    studies = [ 'ctr', 'mg']
    studies = ['mg']
    for DE_type, DE_data in F_DE_data().items():
        if 'day1_11_KNN' not in DE_type: #TODO: remove this
            continue
        DE_data = F_DE_data()[DE_type]
        for study in studies:
            print(f'----------------------{DE_type}, {study} ---------------')
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values.tolist()
            data = process_data(DE_data, study=study, standardize=False)
            n_timepoints = data.shape[0]
            days = time_points()[0:n_timepoints]
            specs = dict(
                i_start=100,
                i_end=500,
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
            batch_GRN(study=study, data=data, DE_type=DE_type, param_unique=param_unique, **specs)
