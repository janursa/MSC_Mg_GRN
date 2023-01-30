import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import df_target, time_points, protnames, GRN_DIR, CALIBRATION_DIR
from scripts.utils import process_data, calibration
from scripts.utils.links import batch_GRN


if __name__ == '__main__':
    method = 'RF'
    param = dict(estimator_t=method)

    # studies = ['combined']
    studies = ['ctr', 'mg', 'combined']
    for study in studies:
        if study == 'combined':
            gene_names = protnames()+['mg']
        else:
            gene_names = protnames()
        specs = dict(
            i_start=0,
            i_end=10,
            param=param,
            gene_names=gene_names,
            time_points=time_points(),
            test_size=0,
            method=method,
            output_dir=GRN_DIR,
        )
        data = process_data(df_target(), study=study, time_points=time_points(), standardize=False)
        data = np.asarray(data)
        print(f'Data shape: (n_samples_time_series, n_genes) = {data.shape}')
        _, param_unique = calibration.retrieve_data(study=study, method=method, output_dir=CALIBRATION_DIR)
        batch_GRN(study=study, data=data, param_unique=param_unique, **specs)
