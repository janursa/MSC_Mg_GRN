import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import time_points, F_DE_data, GRN_DIR, CALIBRATION_DIR
from scripts.utils import process_data, calibration
from scripts.utils.links import batch_GRN


if __name__ == '__main__':
    method = 'RF'
    param = dict(estimator_t=method)
    if not os.path.isdir(os.path.join(GRN_DIR, method)):
        os.makedirs(os.path.join(GRN_DIR, method))

    # studies = ['ctr', 'mg', 'all-in']
    studies = [ 'ctr', 'mg']
    # studies = ['mg']
    DE_types = ['combined_30', 'combined_50']
    for DE_type in DE_types:
    # for DE_type, DE_data in F_DE_data().items():
        DE_data = F_DE_data()[DE_type]
        for study in studies:
            print(f'----------------------{DE_type}, {study} ---------------')
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values.tolist()
            specs = dict(
                i_start=50,
                i_end=100,
                param=param,
                gene_names=gene_names,
                time_points=time_points(),
                method=method,
                output_dir=GRN_DIR,
                verbose=False,

            )
            data = process_data(DE_data, study=study, time_points=time_points(), standardize=False)
            data = np.asarray(data)
            print(f'Data shape: (n_samples_time_series, n_genes) = {data.shape}')
            _, param_unique = calibration.retrieve_data(study=study, method=method, DE_type=DE_type, output_dir=CALIBRATION_DIR)
            batch_GRN(study=study, data=data, DE_type=DE_type, param_unique=param_unique, **specs)
