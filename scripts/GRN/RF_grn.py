import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *

if __name__ == '__main__':
    #- read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
    data_mg = utils.process_data(df_target, study='mg', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #-
    method='RF'
    test_size=0
    param = dict(estimator_t=method)
    specs = dict(
        i_start=0,
        i_end=10,
        param=param,
        gene_names=protnames,
        time_points=time,
        test_size=test_size,
        method=method,
        output_dir=os.path.join(OUTPUT_DIR,'GRN'),
    )

    _, param_unique_ctr = utils.calibration.retreive_data(study='ctr', method=method, OUTPUT_DIR=OUTPUT_DIR)
    _, param_unique_mg = utils.calibration.retreive_data(study='mg', method=method, OUTPUT_DIR=OUTPUT_DIR)

    # - network inference
    utils.links.batch_GRN(study='ctr', data=data_ctr, param_unique=param_unique_ctr, **specs)
    utils.links.batch_GRN(study='mg', data=data_mg, param_unique=param_unique_mg, **specs)
