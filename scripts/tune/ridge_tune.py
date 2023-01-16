"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *

if __name__ == '__main__':
    # - read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=True)
    data_mg = utils.process_data(df_target, study='mg', standardize=True)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #- create the dataaset
    test_size = 1 / len(time)
    #- define settings and tune
    param = dict(estimator_t = 'ridge')
    specs = dict(
        random_state = 0,
        n_jobs=1,
        n_repeat=5
        #     train_flag=True # maximizing training data
    )
    istart = 0
    iend = 5

    utils.calibration.batch_tune('ctr', 'ridge', data_ctr, param,protnames,test_size,time, param_grid_ridge, istart=istart, iend=iend, OUTPUT_DIR=OUTPUT_DIR, **specs)
    utils.calibration.batch_tune('mg', 'ridge', data_mg, param,protnames,test_size,time, param_grid_ridge, istart=istart, iend=iend, OUTPUT_DIR=OUTPUT_DIR, **specs)

with open(os.path.join(output_dir, 'sampled_permts.txt'), 'w') as f:
    print({'permts': sampled_permts}, file=f)