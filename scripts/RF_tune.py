"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
from imports import *

if __name__ == '__main__':
    method = 'RF'
    param_grid = param_grid_RF
    test_size = 0

    istart = 50
    iend = 100
    # - read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
    data_mg = utils.process_data(df_target, study='mg', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')

    # - define settings and tune
    param = dict(estimator_t=method)
    specs = dict(
        #     train_flag=True # maximizing training data
    )
    utils.calibration.batch_tune('ctr', method, data_ctr, param, protnames, test_size, time, param_grid, specs,
                                 istart=istart, iend=iend, OUTPUT_DIR=OUTPUT_DIR)
    utils.calibration.batch_tune('mg', method, data_mg, param, protnames, test_size, time, param_grid, specs,
                                 istart=istart, iend=iend, OUTPUT_DIR=OUTPUT_DIR)

