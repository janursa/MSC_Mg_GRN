from imports import *

if __name__ == '__main__':
    #- read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
    data_mg = utils.process_data(df_target, study='mg', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #-
    method='RF'
    param = dict(estimator_t=method)
    test_size = 0
    istart = 54
    iend = 100
    param_unique_ctr, _ = utils.calibration.read_write_oo('ctr', 'read', method, OUTPUT_DIR=OUTPUT_DIR)
    param_unique_mg, _ = utils.calibration.read_write_oo('mg', 'read', method, OUTPUT_DIR=OUTPUT_DIR)

    # - network inference
    # utils.links.batch_GRN('ctr', method, data_ctr, param, protnames, test_size,time, param_unique_ctr, istart, iend, OUTPUT_DIR)
    utils.links.batch_GRN('mg', method, data_mg, param, protnames, test_size,time, param_unique_mg, istart, iend, OUTPUT_DIR)
