import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *
if __name__ == '__main__':
    #- read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=True)
    data_mg = utils.process_data(df_target, study='mg', standardize=True)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #- params
    param = dict(estimator_t='ridge')
    # estimated_decay_coeffs = utils.estimate_decay_rates([data_ctr, data_mg], [time, time])

    #- network inference
    test_size = 0
    param_unique_ctr,_ = utils.calibration.read_write_oo('ctr', 'read', 'ridge', OUTPUT_DIR=OUTPUT_DIR)
    param_unique_mg,_ = utils.calibration.read_write_oo('mg', 'read', 'ridge', OUTPUT_DIR=OUTPUT_DIR)

    _,_,links_ctr = utils.links.network_inference(data_ctr, protnames, time, test_size, param, param_unique_ctr)
    _,_,links_mg = utils.links.network_inference(data_mg, protnames, time, test_size, param, param_unique_mg)

    utils.links.read_write_links(links=links_ctr, study='ctr', mode='write_links', method='ridge', OUTPUT_DIR=OUTPUT_DIR)
    utils.links.read_write_links(links=links_mg, study='mg', mode='write_links', method='ridge', OUTPUT_DIR=OUTPUT_DIR)

    #- compare to string
    def compare(links):
        top_n = 100
        links_short = utils.links.choose_top_count(links_ctr, n=top_n)
        # links_short = utils.links.choose_top_quantile(links, quantile=.75)
        match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
        print('match: ', match_count)
    compare(links_ctr)


