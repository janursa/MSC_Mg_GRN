import sys
import os
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import df_target, time_points, OUTPUT_DIR, protnames, GRN_DIR
from utils import process_data
from utils.calibration import retreive_data, plot_scores
from utils.links import grn, read_write_links, write_scores
if __name__ == '__main__':
    #- read the data
    data_ctr = process_data(df_target, study='ctr', time_points=time_points(), standardize=False)
    data_mg = process_data(df_target, study='mg', time_points=time_points(), standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #- params
    method='ridge'
    param = dict(estimator_t=method, alpha=1)

    #- network inference
    test_size = 0
    _, param_unique_ctr = retreive_data(study='ctr', method=method, OUTPUT_DIR=OUTPUT_DIR)
    _, param_unique_mg = retreive_data(study='mg', method=method, OUTPUT_DIR=OUTPUT_DIR)
    # param_unique_ctr = None
    # param_unique_mg = None

    _, trainscores_ctr, links_ctr, _, testscores_ctr = grn(data=data_ctr, gene_names=protnames, time_points=time_points(),
                                    test_size=test_size, param=param, param_unique=param_unique_ctr)
    _, trainscores_mg, links_mg, _, testscores_mg = grn(data=data_mg, gene_names=protnames, time_points=time_points(),
                                            test_size=test_size, param=param, param_unique=param_unique_mg)
    # print(np.mean(links_ctr['Weight']))
    #- write to files
    read_write_links(links=links_ctr, study='ctr', mode='write', method=method, output_dir=GRN_DIR)
    read_write_links(links=links_mg, study='mg', mode='write', method=method, output_dir=GRN_DIR)

    output_dir = os.path.join(OUTPUT_DIR, 'GRN', method)
    write_scores(method=method, study='ctr', trainscores=trainscores_ctr, testscores=testscores_mg, output_dir=output_dir)
    write_scores(method=method, study='mg', trainscores=trainscores_mg, testscores=testscores_mg, output_dir=output_dir)


    #- compare to vs_string
    # def compare(links):
    #     # top_n = 100
    #     # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
    #     links_short = utils.links.choose_top_quantile(links, quantile=.75)
    #     match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
    #     print('match: ', match_count)
    # compare(links_ctr)


