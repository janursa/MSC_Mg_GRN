import os
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import df_target, time_points, OUTPUT_DIR, protnames, GRN_DIR
from utils import process_data
from utils.calibration import retrieve_data
from utils.links import grn, read_write_links, write_scores


if __name__ == '__main__':
    method = 'ridge'
    param = dict(estimator_t=method, alpha=1)

    for study in ['ctr', 'mg', 'combined']:
        data = process_data(df_target, study=study, time_points=time_points(), standardize=False)
        data = np.asarray(data)
        print(f'Data shape: (n_samples_time_series, n_genes) = {data.shape}')
        _, param_unique_ctr = retrieve_data(study=study, method=method, OUTPUT_DIR=OUTPUT_DIR)
        _, train_scores, links_ctr, _, test_scores = grn(
            data=data, gene_names=protnames, time_points=time_points(),
            test_size=0, param=param, param_unique=param_unique_ctr
        )
        read_write_links(links=links_ctr, study=study, mode='write', method=method, output_dir=GRN_DIR)
        output_dir = os.path.join(OUTPUT_DIR, 'GRN', method)
        write_scores(method=method, study=study, trainscores=train_scores, testscores=test_scores, output_dir=output_dir)


    #- compare to vs_string
    # def compare(links):
    #     # top_n = 100
    #     # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
    #     links_short = utils.links.choose_top_quantile(links, quantile=.75)
    #     match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
    #     print('match: ', match_count)
    # compare(links_ctr)


