import sys
import os
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import df_target, time_points, protnames, GRN_DIR, CALIBRATION_DIR
from scripts.utils import process_data
from scripts.utils.calibration import retrieve_data, plot_scores
from scripts.utils.links import grn, read_write_links, write_scores

if __name__ == '__main__':

    #- params
    method='ridge'
    param = dict(estimator_t=method)
    #- network inference
    # - create dir
    if not os.path.isdir(os.path.join(GRN_DIR, method)):
        os.makedirs(os.path.join(GRN_DIR, method))
    test_size = 0
    for study in ['ctr', 'mg', 'combined']:
        if study == 'combined':
            gene_names = protnames()+['mg']
        else:
            gene_names = protnames()
        # - read the data
        data = process_data(df_target(), study=study, time_points=time_points(), standardize=False)
        print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
        #- read results of calibration
        _, param_unique = retrieve_data(study=study, method=method, output_dir=CALIBRATION_DIR)
        _, trainscores, links, _, testscores = grn(data=data, gene_names=gene_names, time_points=time_points(),
                                        test_size=test_size, param=param, param_unique=param_unique)
        #- write to file
        read_write_links(links=links, study=study, mode='write', method=method, output_dir=GRN_DIR)

        output_dir = os.path.join(GRN_DIR, method)
        write_scores(method=method, study=study, trainscores=trainscores, testscores=testscores, output_dir=output_dir)


    #- compare to vs_string
    # def compare(links):
    #     # top_n = 100
    #     # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
    #     links_short = utils.links.choose_top_quantile(links, quantile=.75)
    #     match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
    #     print('match: ', match_count)
    # compare(links_ctr)


