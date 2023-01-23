import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import time_points, protnames, OUTPUT_DIR, MAIN_DIR, df_target, GRN_DIR
from utils import process_data
from utils.links import read_write_links

from geneRNI.links import format_links
dir_portia = os.path.join(MAIN_DIR,'..','external/PORTIA-master')
sys.path.insert(0,dir_portia)
import portia as pt
def GRN(data, study):
    dataset = data_process(data, time_points())
    M_bar = pt.run(dataset, method='fast')
    links_df = format_links(M_bar, protnames)
    read_write_links(links=links_df, study=study, mode='write',method='portia', output_dir=GRN_DIR)
    return links_df

def data_process(data, time_points):
#     dataset = Data(gene_names=protnames, ss_data=None, ts_data=[data], time_points=[time],
#                        test_size=0)
    portia_dataset = pt.GeneExpressionDataset()
    exp_id = 1
    for data_i, time_i in zip(data, time_points):
#         X, y = dataset.process_time_series(0)
        portia_dataset.add(pt.Experiment(exp_id, data_i))
        exp_id+=1
    return portia_dataset

if __name__ == '__main__':
    #- read the data
    data_ctr = process_data(df_target, study='ctr', time_points=time_points(), standardize=False)
    data_mg = process_data(df_target, study='mg', time_points=time_points(), standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    # - run
    # decay_coeffs = utils.estimate_decay_rates([data_ctr,data_mg], [time, time])
    links_ctr = GRN(data_ctr, 'ctr')
    links_mg = GRN(data_mg, 'mg')
    #- compare to vs_string
    # def compare(links):
    #     top_n = 150
    #     # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
    #     links_short = choose_top_quantile(links, quantile=.9)
    #     match_count = compare_network_string(links_short.copy(), OUTPUT_DIR)
    #     print('match: ', match_count)
    # compare(links_ctr)