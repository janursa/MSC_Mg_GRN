import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *
dir_portia = os.path.join(MAIN_DIR,'..','external/PORTIA-master')
sys.path.insert(0,dir_portia)
import portia as pt
def GRN(data, study, decay_coeffs):
    dataset = data_process(data, time, decay_coeffs)
    M_bar = pt.run(dataset, method='fast')
    links_df = tools.Links.format(M_bar, protnames)
    utils.links.read_write_links(links=links_df, study=study, mode='write',method='portia', output_dir=OUTPUT_DIR)
    return links_df

def data_process(data, time, decay_coeffs):
#     dataset = Data(gene_names=protnames, ss_data=None, ts_data=[data], time_points=[time],
#                        test_size=0)
    portia_dataset = pt.GeneExpressionDataset()
    exp_id = 1
    for data_i, time_i in zip(data,time):
#         X, y = dataset.process_time_series(0)
        portia_dataset.add(pt.Experiment(exp_id, data_i))
        exp_id+=1
    return portia_dataset

if __name__ == '__main__':
    # - read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
    data_mg = utils.process_data(df_target, study='mg', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    # - run
    decay_coeffs = utils.estimate_decay_rates([data_ctr,data_mg], [time, time])
    links_ctr = GRN(data_ctr, 'ctr', decay_coeffs)
    links_mg = GRN(data_mg, 'mg', decay_coeffs)
    #- compare to string
    def compare(links):
        top_n = 150
        # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
        links_short = utils.links.choose_top_quantile(links, quantile=.9)
        match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
        print('match: ', match_count)
    compare(links_ctr)