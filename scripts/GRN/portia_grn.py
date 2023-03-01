import sys
import os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import time_points, MAIN_DIR, F_DE_data, GRN_DIR
from scripts.utils import process_data

from geneRNI.links import format_links

dir_portia = os.path.join(MAIN_DIR, '..', 'external/PORTIA-master')
sys.path.insert(0, dir_portia)
import portia as pt

if __name__ == '__main__':
    #- create dir
    if not os.path.isdir(os.path.join(GRN_DIR, 'portia')):
        os.makedirs(os.path.join(GRN_DIR, 'portia'))
    for DE_data_type, DE_data in F_DE_data().items():
        for study in ['ctr', 'mg']:
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values
            data = process_data(DE_data, study=study, standardize=False)
            n_timepoints = data.shape[0]
            days = time_points()[0:n_timepoints]
            data = np.asarray(data)
            # print(f'Data shape: (n_samples_time_series, n_genes) = {data.shape}')
            # - run
            # decay_coeffs = utils.estimate_decay_rates([data_ctr,data_mg], [time, time])
            #- create portia dataset
            portia_dataset = pt.GeneExpressionDataset()
            for exp_id, (data_i, time_i) in enumerate(zip(data, days)):
                portia_dataset.add(pt.Experiment(exp_id, data_i))
            #- GRN inference
            M_bar = pt.run(portia_dataset, method='fast')
            links_df = format_links(M_bar, gene_names)
            links_df.to_csv(os.path.join(GRN_DIR, 'portia', f'links_{DE_data_type}_{study}.csv'), index=False)
            #- compare to model_selection
            # def compare(links):
            #     top_n = 150
            #     # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
            #     links_short = choose_top_quantile(links, quantile=.9)
            #     match_count = compare_network_string(links_short.copy(), OUTPUT_DIR)
            #     print('match: ', match_count)
            # compare(links_ctr)