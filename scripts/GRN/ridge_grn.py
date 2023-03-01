import sys
import os
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import F_DE_data, time_points, GRN_DIR, CALIBRATION_DIR
from scripts.utils import process_data
from scripts.utils.calibration import retrieve_data, plot_scores
from scripts.utils.links import grn

if __name__ == '__main__':
    verbose = False
    #- params
    method='ridge'
    param = dict(estimator_t=method)
    #- network inference
    # - create dir
    if not os.path.isdir(os.path.join(GRN_DIR, method)):
        os.makedirs(os.path.join(GRN_DIR, method))

    for DE_data_type, DE_data in F_DE_data().items():
        for study in ['ctr', 'mg']:
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values
            # - read the data
            data = process_data(DE_data, study=study, standardize=False)
            n_timepoints = data.shape[0]
            days = time_points()[0:n_timepoints]
            if verbose:
                print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
            #- read results of calibration
            _, param_unique = retrieve_data(study=study, DE_type=DE_data_type, method=method, output_dir=CALIBRATION_DIR)
            _, trainscores, links, _, testscores = grn(data=data, gene_names=gene_names, time_points=days,
                                             param=param, param_unique=param_unique)
            #- write to file
            links.to_csv(os.path.join(GRN_DIR, method, f'links_{DE_data_type}_{study}.csv'), index=False)
            # TODO: test scores are not calculated
            # np.savetxt(os.path.join(GRN_DIR, method, f'testscores_{DE_data_type}_{study}.csv'), testscores)
            np.savetxt(os.path.join(GRN_DIR, method, f'trainscores_{DE_data_type}_{study}.csv'), trainscores)


    #- compare to model_selection
    # def compare(links):
    #     # top_n = 100
    #     # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
    #     links_short = utils.links.choose_top_quantile(links, quantile=.75)
    #     match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
    #     print('match: ', match_count)
    # compare(links_ctr)


