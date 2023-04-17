import sys
import os
import numpy as np
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_DE_data, time_points, GRN_DIR, CALIBRATION_DIR
from common_tools import read_write_data
from common_tools.calibration import retrieve_data, plot_scores
from common_tools.links import run_generni

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    args, remaining_args = parser.parse_known_args()

    studies = args.studies

    verbose = False
    #- params
    method='ridge'
    param = dict(estimator_t=method)
    #- network inference
    # - create dir
    if not os.path.isdir(os.path.join(GRN_DIR, method)):
        os.makedirs(os.path.join(GRN_DIR, method))

    for DE_type, DE_data in F_DE_data().items():
        for study in studies:
            if study == 'all-in':
                gene_names = DE_data['Protein'].values.tolist() + ['mg']
            else:
                gene_names = DE_data['Protein'].values
            # - read the data
            data = data = read_write_data(mode='read', tag=f'{DE_type}_{study}')
            n_timepoints = data.shape[0]
            days = time_points()[0:n_timepoints]
            if verbose:
                print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
            #- read results of calibration
            _, param_unique = retrieve_data(study=study, DE_type=DE_type, method=method, output_dir=CALIBRATION_DIR)
            _, trainscores, links, _, testscores = run_generni(data=data, gene_names=gene_names, time_points=days,
                                             param=param, param_unique=param_unique)
            #- write to file
            links.to_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index=False)
            # np.savetxt(os.path.join(GRN_DIR, method, f'testscores_{DE_type}_{study}.csv'), testscores)
            np.savetxt(os.path.join(GRN_DIR, method, f'trainscores_{DE_type}_{study}.csv'), trainscores)


