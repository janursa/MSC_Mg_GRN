import sys
import os
import numpy as np
import yaml
from importlib.resources import open_text

from proteomics_MSC.imports import GRN_DIR, CALIBRATION_DIR
from proteomics_MSC.common_tools import read_write_data, calibration, F_DE_data
from proteomics_MSC.common_tools.links import batch_run_generni, run_generni


with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    studies = config['studies']
    time_points = config['time_points']

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
            days = time_points[0:n_timepoints]
            if verbose:
                print('Data shape:', np.array(data).shape, '(n_samples_time_series*n_genes)')
            #- read results of calibration
            _, param_unique = calibration.retrieve_data(study=study, DE_type=DE_type, method=method, output_dir=CALIBRATION_DIR)
            _, trainscores, links, _, testscores = run_generni(data=data, gene_names=gene_names, time_points=days,
                                             param=param, param_unique=param_unique)
            #- write to file
            to_save = os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv')
            links.to_csv(to_save, index=False)
            print(f'output -> {to_save}')
            # np.savetxt(os.path.join(GRN_DIR, method, f'testscores_{DE_type}_{study}.csv'), testscores)
            np.savetxt(os.path.join(GRN_DIR, method, f'trainscores_{DE_type}_{study}.csv'), trainscores)


