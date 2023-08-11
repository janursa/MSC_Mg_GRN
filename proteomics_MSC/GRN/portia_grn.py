import sys
import os
import numpy as np

from proteomics_MSC.imports import GRN_DIR, CALIBRATION_DIR
from proteomics_MSC.common_tools import read_write_data, calibration, F_DE_data, F_DE_protnames
from proteomics_MSC.common_tools.links import run_portia

import yaml
from importlib.resources import open_text
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    studies = config['studies']
    time_points = config['time_points']

    if not os.path.isdir(os.path.join(GRN_DIR, 'portia')):
        os.makedirs(os.path.join(GRN_DIR, 'portia'))
    for DE_type, DE_proteins in F_DE_protnames().items():
        for study in studies:
            data = read_write_data(mode='read', tag=f'{DE_type}_{study}')
            links_df = run_portia(data, DE_proteins)
            to_save = os.path.join(GRN_DIR, 'portia', f'links_{DE_type}_{study}.csv')
            links_df.to_csv(to_save, index=False)
            print(f'output -> {to_save}')
    print('GRN inference using Portia completed')

