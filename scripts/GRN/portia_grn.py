import sys
import os
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import F_DE_protnames, GRN_DIR
from common_tools import read_write_data
from common_tools.links import run_portia


if __name__ == '__main__':
    #- arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    args, remaining_args = parser.parse_known_args()

    studies = args.studies

    if not os.path.isdir(os.path.join(GRN_DIR, 'portia')):
        os.makedirs(os.path.join(GRN_DIR, 'portia'))
    for DE_type, DE_proteins in F_DE_protnames().items():
        for study in studies:
            data = read_write_data(mode='read', tag=f'{DE_type}_{study}')
            links_df = run_portia(data, DE_proteins)
            links_df.to_csv(os.path.join(GRN_DIR, 'portia', f'links_{DE_type}_{study}.csv'), index=False)
    print('GRN inference using Portia completed')

