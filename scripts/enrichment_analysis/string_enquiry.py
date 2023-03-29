import sys
import os
import argparse
import requests
from typing import List
from urllib.parse import urljoin

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import ENRICH_DIR, F_DE_protnames
from utils.enrich_analysis import string_functions_analysis, string_network_analysis



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--string_api', nargs='+', default="https://version-11-5.string-db.org/api", dest='string_api')
    parser.add_argument('--do_func_analysis', default=True, dest='do_func_analysis')
    parser.add_argument('--do_network_analysis', default=True, dest='do_network_analysis')
    args, remaining_args = parser.parse_known_args()

    string_api = args.string_api
    do_func_analysis = args.do_func_analysis
    do_network_analysis = args.do_network_analysis

    if not os.path.isdir(ENRICH_DIR):
        os.makedirs(ENRICH_DIR)

    #- functional enrichment
    if do_func_analysis:
        for DE_type, DE_proteins in F_DE_protnames().items():
            df_enrich = string_functions_analysis(DE_proteins, string_api)
            df_enrich.to_csv(os.path.join(ENRICH_DIR, f'enrichment_all_{DE_type}.csv'), index=False)
        print('functional analysis is completed')

    #- network enrichment
    if do_network_analysis:
        for DE_type, DE_proteins in F_DE_protnames().items():
            df_enrich = string_network_analysis(DE_proteins, string_api)
            df_enrich.to_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index=False)
        print('network analysis is completed')

