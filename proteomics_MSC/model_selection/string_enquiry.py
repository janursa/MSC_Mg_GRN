import sys
import os
import argparse
import requests
from typing import List
from urllib.parse import urljoin

import numpy as np
import pandas as pd

from proteomics_MSC.imports import STRING_DIR
from proteomics_MSC.common_tools import F_DE_protnames
from proteomics_MSC.common_tools.model_selection import string_functions_analysis

if __name__ == '__main__':

    if not os.path.isdir(STRING_DIR):
        os.makedirs(STRING_DIR)

    #- functional enrichment
    for DE_type, DE_proteins in F_DE_protnames().items():
        df = string_functions_analysis(DE_proteins)
        to_save = os.path.join(STRING_DIR, f'links_{DE_type}.csv')
        df.to_csv(to_save, index=False)
        print(f'output -> {to_save}')

    #- network enrichment
    if False:
        for DE_type, DE_proteins in F_DE_protnames().items():
            df_enrich = string_network_analysis(DE_proteins)
            df_enrich.to_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index=False)
        print('network analysis is completed')

