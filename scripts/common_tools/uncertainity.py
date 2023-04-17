"""
    This class handles the functions for the uncertainity_analysis analysis
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple, Callable

import pandas as pd
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from common_tools import serif_font
from common_tools.VSA import role_analysis, RolePlot

def add_noise(data_org, noise_type, noise_intensity):

    if noise_type == 'mp':
        noisy_data = _add_mpnoise(data_org, std=noise_intensity)
    else:
        noisy_data = _add_adnoise(data_org, std=noise_intensity)

    return noisy_data
def create_noisy_data(n_repeat, **kwargs) -> List[pd.DataFrame]:
    """Create n number of random noise"""
    noisy_sets = []
    for i in range(n_repeat):
        noisy_data = add_noise(**kwargs)
        noisy_sets.append(noisy_data)
    return noisy_sets
def create_noisy_links(data_org:pd.DataFrame, n_repeat:int, noise_type:str,
                       noise_intensity:float,  gene_names: List[str], grn_function: Callable[[pd.DataFrame, List[str]], pd.DataFrame]) -> List[pd.DataFrame]:
    """
    add noise to the data, run grn, and save it. If file already exist, skip it unless forced.
    Collect the links and retunt a list.
    """

    # - create n noisy dataset
    noisy_datasets = create_noisy_data(data_org=data_org, noise_type=noise_type, noise_intensity=noise_intensity, n_repeat=n_repeat)
    # - run grn for each noisy link
    noisy_links = []
    with tqdm(total=n_repeat, desc=f'Run GRN for noisy data: {noise_type} {noise_intensity}') as pbar:
        for data_i, data in enumerate(noisy_datasets):
            # - run GRN
            links = grn_function(data, gene_names)
            noisy_links.append(links)
            pbar.update(1)
    return noisy_links



def _add_mpnoise(data: np.ndarray, std=.3) -> np.array:
    """ Multiplicative noise
         Creates n_noisy_datasets data
    """
    data_std = np.std(data)
    applied_std = std * np.array(data_std)

    frac_noise = np.random.rand()
    noise_mask = np.random.choice([0, 1], size=data.shape, p=[1 - frac_noise, frac_noise])
    rand_values = np.random.normal(loc=1, scale=applied_std, size=data.shape)
    noisy_data = data * (1 + (rand_values - 1) * noise_mask)


    return noisy_data

def _add_adnoise(data: np.ndarray, std=.05) -> np.array:
    """ Additive noise
             Creates n_relica noised links
    """
    data_std = np.std(data)
    applied_std = data_std * std
    rand_values = np.random.normal(loc=1, scale=applied_std, size=data.shape)
    noisy_data = rand_values + data

    return noisy_data

