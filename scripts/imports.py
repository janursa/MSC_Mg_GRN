import os
import sys
import importlib
import pathlib


MAIN_DIR = os.path.join(pathlib.Path(__file__).parent.resolve(), '..')
sys.path.insert(0, MAIN_DIR)

import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import copy

from sklearn import preprocessing
from scripts.utils import time
import matplotlib.pyplot as plt
import json
import random
#- define main directories

OUTPUT_DIR = os.path.join(MAIN_DIR, 'results')
CALIBRATION_DIR = os.path.join(OUTPUT_DIR, 'calibration')
#- add paths
sys.path.insert(0,MAIN_DIR)
geneRNI_dir = os.path.join(MAIN_DIR,'..','geneRNI')
sys.path.insert(0, geneRNI_dir)
#- local imports
from scripts import utils
from geneRNI import geneRNI, tools, search_param
from geneRNI.data import Data
#- import training data
df_target = pd.read_csv(os.path.join(MAIN_DIR,'results/postprocess/DE_data.csv'))
#- DE protnames
protnames = list(df_target['Protein'].values)

#- funcs
def retreive_grn(study, method):
    if method == 'portia' or method=='ridge':
        links = utils.Links.read_write_links(study=study, mode='read',method=method,OUTPUT_DIR=OUTPUT_DIR)
    else:
        links = pd.read_pickle(os.path.join(OUTPUT_DIR,'GRN', f'links_{study}_{method}.csv'))
    return links

param_grid = dict(decay_coeff=np.arange(0,1,.05), min_samples_leaf=np.arange(1,5,1), max_depth=np.arange(5,34,1))
