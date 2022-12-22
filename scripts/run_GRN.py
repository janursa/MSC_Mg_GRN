import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import importlib
import sys
import copy
import os
from sklearn import preprocessing
import matplotlib.pyplot as plt
import json
import random
import importlib
#- define main directories
MAIN_DIR = 'C:/Users/nourisa/Downloads/testProjs/proteomics_MSC/'
OUTPUT_DIR = os.path.join(MAIN_DIR, 'results')
CALIBRATION_DIR = os.path.join(OUTPUT_DIR, 'calibration')
#- add paths
sys.path.insert(0,MAIN_DIR)
geneRNI_dir = os.path.join(MAIN_DIR,'..','geneRNI')
sys.path.insert(0, geneRNI_dir)
#- local imports
from scripts.utils import time
from scripts import utils
from geneRNI import geneRNI, tools, search_param
#- import training data
df_target = pd.read_csv('results/postprocess/DE_data.csv')
#- DE protnames
protnames = list(df_target['Protein'].values)
#- calibration 
param_grid = dict(alpha=np.arange(0,1,.05), min_samples_split=np.arange(2,5,1), min_samples_leaf=np.arange(1,5,1),max_depth=np.arange(10,20,1))
def GRN(data, study='ctr', i=0):
#     print(f'----GRN for {study}-----')
    Xs, ys = tools.Data.process_time_series(TS_data=[data], time_points=[time],  gene_names=protnames)
    # read the outputs of tunning
    best_params, _ = utils.read_write_oo(study, mode='read', OUTPUT_DIR=OUTPUT_DIR)
    # run the network inference
    param = dict(estimator_t='RF')

    ests, train_scores, links_df, oob_scores, test_scores = \
        geneRNI.network_inference(Xs, ys, gene_names=protnames, 
                                  param=param, param_unique=best_params, 
                                  Xs_test=None, ys_test=None, verbose=False
                                  )
    utils.read_write_links(links=links_df, study=study, mode='write',i=i,OUTPUT_DIR=OUTPUT_DIR)
def run_batch(data, study,istart, iend):
    for i in range(istart, iend):
        GRN(data, study, i)
        print(f'{i} finished')
data_ctr = np.genfromtxt(os.path.join(OUTPUT_DIR,'data', 'data_ctr.csv'),  delimiter=',')
data_mg = np.genfromtxt(os.path.join(OUTPUT_DIR,'data', 'data_mg.csv'),  delimiter=',')

# run_batch(data_ctr,'ctr',961, 1000) 
run_batch(data_ctr,'mg',869, 1000) 