"""

"""
import os
import sys
import pathlib
import json
import pandas as pd
import numpy as np

MAIN_DIR = os.path.join(pathlib.Path(__file__).parent.resolve(), '..')
sys.path.insert(0, MAIN_DIR)

OUTPUT_DIR = os.path.join(MAIN_DIR, 'results')
STATISTICAL_ANALYSIS_DIR =  os.path.join(MAIN_DIR, 'statistical_analysis')
CALIBRATION_DIR = os.path.join(OUTPUT_DIR, 'calibration')
GRN_DIR = os.path.join(OUTPUT_DIR, 'GRN')
DATA_DIR = os.path.join(OUTPUT_DIR, 'data')
VSA_DIR = os.path.join(OUTPUT_DIR, 'VSA')
ENRICH_DIR = os.path.join(OUTPUT_DIR, 'enrichment_analysis')
VS_STRING_DIR = os.path.join(OUTPUT_DIR, 'vs_string')

geneRNI_dir = os.path.join(MAIN_DIR,'..','geneRNI')
sys.path.insert(0, geneRNI_dir)

from geneRNI.models import get_estimator_wrapper

#- import training data
def F_DE_data():
    """
    returns DE_data that is a dict of df
    """
    # return pd.read_csv(os.path.join(DATA_DIR,'DE_data.csv'))
    with open(os.path.join(DATA_DIR,'DE_data.csv')) as f:
        data = json.load(f)
    data = {ky: pd.read_json(df_string) for ky, df_string in data.items()}
    return data
def F_DE_protiens():
    """
    returns DE_proteins that is a dict of list
    """
    with open(os.path.join(DATA_DIR,'DE_protnames.txt')) as f:
        data = eval(f.read())
    return data
#- DE protnames
# def protnames():
#     return list(df_target()['Protein'].values)


def param_grid_RF(n_proteins):
    return dict(
        decay_coeff=np.arange(0,1,.05),
        # min_samples_leaf=np.arange(1,5,1),
        # max_depth=np.arange(5,34,1),
        max_features=np.arange(int(np.sqrt(n_proteins)),n_proteins,1),
    )
def param_grid_ridge():
    param_grid_ridge_ = get_estimator_wrapper('ridge').get_grid_parameters()
    return {**param_grid_ridge_,'decay_coeff':np.arange(0,1,.05)}

def time_points():
    """
    Measurement time points
    """
    return [1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 21]

