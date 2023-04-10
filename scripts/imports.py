"""

"""
import os
import sys
import pathlib
import json
import pandas as pd
import numpy as np
from pathlib import Path

MAIN_DIR = os.path.join(pathlib.Path(__file__).parent.resolve(), '..')
sys.path.insert(0, MAIN_DIR)

OUTPUT_DIR = os.path.join(MAIN_DIR, 'results')
STATISTICAL_ANALYSIS_DIR =  os.path.join(MAIN_DIR, 'statistical_analysis')
CALIBRATION_DIR = os.path.join(OUTPUT_DIR, 'calibration')
GRN_DIR = os.path.join(OUTPUT_DIR, 'GRN')
DATA_DIR = os.path.join(OUTPUT_DIR, 'data')
VSA_DIR = os.path.join(OUTPUT_DIR, 'VSA')
ENRICH_DIR = os.path.join(OUTPUT_DIR, 'enrichment_analysis')
MODELSELECTION_DIR = os.path.join(OUTPUT_DIR, 'model_selection')
VSA_NOISE_DIR = Path(VSA_DIR)/'noise'
GRN_VISUALIZE_DIR = Path(GRN_DIR)/'visualize'
PLOT_WEIGHTS_DIR = Path(GRN_DIR)/'plot_weights'
RANDOM_MODELS_DIR = Path(MODELSELECTION_DIR) / 'baseline_scores'

from geneRNI.models import get_estimator_wrapper

study_colors = ['lightblue', 'pink'] #color map for ctr and sample, respectively

top_quantiles = np.linspace(.75, .9, 10) #used to cut-off links for epr calculation

def F_selected_models():
    return np.loadtxt(os.path.join(MODELSELECTION_DIR, f'selected_models.txt'), dtype=str, delimiter=",")

def F_model_name_2_method_and_DE_type(model_name):
    model_name_parts = model_name.split('_')
    method = model_name_parts[-1]
    DE_type = '_'.join(model_name_parts[0:2])
    return method, DE_type
def F_DE_data():
    """
    returns DE_data that is a dict of df
    """
    # return pd.read_csv(os.path.join(DATA_DIR,'DE_data.csv'))
    with open(os.path.join(DATA_DIR,'DE_data.csv')) as f:
        data = json.load(f)
    data = {ky: pd.read_json(df_string) for ky, df_string in data.items()}
    return data
def F_DE_protnames():
    """
    returns DE_proteins that is a dict of list. It contains protnames
    """
    with open(os.path.join(DATA_DIR,'DE_protnames.txt')) as f:
        data = eval(f.read())
    return data
def F_protnames_to_genenames():
    with open(Path(DATA_DIR)/ 'protnames_to_genenames.json', 'r') as ff:
        map_protnames_genenames = json.loads(ff.read())
    return map_protnames_genenames
def F_DE_genenames():
    """
        returns DE_proteins that is a dict of list. It contains genenames
    """
    map_protnames_genenames = F_protnames_to_genenames()
    DE_genenames = {}
    for DE_name, protnames in F_DE_protnames().items():
        genenames = [map_protnames_genenames[protname] for protname in protnames]
        DE_genenames[DE_name] = genenames
    return DE_genenames

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

