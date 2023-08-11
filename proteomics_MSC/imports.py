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

OUTPUT_DIR = f'{MAIN_DIR}/results'
DATA_DIR = f'{OUTPUT_DIR}/data'
STATISTICAL_ANALYSIS_DIR =  f'{MAIN_DIR}/statistical_analysis'

CALIBRATION_DIR = f'{OUTPUT_DIR}/calibration'
GRN_DIR = f'{OUTPUT_DIR}/GRN'

MODELSELECTION_DIR = f'{OUTPUT_DIR}/model_selection'
MODELSELECTION_BASELINE_DIR = f'{MODELSELECTION_DIR}/baseline_models'
STRING_DIR = f'{MODELSELECTION_DIR}/string'
RANDOM_REGULATORS_DIR = f'{MODELSELECTION_DIR}/random_regulators'

VSA_DIR = f'{OUTPUT_DIR}/VSA'
VSA_VERBOSE_DIR = f'{VSA_DIR}/verbose'
VSA_ROBUSTNESS_DIR = f'{VSA_DIR}/robustness'

DIVERGENCE_ANALYSIS_DIR = f'{OUTPUT_DIR}/divergence_analysis'
DIVERGENCE_ANALYSIS_ROBUSTNESS_DIR = f'{DIVERGENCE_ANALYSIS_DIR}/robustness'

TARGETS_DIR = f'{OUTPUT_DIR}/target_genes'

original_omics_file = f'{MAIN_DIR}/data/original_omics.xlsx'
processed_omics_file = f'{OUTPUT_DIR}/data/edited_data.csv'

DE_protnames_file = Path(DATA_DIR)/'DE_protnames.txt'
DE_data_file = Path(DATA_DIR)/'DE_data.csv'
DE_common_proteins_file = Path(DATA_DIR)/'common_proteins.csv'

