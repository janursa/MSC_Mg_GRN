import os
from subprocess import run
from pathlib import Path
import numpy as np
from scripts.imports import MAIN_DIR, DATA_DIR, STATISTICAL_ANALYSIS_DIR, ENRICH_DIR, CALIBRATION_DIR, GRN_DIR, \
    MODELSELECTION_DIR, VSA_DIR, VSA_NOISE_DIR, F_selected_models, F_model_name_2_method_and_DE_type, \
    F_DE_data,RANDOM_REGULATORS_DIR, MODELSELECTION_BASELINE_DIR, NETWORK_ANALYSIS_DIR

def selected_models():
    if os.path.exists(Path(MODELSELECTION_DIR)/ 'selected_models.txt'):
        return F_selected_models()
    else:
        return []