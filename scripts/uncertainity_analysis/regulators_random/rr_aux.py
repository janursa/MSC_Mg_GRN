import numpy as np
import os
from pathlib import Path
from imports import RANDOM_REGULATORS_DIR
from typing import List
def retrieve_random_scores(selected_model:str) -> List[float]:
    with open(Path(RANDOM_REGULATORS_DIR) / f'ep_scores_random_{selected_model}.csv', 'r') as ff:
        ep_scores_random = np.loadtxt(ff, delimiter=',')
    return ep_scores_random
def save_random_scores(selected_model:str, ep_score_stack:List[float]) -> None:
    # - create directories
    if not os.path.isdir(RANDOM_REGULATORS_DIR):
        os.makedirs(RANDOM_REGULATORS_DIR)
    with open(Path(RANDOM_REGULATORS_DIR) / f'ep_scores_random_{selected_model}.csv', 'w') as ff:
        np.savetxt(ff, ep_score_stack, delimiter=',')