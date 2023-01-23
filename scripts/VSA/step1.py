"""
    Conduct vester's sensitivity analysis for non noised links
"""
import sys
import os
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, protnames, VSA_DIR
from scripts.utils.links import read_write_links
from scripts.utils.VSA import VestersSA, VSA_plot, role_change

def retreive_grn(method):
    links_ctr= read_write_links(method=method, study='ctr', mode='read', output_dir=GRN_DIR)
    links_mg= read_write_links(method=method, study='mg', mode='read', output_dir=GRN_DIR)
    return links_ctr, links_mg
if __name__ == '__main__':
    method = 'portia'
    links_ctr, links_sample = retreive_grn(method)
    oo_vsa_ctr = VestersSA(links_ctr, protnames)
    oo_vsa_sample = VestersSA(links_sample, protnames)
    # - plot 1: seperate windows for ctr and sample
    VSA_plot.plot_ctr_vs_sample(oo_vsa_ctr, oo_vsa_sample, protnames, save_dir=VSA_DIR)

    # - plot 2: change in the role from critical to anything or the opposite:
    df_role_change = role_change(oo_vsa_ctr, oo_vsa_sample, target_role=3)
    df_role_change.to_csv(os.path.join(VSA_DIR, f'role_change_{method}.csv'), index=False)
    VSA_plot.plot_role_change(df_role_change, save_dir=VSA_DIR)

