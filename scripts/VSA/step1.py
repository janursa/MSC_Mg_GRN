"""
    Conduct vester's sensitivity analysis for non noised links
"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, VSA_DIR, F_DE_protiens
from scripts.utils.VSA import VestersSA, VSA_plot, role_change
from scripts.utils import make_title_pretty, comic_font

def retreive_grn(method, DE_type):
    links_stack = []
    for study in ['ctr','mg']:
        links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
        links_stack.append(links)
    return links_stack
if __name__ == '__main__':
    comic_font()
    #- create dir
    if not os.path.isdir(VSA_DIR):
        os.makedirs(VSA_DIR)

    methods = ['RF', 'portia', 'ridge']
    # methods = ['RF'] #TODO: remove this
    combined_plot = True  # to plot all DE_types in one fig. only for role change

    for method in methods:
        if combined_plot:
            ncols = 2
            nrows = 3
            fig_rl, axes_rl = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(3 * ncols, 3 * nrows))

        if not os.path.isdir(os.path.join(VSA_DIR,method)):
            os.makedirs(os.path.join(VSA_DIR,method))
        for idx, (DE_type, DE_proteins) in enumerate(F_DE_protiens().items()):
            links_ctr, links_sample = retreive_grn(method, DE_type)

            oo_vsa_ctr = VestersSA(links_ctr, DE_proteins)
            oo_vsa_sample = VestersSA(links_sample, DE_proteins)
            # -------- plot 1: ctr vs sample --------------
            if True:
                fig = VSA_plot.plot_ctr_vs_sample(oo_vsa_ctr, oo_vsa_sample, DE_proteins)
                # fig.savefig(os.path.join(VSA_DIR, method, f'roles_{DE_type}.pdf', bbox_extra_artists=(ll,), bbox_inches='tight')
                fig.savefig(os.path.join(VSA_DIR, method, f'roles_{DE_type}.pdf'),  bbox_inches='tight')
                fig.savefig(os.path.join(VSA_DIR, method, f'roles_{DE_type}.png'), dpi=300, transparent=True, bbox_inches='tight')

            # -------- plot 2: from ctr to sample -----
            df_role_change = role_change(oo_vsa_ctr, oo_vsa_sample, target_role=3)
            if True: # to write role change to files
                df_role_change.to_csv(os.path.join(VSA_DIR, method, f'role_change_{DE_type}.csv'), index=False)

            if combined_plot:
                i = int(idx/ncols)
                j = idx%ncols
                ax = axes_rl[i][j]
                VSA_plot.plot_role_change(df_role_change, ax=ax, title=make_title_pretty(DE_type))
            else:
                fig = VSA_plot.plot_role_change(df_role_change)
                fig.savefig(os.path.join(VSA_DIR, method, f'role_change_{DE_type}.png'), dpi=300, transparent=True)
                fig.savefig(os.path.join(VSA_DIR, method, f'role_change_{DE_type}.pdf'))

        if combined_plot:
            fig_rl.savefig(os.path.join(VSA_DIR, method, 'role_changes.png'), dpi=300, transparent=True)
            fig_rl.savefig(os.path.join(VSA_DIR, method, 'role_changes.pdf'))