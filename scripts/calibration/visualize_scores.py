"""
Reads the results of hyper param tuning (RF), 
pool them (because there are multiple runs), average them, and output them to files.
"""
import sys
import os
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import CALIBRATION_DIR, F_DE_data
from scripts.utils import calibration, serif_font, make_title_pretty

def visualize_scores(method, ylabel, ylim):
    # ----------  individual plots for each DE_type---------------
    if False:
        calibration.plot_oo(method, CALIBRATION_DIR, ylabel=ylabel)

    # ----------  combined plot-------------------
    if True:
        DE_types = F_DE_data().keys()
        studies = ['ctr', 'mg', 'all-in']

        serif_font()
        ncols = 2
        nrows = int(len(DE_types) / ncols)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.5 * ncols, 2.5 * nrows))
        for idx, DE_type in enumerate(DE_types):
            best_scores_stack = []
            for study in studies:
                best_scores, best_params = calibration.retrieve_data(
                    study, method=method, DE_type=DE_type, output_dir=CALIBRATION_DIR)
                best_scores_stack.append(best_scores)
            # - plot
            i = int(idx / ncols)
            j = idx % ncols
            ax = axes[i][j]
            ax.set_title(make_title_pretty(DE_type))
            ax.set_ylim(ylim)
            calibration.plot_scores(best_scores_stack, xtickslabels=studies, ylabel=ylabel, ax=ax)

        fig.savefig(os.path.join(CALIBRATION_DIR, method, 'scores.png'), dpi=300, transparent=True, facecolor='white')
        fig.savefig(os.path.join(CALIBRATION_DIR, method, 'scores.pdf'), facecolor='white')


if __name__ == '__main__':
    visualize_scores(method='ridge', ylabel='Leave_one_out score', ylim=[-4, 1.25])
    visualize_scores(method='RF', ylabel='OOB score', ylim=[-.7, 1.25])




