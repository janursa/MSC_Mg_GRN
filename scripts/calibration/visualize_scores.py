"""
Reads the results of hyper param tuning (RF), 
pool them (because there are multiple runs), average them, and output them to files.
"""
import sys
import os
import matplotlib.pyplot as plt
import matplotlib
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import CALIBRATION_DIR, F_DE_data, make_title_pretty
from utils import calibration, serif_font

def visualize_scores(method, ylabel, xticks, ylim):
    serif_font()
    # matplotlib.rcParams.update({'font.size': 12})
    DE_types = F_DE_data().keys()
    # studies = ['ctr', 'mg', 'all-in']
    studies = ['ctr', 'mg']

    ncols = 4
    nrows = int(len(DE_types) / ncols)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(1.6 * ncols, 1.8 * nrows))
    for idx, DE_type in enumerate(DE_types):
        best_scores_stack = []
        for study in studies:
            best_scores, best_params = calibration.retrieve_data(
                study, method=method, DE_type=DE_type, output_dir=CALIBRATION_DIR)
            best_scores_stack.append(best_scores)
        # - plot
        ax = axes[idx]
        title = DE_type
        title_parts = title.split('_')
        title = '_'.join([title_parts[0]+'\n',title_parts[1]])
        ax.set_title(make_title_pretty(title))
        bplot = ax.boxplot(best_scores_stack, notch=True, patch_artist=True, meanline=True, widths=.5)
        ax.set_ylabel(ylabel)
        ax.set_yticks(xticks)
        ax.set_yticklabels(xticks)
        # if idx == 0:
        #     ax.set_ylabel(ylabel)
        #     ax.set_yticks(xticks)
        #     ax.set_yticklabels(xticks)
        #
        # else:
        #     ax.set_yticks(xticks)
        #     ax.set_yticks([])
        #     ax.set_ylabel('')
        ax.set_ylim(ylim)
        ax.set_xticks(range(1, len(studies) + 1))
        ax.set_xticklabels(studies, rotation=0)
        ax.axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=.5)
        ax.set_ymargin(.1)
        ax.set_xmargin(.15)
        # - face colors
        colors = ['orange', 'lightgreen', 'lightblue', 'cyan']
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_edgecolor('black')
            patch.set_alpha(1)
    fig.savefig(os.path.join(CALIBRATION_DIR, method, 'scores.png'), dpi=300, transparent=True, facecolor='white')
    fig.savefig(os.path.join(CALIBRATION_DIR, method, 'scores.pdf'), facecolor='white')


if __name__ == '__main__':
    visualize_scores(method='ridge', ylabel='LOO score', xticks=[-2, -1, 0, 1],  ylim= [-2.94, 1.3])
    visualize_scores(method='RF', ylabel='OOB score', xticks=[-1, 0, 1], ylim= [-1.2, 1.2])






