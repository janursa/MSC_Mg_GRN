"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os
import matplotlib.pyplot as plt

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import GRN_DIR, F_DE_data
from scripts.utils import calibration, serif_font, make_title_pretty
from scripts.utils.links import retrieve_scores

def visualize_scores(method, ylabel, ylim):
    #-- read the data
    DE_types = F_DE_data().keys()
    studies = ['ctr', 'mg', 'all-in']
    trainscores_stack_stack = []
    testscores_stack_stack = []
    for DE_type in DE_types:
        trainscores_stack = []
        testscores_stack = []
        for study in studies:
            trainscores, testscores = retrieve_scores(method=method, study=study, DE_type=DE_type, output_dir=GRN_DIR)
            trainscores_stack.append(trainscores)
            testscores_stack.append(testscores)
        trainscores_stack_stack.append(trainscores_stack)
        testscores_stack_stack.append(testscores_stack)
    #---- plot test scores and train scores seperately
    serif_font()
    ncols = 2
    nrows = int(len(DE_types) / ncols)
    def plot_seperately(scores_stack, ylim):
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(2.5 * ncols, 2.5 * nrows))
        for idx, DE_type in enumerate(DE_types):
            scores = scores_stack[idx]
            i = int(idx / ncols)
            j = idx % ncols
            ax = axes[i][j]
            ax.set_title(make_title_pretty(DE_type))
            ax.set_ylim(ylim)
            calibration.plot_scores(scores, xtickslabels=studies, ylabel=ylabel, ax=ax)
        return fig

    fig = plot_seperately(testscores_stack_stack, ylim=ylim)
    fig.savefig(os.path.join(GRN_DIR, method, 'testscores.png'), dpi=300, transparent=True, facecolor='white')
    fig.savefig(os.path.join(GRN_DIR, method, 'testscores.pdf'), facecolor='white')

    fig = plot_seperately(trainscores_stack_stack, ylim=[0,1.1])
    fig.savefig(os.path.join(GRN_DIR, method, 'trainscores.png'), dpi=300, transparent=True, facecolor='white')
    fig.savefig(os.path.join(GRN_DIR, method, 'trainscores.pdf'), facecolor='white')



if __name__ == '__main__':
    visualize_scores(method='RF', ylabel='Scores', ylim=[-1,1.1])
    # main(method='ridge')
