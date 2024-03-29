import sys
import os
import numpy as np
import yaml
from importlib.resources import open_text
import matplotlib.pyplot as plt

from proteomics_MSC.imports import CALIBRATION_DIR
from proteomics_MSC.common_tools import F_DE_data, make_title_pretty, serif_font, calibration
with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

def get_and_visualize_scores(method, ylabel, yticks):
    serif_font()
    DE_types = F_DE_data().keys()
    fig, ax = plt.subplots(nrows=1, ncols=1, tight_layout=True, figsize=(4 , 2))
    best_scores_stack = []

    for idx, DE_type in enumerate(DE_types):
        for study in studies:
            best_scores, best_params = calibration.retrieve_data(
                study, method=method, DE_type=DE_type, output_dir=CALIBRATION_DIR)
            best_scores_stack.append(best_scores.tolist())
    best_scores_stack = [item[0:10] for item in best_scores_stack]
    # aa
    # - x ticks
    x_ticks_locations = [i + 1 + (i // 4) for i in range(len(best_scores_stack))]
    ax.set_xticks(x_ticks_locations)
    ax.set_xticklabels(np.tile(studies, 4), rotation=0)
    ax.axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=.5)
    ax.set_xmargin(.03)
    # - plot
    bplot = ax.violinplot(np.asarray(best_scores_stack).T, positions=x_ticks_locations, showmeans=True, showextrema=False)
    # - y ticks
    ax.set_ylabel(ylabel)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ymin, ymax = ax.get_ylim()
    margin = 0.4 * (ymax - ymin)
    ymax_new = ymax + margin
    ax.set_ylim(ymin, ymax_new)
    # - face colors
    colors = np.tile(np.repeat(group_colors, 2),2)
    for patch, color in zip(bplot['bodies'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('grey')
        patch.set_alpha(1)
    # - annotate mean scores
    for i, body in enumerate(bplot['bodies']):
        y_max = np.max(body.get_paths()[0].vertices[:, 1])
        mean = np.mean(best_scores_stack[i])
        xy = (x_ticks_locations[i], y_max)
        mean_value =  np.round(mean, 1) if np.round(mean, 1) != -0.0 else 0.0
        ax.annotate(mean_value, xy=xy,
                    ha='center',
                    va='bottom',
                    verticalalignment='baseline',
                    textcoords='offset points',
                    fontsize=9,
                    xytext=(0, 3)
                    )

    fig.savefig(os.path.join(CALIBRATION_DIR, method, 'scores.png'), dpi=300, transparent=True, facecolor='white')
    fig.savefig(os.path.join(CALIBRATION_DIR, method, 'scores.pdf'), facecolor='white')


if __name__ == '__main__':

    studies =  config['studies']
    group_colors = ["#1E90FF", "#FFA07A"]  # for imputation types

    get_and_visualize_scores(method='ridge', ylabel=r'$\mathrm{R}^2$ score (LOO)', yticks=[-4, -2, 0, 1])
    get_and_visualize_scores(method='RF', ylabel=r'$\mathrm{R}^2$ score (OOB)', yticks=[-1, 0, 1])

    # - plot the legends for the groups
    fig, ax = plt.subplots(nrows=1, ncols=1, tight_layout=True, figsize=(4, 2))
    handles = []
    labels = ['MinProb', 'KNN']
    for i, color in enumerate(group_colors):
        handles.append(ax.scatter([], [], marker='o', label=labels[i], color=color,
                                  s=30, alpha=1))
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1, handles=handles, title='', fancybox=False,
              frameon=False, prop={'size': 10}, title_fontproperties={'size': 9, 'weight': 'bold'})
    to_save = os.path.join(CALIBRATION_DIR, 'RF', 'scores_legend.png')
    fig.savefig(to_save, dpi=300, transparent=True, facecolor='white')
    print(f'output -> {to_save}')




