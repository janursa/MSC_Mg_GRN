"""
    Sets of functions useful for hyperparameter tuning
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
import utils

from ._imports import *
def main(study, method, data, gene_names, time_points, param_grid, OUTPUT_DIR, i_start, i_end, test_size, **specs):
    """
        Interface to search function of geneRNI
    """
    print(f'----Tuning for method: {method}, study: {study}')
    dataset = Data(gene_names=gene_names, ss_data=None, ts_data=[data], time_points=[time_points])
    DIR = os.path.join(OUTPUT_DIR, 'calibration', method, study)
    search_param.rand_search(dataset,
                             param_grid=param_grid,
                             output_dir=DIR,
                             i_start=i_start,
                             i_end=i_end,
                             **specs)
    best_scores, best_params = search_param.pool(DIR, n_repeat=i_end)
    #- save
    np.save(os.path.join(OUTPUT_DIR, 'calibration', method, f'best_scores_{study}.npy'), best_scores)
    np.save(os.path.join(OUTPUT_DIR, 'calibration', method, f'best_params_{study}.npy'), best_params)

def plot_scores_pool(bestscores_pool_ctr, bestscores_pool_sample, xticks_labels):
    """plots scores as a box plot for a set"""
    fig, axes = plt.subplots(2, 1, tight_layout=True, figsize=(10, 6))
    data_s = [bestscores_pool_ctr, bestscores_pool_sample]
    titles = ['ctr', 'mg']
    for i, data in enumerate(data_s):
        if data is None:
            continue
        ax = axes[i]
        # ax.boxplot(np.array(data).T)
        bplot = ax.violinplot(data, showmeans=True, showextrema=False)
        ax.set_ylabel('Testing score')
        ax.set_xlabel('')
        ax.set_title(titles[i])
        ax.set_xticks(range(1, len(xticks_labels) + 1))
        ax.set_xticklabels(xticks_labels, rotation=90)
        ax.set_ymargin(.1)
        ax.set_xmargin(.02)
        # colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan','grey']
        for patch in bplot['bodies']:
            patch.set_facecolor('lightgreen')
            patch.set_edgecolor('black')
            patch.set_alpha(1)



    return fig
def plot_bestparams_pool(data, priors, xticks_labels):
    """
        Plots boxplot for indivual param in seperate window. In each window, the variation in
        best params are given for each gene seperately, obtained during different runs of tunning.
    """
    nrows = len(priors)
    ncols = 1
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(10 * ncols, 9))

    for i, key in enumerate(priors.keys()):

        ax = axes[i]
        pool = [item[key] for item in data]
        bplot = ax.violinplot(pool, showmeans=True, showextrema=False)
        ax.set_title(key)
        ax.set_ylabel('Quantity')
        ax.set_xticks(range(1, len(xticks_labels) + 1))
        ax.set_xticklabels(xticks_labels, rotation=90)
        ax.set_ymargin(.1)
        ax.set_xmargin(.02)
        # colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan','grey']
        for patch in bplot['bodies']:
            patch.set_facecolor('lightgreen')
            patch.set_edgecolor('black')
            patch.set_alpha(1)
    return fig
def plot_scores(data_ctr, data_sample, xlabel='', ylabel='OOB score'):
    """plots oob scores as a box plot for ctr and mg side by side"""
    utils.serif_font()
    fig, axes = plt.subplots(1, 1, tight_layout=True, figsize=(2.5, 3),
                             # gridspec_kw={'width_ratios': [2, 2]}
                             )
    data_s = [data_ctr, data_sample]
    labels = ['ctr', 'mg']
    ax = axes
    bplot = ax.boxplot(data_s, notch=True, widths=[.5, .5], patch_artist=True, meanline=True)
    # bplot = ax.violinplot(data_s, showmeans=True, showextrema=True, bootstrap=True
    #     )
    ax.set_ylabel(ylabel)
    # ax.set_title('(A)')
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=0)
    ax.axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=.5)
    ax.set_ymargin(.1)
    ax.set_xmargin(.15)
    # - face colors
    colors = ['pink', 'lightblue', 'lightgreen', 'cyan']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(1)

    # colors = ['black' for i in range(len(labels))]
    # tags = ['cbars']
    # for tag in tags:
    #     bplot[tag].set_color(colors)
    #     bplot[tag].set_decay_coeff(.5)
    return fig
def plot_bestparams(data_ctr, data_sample, priors):
    """
        Plots boxplot for indivual param in seperate window.
    """
    utils.serif_font()
    fig, axes = plt.subplots(1, 2, tight_layout=True, figsize=(7, 4))
    datas = [data_ctr, data_sample]
    titles = ['ctr', 'mg']

    def normalize(xx, priors):
        xx = {key: np.array(list(set(values))) for key, values in xx.items()}
        return {key: (values - min(priors[key])) / (max(priors[key]) - min(priors[key])) for key, values in
                xx.items()}

    for i, data in enumerate(datas):
        data = {key: [item[key] for item in data] for key in data[0].keys()}
        data = normalize(data, priors)
        bplot = axes[i].boxplot(data.values(), notch=False, widths=.4, patch_artist=True, meanline=True)

        axes[i].set_title(titles[i])
        axes[i].set_ylabel('Normalized quantity')
        axes[i].set_xticklabels(data.keys(), rotation=45)
        axes[i].set_ymargin(.15)
        axes[i].set_xmargin(.15)
        colors = ['c', 'olive', 'lightgreen', 'g']
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
    return fig
def retreive_data(study, method, OUTPUT_DIR):
    """
        Reads results of pooling
    """
    best_scores = np.load(os.path.join(OUTPUT_DIR, 'calibration', method, f'best_scores_{study}.npy'))
    best_params = np.load(os.path.join(OUTPUT_DIR, 'calibration', method, f'best_params_{study}.npy'),allow_pickle=True)
    return best_scores, best_params
def plot_oo(method, priors, protnames, OUTPUT_DIR):
    """
     Plots a series of graphs for best params and best scores (individual protein and combined)
    """
    dir = os.path.join(OUTPUT_DIR, 'calibration', method)
    best_scores_ctr, best_params_ctr = retreive_data(
        'ctr', method=method, OUTPUT_DIR=OUTPUT_DIR)
    best_scores_mg, best_params_mg = retreive_data(
        'mg', method=method, OUTPUT_DIR=OUTPUT_DIR)
    # - pool score
    # fig = plot_scores_pool(scores_pool_ctr, scores_pool_mg, protnames)
    # fig.savefig(os.path.join(dir, 'scores_pool.png'), dpi=300, transparent=True,
    #             facecolor='white')

    # - best param
    # fig = plot_bestparams_pool(bestparams_pool_ctr, priors, protnames)
    # fig.savefig(os.path.join(dir, 'bestparams_pool_ctr.png'), dpi=300, transparent=True, facecolor='white')
    #
    # fig = plot_bestparams_pool(bestparams_pool_mg, priors, protnames)
    # fig.savefig(os.path.join(dir, 'bestparams_pool_mg.png'), dpi=300, transparent=True, facecolor='white')

    #- best score mean
    fig = plot_scores(best_scores_ctr, best_scores_mg)
    fig.savefig(os.path.join(dir, 'besttrainscores.png'), dpi=300, transparent=True, facecolor='white')
    # # - best param
    # fig = plot_bestparams(best_params_ctr, best_params_mg, priors=priors)
    # fig.savefig(os.path.join(dir, 'bestparams.png'), dpi=300, transparent=True, facecolor='white')
