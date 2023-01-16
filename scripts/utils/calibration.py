"""
    Sets of functions useful for hyperparameter tuning
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
import utils

from ._imports import *
def read_write_oo(study='ctr', mode='read', method='RF', bestparams=None, besttestscores=None, besttrainscores=None, i=None, OUTPUT_DIR=''):
    '''
        Read and writes calibration results (best params and best scores) to file
    '''
    assert (mode in ['read', 'write'])
    assert (study in ['ctr', 'mg', 'combined'])
    DIR = os.path.join(OUTPUT_DIR, 'calibration', str(method))
    if i is None:
        FILE = os.path.join(DIR, f'oo_{study}.txt')
    else:
        FILE = os.path.join(DIR, 'pool', f'oo_{study}_{i}.txt')
    if mode == 'write':
        with open(FILE, 'w') as f:
            print({'bestparams': bestparams, 'besttestscores': besttestscores, 'besttrainscores':besttrainscores}, file=f)
    elif mode == 'read':
        with open(FILE, 'r') as f:
            oo = eval(f.read())
        return oo['bestparams'], oo['besttrainscores'], oo['besttestscores']
def main(repeat_n, study, method, data, gene_names, time_points, test_size, param, param_grid, **specs):
    """
        Interface to search function of geneRNI
    """
    print(f'----Tuning for method: {method} study: {study}')
    dataset = Data(gene_names=gene_names, ss_data=None, ts_data=[data], time_points=[time_points],
                       test_size=test_size)
    DIR = os.path.join(OUTPUT_DIR, 'calibration', method, study)
    #TODO: create the DIR if it doesnt exist
    return search_param.rand_search(dataset,
                                    param=param,
                                    param_grid=param_grid,
                                    repeat_n=repeat_n,
                                    output_dir=DIR,
                                    **specs)
def batch_tune(study, method, data, param, gene_names, test_size, time_points, param_grid, istart=None, iend=None,OUTPUT_DIR='', **specs):


    # read_write_oo(study, method=method, mode='write', bestparams=best_params, besttestscores=best_testscores, besttrainscores=best_trainscores,
    #                     OUTPUT_DIR=OUTPUT_DIR, i=i)
def pool_oo(study, method, n=100, OUTPUT_DIR=''):
    '''
        Read the results of hyper parameter tunning and pool them.
    '''
    best_paramss = []
    best_trainscoress = []
    best_testscoress = []
    for i in range(n):
        best_params, best_trainscores, best_testscores = read_write_oo(study, method=method, mode='read', OUTPUT_DIR=OUTPUT_DIR, i=i)
        best_paramss.append(best_params)
        best_trainscoress.append(best_trainscores)
        best_testscoress.append(best_testscores)

    # - process best scores
    besttrainscores_pool = np.array(best_trainscoress).T.tolist()  # n_genes*n
    besttrainscores = np.mean(besttrainscores_pool, axis=1).tolist()  # n_genes

    besttestscores_pool = np.array(best_testscoress).T.tolist()  # n_genes*n
    besttestscores = np.mean(besttestscores_pool, axis=1).tolist()  # n_genes

    # - process best params
    bestparams_pool = []  # n_gene*best_params (n values)
    for prot_i in range(len(best_paramss[0])):
        bestparams = {}
        for key in best_paramss[0][0].keys():
            vector = []
            for item in best_paramss:
                vector.append(item[prot_i][key])
            bestparams[key] = vector
        bestparams_pool.append(bestparams)
    # - average best params for each gene
    bestparams = [{key: np.mean(vector) for key, vector in bestparams.items()} for bestparams in bestparams_pool]
    # - make some param int
    bestparams = [{key: int(value) if (key == 'max_depth' or key == 'min_samples_leaf' or  key == 'max_features') else value for key, value in
                   bestparam.items()} for bestparam in bestparams]

    with open(os.path.join(OUTPUT_DIR, 'calibration', method, f'oo_{study}.txt'), 'w') as f:
        print({'bestparams': bestparams, 'besttrainscores': besttrainscores, 'besttestscores': besttestscores}, file=f)

    with open(os.path.join(OUTPUT_DIR, 'calibration', method, f'oo_pool_{study}.txt'), 'w') as f:
        print({'bestparams_pool': bestparams_pool, 'besttrainscores_pool': besttrainscores_pool,'besttestscores_pool': besttestscores_pool}, file=f)
def retreive_data(study, method, read_dir):
    """
        Reads results of pooling (#process_pool)
    """
    with open(os.path.join(read_dir, f'oo_{study}.txt')) as f:
        oo = eval(f.read())
        bestparams, besttrainscores, besttestscores  = oo['bestparams'], oo['besttrainscores'], oo['besttestscores']
    with open(os.path.join(read_dir, f'oo_pool_{study}.txt')) as f:
        oo = eval(f.read())
        bestparams_pool, besttrainscores_pool, besttestscores_pool = oo['bestparams_pool'], oo['besttrainscores_pool'], oo['besttestscores_pool']
    return bestparams,  besttrainscores, besttestscores, bestparams_pool, besttrainscores_pool, besttestscores_pool
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
def plot_oo(method, priors, protnames, OUTPUT_DIR):
    """
     Plots a series of graphs for best params and best scores (individual protein and combined)
    """
    dir = os.path.join(OUTPUT_DIR, 'calibration', method)
    bestparams_ctr, besttrainscores_ctr, besttestscores_ctr, bestparams_pool_ctr, besttrainscores_pool_ctr, besttestscores_pool_ctr = retreive_data(
        'ctr', method=method, read_dir=dir)
    bestparams_mg, besttrainscores_mg, besttestscores_mg, bestparams_pool_mg, besttrainscores_pool_mg, besttestscores_pool_mg = retreive_data('mg',
                                                                                                           method=method,
                                                                                                           read_dir=dir)

    # - best pool train score
    fig = plot_scores_pool(besttrainscores_ctr, besttrainscores_mg, protnames)
    fig.savefig(os.path.join(dir, 'besttrainscores_pool.png'), dpi=300, transparent=True,
                facecolor='white')
    # - best pool train score
    fig = plot_scores_pool(besttestscores_ctr, besttestscores_mg, protnames)
    fig.savefig(os.path.join(dir, 'besttestscores_pool.png'), dpi=300, transparent=True,
                facecolor='white')
    # - best param pool
    fig = plot_bestparams_pool(bestparams_pool_ctr, priors, protnames)
    fig.savefig(os.path.join(dir, 'bestparams_pool_ctr.png'), dpi=300, transparent=True, facecolor='white')

    fig = plot_bestparams_pool(bestparams_pool_mg, priors, protnames)
    fig.savefig(os.path.join(dir, 'bestparams_pool_mg.png'), dpi=300, transparent=True, facecolor='white')

    # - best train score mean
    fig = plot_scores(besttrainscores_ctr, besttrainscores_mg)
    fig.savefig(os.path.join(dir, 'besttrainscores.png'), dpi=300, transparent=True, facecolor='white')
    # - best test score mean
    fig = plot_scores(besttestscores_ctr, besttestscores_mg)
    fig.savefig(os.path.join(dir, 'besttestscores.png'), dpi=300, transparent=True, facecolor='white')
    # - best param mean
    fig = plot_bestparams(bestparams_ctr, bestparams_mg, priors=priors)
    fig.savefig(os.path.join(dir, 'bestparams.png'), dpi=300, transparent=True, facecolor='white')