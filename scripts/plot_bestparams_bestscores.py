"""
Reads the results of processed hyper param tuning (pooled and averaged) and plot them.
"""
import os
import pathlib
import sys
dir_main = os.path.join(pathlib.Path(__file__).parent.resolve(), '..')
sys.path.insert(0, dir_main)  # TODO: not recommended (let's make a setup.py file instead)

from scripts.imports import *

utils.serif_font()

#- read pooled and average data of best params and best scores
def retreive_data(study):
    with open(os.path.join(OUTPUT_DIR,'calibration', f'oo_{study}.txt')) as f:
        oo = eval(f.read())
        bestparams, bestscores = oo['bestparams'], oo['bestscores']
    with open(os.path.join(OUTPUT_DIR,'calibration', f'oo_pool_{study}.txt')) as f:
        oo = eval(f.read())
        bestparams_pool, bestscores_pool = oo['bestparams_pool'], oo['bestscores_pool']
    return bestparams, bestscores, bestparams_pool, bestscores_pool

bestparams_ctr, bestscores_ctr, bestparams_pool_ctr, bestscores_pool_ctr = retreive_data('ctr')
bestparams_mg, bestscores_mg, bestparams_pool_mg, bestscores_pool_mg = retreive_data('mg')

#- plots
def plot_bestscores_pool():
    """plots scores as a box plot for a set"""
    fig, axes = plt.subplots(2, 1, tight_layout=True, figsize=(10,6))
    data_s = [bestscores_pool_ctr, bestscores_pool_mg]
    titles = ['ctr','mg']
    for i, data in enumerate(data_s):
        if data is None:
            continue
        ax = axes[i]
        # ax.boxplot(np.array(data).T)
        bplot=ax.violinplot(data, showmeans=True, showextrema=False)
        ax.set_ylabel('Testing score')
        ax.set_xlabel('')
        ax.set_title(titles[i])
        ax.set_xticks(range(1, len(protnames) + 1))
        ax.set_xticklabels(protnames,rotation=90)
        ax.set_ymargin(.1)
        ax.set_xmargin(.02)
        # colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan','grey']
        for patch in bplot['bodies']:
            patch.set_facecolor('lightgreen')
            patch.set_edgecolor('black')
            patch.set_decay_coeff(1)
    return fig

def plot_bestparams_pool(study):
    """
        Plots boxplot for indivual param in seperate window. In each window, the variation in 
        best params are given for each gene seperately, obtained during different runs of tunning.
    """
    nrows = len(param_grid)
    ncols = 1
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(10*ncols,9))
    if study == 'ctr':
        data = bestparams_pool_ctr
    else:
        data = bestparams_pool_mg
    for i, key in enumerate(param_grid.keys()):

        ax = axes[i]
        pool = [item[key] for item in data]
        bplot = ax.violinplot(pool, showmeans=True, showextrema=False)
        ax.set_title(key)
        ax.set_ylabel('Quantity')
        ax.set_xticks(range(1, len(protnames) + 1))
        ax.set_xticklabels(protnames,rotation=90)
        ax.set_ymargin(.1)
        ax.set_xmargin(.02)
        # colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan','grey']
        for patch in bplot['bodies']:
            patch.set_facecolor('lightgreen')
            patch.set_edgecolor('black')
            patch.set_decay_coeff(1)
    return fig
def plot_bestscores(data_ctr, data_sample, xlabel=''):

    """plots scores as a box plot for ctr and mg side by side"""
    serif_font()
    fig, axes = plt.subplots(1, 1, tight_layout=True, figsize=(2.5,3), 
        # gridspec_kw={'width_ratios': [2, 2]}
        )
    data_s = [data_ctr, data_sample]
    labels = ['ctr','mg']
    ax = axes
    bplot = ax.boxplot(data_s, notch=True, widths =[.5,.5], patch_artist=True, meanline=True)
    # bplot = ax.violinplot(data_s, showmeans=True, showextrema=True, bootstrap=True
    #     )
    ax.set_ylabel('OOB Score')
    # ax.set_title('(A)')
    ax.set_xticks(range(1,len(labels)+1))
    ax.set_xticklabels(labels, rotation=0)
    ax.axhline(0,color='red', linestyle='--',linewidth=1.5, decay_coeff=.5)
    ax.set_ymargin(.1)
    ax.set_xmargin(.15)
    #- face colors
    colors = ['pink', 'lightblue', 'lightgreen', 'cyan']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_decay_coeff(1)

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
    serif_font()
    fig, axes = plt.subplots(1, 2, tight_layout=True, figsize=(7,4))
    datas = [data_ctr,data_sample]
    titles = ['ctr', 'mg']
    def normalize(xx, priors):
        xx = {key: np.array(list(set(values))) for key, values in xx.items()}
        return {key: (values - min(priors[key])) / (max(priors[key]) - min(priors[key])) for key, values in
                    xx.items()}
    
    for i, data in enumerate(datas):
        data = {key:[item[key] for item in data] for key in data[0].keys()}
        data = normalize(data, priors)
        bplot = axes[i].boxplot(data.values(), notch=False, widths =.4, patch_artist=True, meanline=True)

        axes[i].set_title(titles[i])
        axes[i].set_ylabel('Normalized quantity')
        axes[i].set_xticklabels(data.keys(), rotation=45)
        axes[i].set_ymargin(.15)
        axes[i].set_xmargin(.15)
        colors = ['c', 'olive', 'lightgreen','g']
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)


    return fig


if __name__ == '__main__':
    #- best score pool
    fig = plot_bestscores_pool()
    fig.savefig(os.path.join(OUTPUT_DIR,'calibration','bestscores_pool.png'), dpi=300, transparent=True, facecolor='white')
    #- best param pool
    fig = plot_bestparams_pool('ctr')
    fig.savefig(os.path.join(OUTPUT_DIR,'calibration','bestparams_pool_ctr.png'), dpi=300, transparent=True, facecolor='white')

    fig = plot_bestparams_pool('mg')
    fig.savefig(os.path.join(OUTPUT_DIR,'calibration','bestparams_pool_mg.png'), dpi=300, transparent=True, facecolor='white')

    #- best score mean
    fig = plot_bestscores(bestscores_ctr, bestscores_mg)
    fig.savefig(os.path.join(OUTPUT_DIR,'calibration','bestscores.png'), dpi=300, transparent=True, facecolor='white')
    #- best param mean
    fig = plot_bestparams(bestparams_ctr, bestparams_mg, priors=param_grid)
    fig.savefig(os.path.join(OUTPUT_DIR,'calibration','bestparams.png'), dpi=300, transparent=True, facecolor='white')
