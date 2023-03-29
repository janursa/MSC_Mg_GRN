"""
    Obtains match distribution for two quantiles of .75 and .9
"""
import sys
import os
import numpy as np
import scipy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import random


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import MODELSELECTION_DIR, GRN_DIR, ENRICH_DIR, F_DE_data, F_DE_protnames, CALIBRATION_DIR
from scripts.utils import make_title_pretty
from scripts.utils.links import compare_network_string, format_links_string, normalize_links
from scripts.utils import  calibration
from geneRNI.evaluation import precision_recall_curve
def determine_sig_signes(datas):
    """ to test sig distribtion from noise links
    Conducts t test to determine whether datas[1:] are significantly different than datas[0], which is ctr
    Datas: Tuple(DataFrame), e.g. [ctr, RF, Ridge, Portia]
    """
    ctr = datas[0] #random
    #- determine p values: compared to ctr
    pvalues = np.array([])
    for data in datas[1:]:
        s, p = scipy.stats.ttest_ind(data, ctr)
        pvalues = np.append(pvalues, p)
    #- determine whether mean distribution is higher than ctr: we only plot higher ones
    increase_flags = np.array((np.mean(datas[1:], axis=1) - np.mean(ctr))>0)
    #- use p values with value flags
    def define_sign(p):
        if p:
            sign = r'$*$'
        else:
            sign=''
        return sign
    flags = (pvalues<0.05)*increase_flags
    sig_signs = ['']+[define_sign(flag) for flag in flags]
    return sig_signs

def plot_recall_dist(ax, data_stack, labels, sig_signs):
    matplotlib.rcParams.update({'font.size': 12})

    # fig, ax = plt.subplots(1, 1, tight_layout=True, figsize=(4.7,3.5),
    #     )

    bplot = ax.violinplot(data_stack, showmeans=True, showextrema=False)

    ax.set_ylabel('Number of matched interactions')
    ax.set_xticks(list(range(1,len(labels)+1)))
    ax.set_xticklabels(labels,rotation=0)
    ax.set_ymargin(.25)
    #- face colors
    colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan','grey']
    for patch, color in zip(bplot['bodies'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(1)
    #- plot sig
    xs = ax.get_xticks()
    ys = np.max(data_stack, axis=1)
    for i, sign in enumerate(sig_signs):
        if sign != '':
            ax.annotate(sign, xy=(xs[i],ys[i]),
                ha='center',
                va='bottom',
                )

def create_random_links(links_assembly, n=1000):
    #- TODO: create matrix (using protname)
    links_assembly = [normalize_links(links) for links in links_assembly]
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j] #flatten
    sample_links =  links_assembly[0]
    random_links = pd.DataFrame({key:sample_links[key] for key in ['Regulator', 'Target']})
    weightpoolvector = []
    for i in range(len(links_assembly[0])):
        weightpoolvector.append(random.sample(weights, n))
    random_links['Weight']= np.mean(weightpoolvector, axis=1)
    random_links_pool = random_links.copy()
    random_links_pool['WeightPool'] = weightpoolvector
    return random_links, random_links_pool
if __name__ == '__main__':

    if not os.path.isdir(MODELSELECTION_DIR):
        os.makedirs(MODELSELECTION_DIR)

    arbitrary_dir = os.path.join(GRN_DIR, 'arbitrary')
    if not os.path.isdir(arbitrary_dir):
        os.makedirs(arbitrary_dir)


    n_repeat = 100
    methods = ['RF', 'ridge', 'portia', 'arbitrary']
    methods_ordered = ['arbitrary', 'RF', 'ridge', 'portia'] #the order in plot
    DE_types = F_DE_data().keys() # different models, early late combined 30 50 -> 6 combination
    study = 'ctr'

    #- retreieve the string links
    links_string_dict = {}
    for DE_type in DE_types:
        links_string: pd.DataFrame = pd.read_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index_col=False)
        links_string_dict[DE_type] = links_string

    #- for plot
    ncols = 2
    nrows = int(len(DE_types) / ncols)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True, figsize=(4 * ncols, 3 * nrows))
    for idx, DE_type in enumerate(DE_types):
        protnames = F_DE_protnames()[DE_type]
        #- retreive or create the links
        links_stack = []
        for method in methods:
            if method in [ 'portia', 'ridge']:
                # - single weight
                links = pd.read_csv(os.path.join(GRN_DIR, method, f'links_{DE_type}_{study}.csv'), index_col=False)
                links_stack.append(links)
            elif method == 'RF':
                links_pool = pd.read_pickle(os.path.join(GRN_DIR, method, f'links_pool_{DE_type}_{study}.csv'))
                links_stack.append(links_pool)
            elif method == 'arbitrary':
                links, links_pool = create_random_links(links_stack, n=n_repeat)
                links.to_csv(os.path.join(arbitrary_dir, f'links_{DE_type}_{study}.csv'), index=False)
                links_pool.to_pickle(os.path.join(arbitrary_dir, f'links_pool_{DE_type}_{study}.csv'))
                links_stack.append(links_pool)
        #- create filter masks based on fit scores
        filter_masks = []
        for i_method, method in enumerate(methods):
            links = links_stack[i_method]
            if method in ['ridge','RF']:
                # best_scores, best_params = calibration.retrieve_data(study,method, DE_type, CALIBRATION_DIR)
                # protnames_left = np.asarray(protnames)[np.asarray(best_scores>0)]
                # mask = links['Target'].isin(protnames_left)
                # filter_masks.append(mask)
                filter_masks.append(None)
            else:
                filter_masks.append(None)

        #- calculate recall
        golden_links = format_links_string(links_string_dict[DE_type], protnames)
        mc_dist_list = []
        for i_method, method in enumerate(methods):
            links = links_stack[i_method]
            if method in ['portia', 'ridge']:
                precision, recall, thresholds = precision_recall_curve(links, golden_links, filter_masks[i_method])
                mc = round(np.mean(recall[:-1]), 2)
                mc_dist = np.repeat(mc, n_repeat)
            elif method in ['RF', 'arbitrary']:
                mc_dist = []
                weightpool = np.array(links['WeightPool'].values.tolist()).T
                for weight in weightpool[0:n_repeat]:
                    links['Weight'] = weight
                    precision, recall, thresholds = precision_recall_curve(links, golden_links, filter_masks[i_method])
                    mc = round(np.mean(recall[:-1]), 2)
                    mc_dist.append(mc)

            mc_dist_list.append(mc_dist)
        #- change the position of arbit to the first one
        mc_dist_list.insert(0, mc_dist_list[-1])
        del mc_dist_list[-1]


        sig_signs = determine_sig_signes(mc_dist_list)
        i = int(idx/ncols)
        j = idx%ncols
        ax = axes[i][j]
        ax.set_title(make_title_pretty(DE_type))
        plot_recall_dist(ax=ax,data_stack=mc_dist_list, labels=methods_ordered, sig_signs=sig_signs)

    fig.savefig(os.path.join(MODELSELECTION_DIR, f'evaluation.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(MODELSELECTION_DIR, f'evaluation.pdf'))






