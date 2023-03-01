"""
    Sets of functions useful for processing network inference
"""
import sys
import os
import statistics
import scipy
import matplotlib
import matplotlib.pyplot as plt
import copy
import numpy as np
import random
import pandas as pd
import json


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from scripts.imports import *
from scripts.utils import create_check_dir, serif_font

from geneRNI.core import network_inference
from geneRNI.data import Data

def batch_GRN(study, method, DE_type, i_start, i_end, output_dir, **specs):
    """
        Runs network inference for a number of times
    """
    #- create/check dirs for method, study, (links and scores)
    DIR_METHOD = create_check_dir(output_dir, method)
    DIR_DE_type = create_check_dir(DIR_METHOD, DE_type)
    DIR_STUDY = create_check_dir(DIR_DE_type, study)
    DIR_LINKS = create_check_dir(DIR_STUDY, 'links')
    DIR_TESTSCORES = create_check_dir(DIR_STUDY, 'testscores')
    DIR_TRAINSCORES = create_check_dir(DIR_STUDY, 'trainscores')
    #- run iteration
    for i in range(i_start, i_end):
        print(f'----GRN for {study} iteration {i}-----')
        ests, train_scores, links_df, oob_scores, test_scores = grn(**specs)
        if method=='RF':
            test_scores = oob_scores
        #- save
        np.savetxt(os.path.join(DIR_TESTSCORES, f'data_{i}.csv'), test_scores)
        np.savetxt(os.path.join(DIR_TRAINSCORES, f'data_{i}.csv'), train_scores)
        links_df.to_csv(os.path.join(DIR_LINKS, f'data_{i}.csv'), index=False)

    #------ pools iteration results such links and scores ---
    testscores_stack = []
    trainscores_stack = []
    links_stack = []
    for i in range(i_end):
        # - retreive
        links = pd.read_csv(os.path.join(DIR_LINKS, f'data_{i}.csv'), index_col=False)
        testscores = np.genfromtxt(os.path.join(DIR_TESTSCORES, f'data_{i}.csv'))
        trainscores = np.genfromtxt(os.path.join(DIR_TRAINSCORES, f'data_{i}.csv'))
        # - stack them
        testscores_stack.append(testscores)
        trainscores_stack.append(trainscores)
        links_stack.append(links)
    # - pool the weights and  create average links df
    ws_pool = np.array([ll['Weight'].values for ll in links_stack]).T
    ws_mean = np.mean(ws_pool, axis=1)
    links = pd.DataFrame()
    links.loc[:, ['Regulator', 'Target']] = links_stack[0].loc[:, ['Regulator', 'Target']]
    links_mean = links.copy()
    links_pool = links.copy()
    links_mean['Weight'] = ws_mean
    links_pool['WeightPool'] = list(ws_pool)
    links_pool['Weight'] = ws_mean
    # - average pool scores
    testscores_mean = np.mean(testscores_stack, axis=1)
    trainscores_mean = np.mean(trainscores_stack, axis=1)
    # - save
    links_mean.to_csv(os.path.join(DIR_METHOD, f'links_{DE_type}_{study}.csv'))
    links_pool.to_pickle(os.path.join(DIR_METHOD, f'links_pool_{DE_type}_{study}.csv'))

    np.savetxt(os.path.join(DIR_DE_type, f'testscores_{study}.csv'), testscores)
    np.savetxt(os.path.join(DIR_DE_type, f'trainscores_{study}.csv'), trainscores)

def retrieve_scores(study, method, DE_type, output_dir):
    """
        Retreiev train and test scores
    """
    testscores = np.genfromtxt(os.path.join(output_dir, method, DE_type, f'testscores_{study}.csv'))
    trainscores = np.genfromtxt(os.path.join(output_dir, method, DE_type, f'trainscores_{study}.csv'))

    return trainscores, testscores
def grn(data, time_points, gene_names, **specs):
    dataset = Data(ts_data=[data], ss_data=None, time_points=[time_points], gene_names=gene_names)
    return network_inference(dataset, gene_names=gene_names, **specs)

def compare_network_string(DE_type, links, top_quantile, enrich_output_dir, verbose=False) -> int:
    '''
        Compare extracted links by GRN to those suggested by model_selection. Returns number of match links.
    '''
    links_short = choose_top_quantile(links, quantile=top_quantile)
    links_string = pd.read_csv(os.path.join(enrich_output_dir, f'network_{DE_type}.csv'), index_col=False)
    links_string.rename(columns={'node1':'Regulator','node2':'Target','score':'Weight'}, inplace=True)
    links_string = links_string.loc[:,['Regulator','Target','Weight']]
    if verbose:
        print(f'Number of model_selection links: {len(links_string)}')
        print(f'Number of extracted links: {len(links_short)}')
    if len(links_string)>len(links_short):
        print('Extracted links cannot be lesser than golden links')
        print('Golden links are shortened')
        links_string.sort_values('Weight', ascending=False)
        links_string = links_string.head(len(links_short))
        # raise ValueError('Extracted links cannot be lesser than golden links')

    #- find model_selection links in the extracted links, label them and add stringweight
    links_short['inString'] = False
    for reg, target, weight in zip(links_string['Regulator'].values, links_string['Target'].values, links_string['Weight'].values):
        links_short.loc[(links_short['Regulator']==reg) & (links_short['Target']==target),'inString'] = True
        links_short.loc[(links_short['Regulator']==reg) & (links_short['Target']==target),'StringWeight'] = weight
    n = links_short['inString'].values.tolist().count(True)
    n_norm = n/len(links_string)

    return n_norm
def compare_network_string(DE_type, links, top_quantile, enrich_output_dir, verbose=False) -> int:
    '''
        Compare extracted links by GRN to those suggested by model_selection. Returns number of match links.
    '''
    links_short = choose_top_quantile(links, quantile=top_quantile)
    links_string = pd.read_csv(os.path.join(enrich_output_dir, f'network_{DE_type}.csv'), index_col=False)

    links_short['inString'] = False
    for reg, target, weight in zip(links_string['Regulator'].values, links_string['Target'].values, links_string['Weight'].values):
        links_short.loc[(links_short['Regulator']==reg) & (links_short['Target']==target),'inString'] = True
        links_short.loc[(links_short['Regulator']==reg) & (links_short['Target']==target),'StringWeight'] = weight
    n = links_short['inString'].values.tolist().count(True)
    n_norm = n/len(links_string)

    return n_norm

def _read_write_links(method, study, mode:str, links=None, output_dir='') -> pd.DataFrame:
    '''
        Read write links extracted from GRN 
    '''
    assert(study in ['ctr','mg','all-in','random'])
    assert(mode in ['read', 'write'])
    #- determine file location
    DIR = os.path.join(output_dir, method)
    FILE = os.path.join(DIR, f'links_{study}.csv')

    if mode=='read':
        return pd.read_csv(FILE, index_col=False)
    else:
        assert(links is not None)
        links.to_csv(FILE, index=False)


def pool_links(study, protnames, output_dir, n, method='') -> pd.DataFrame:
    '''
        Read the links from different GRN runs and create a pool of weights. It also creates the average weights as a database.

    '''
    # def fake_it(links):
    #     links_f = copy.deepcopy(links)
    #     n = len(links_f.loc[:,'Weight'])
    #     mean = np.mean(links_f.loc[:,'Weight'])
    #     links_f.loc[:,'Weight'] = links_f.loc[:,'Weight']+[mean*random.random() for i in range(n)]
    #     return links_f
    # links_temp = read_write_links(study=study, mode='read', OUTPUT_DIR=output_dir)
    # links_list = [fake_it(links_temp) for i in range(n)]
    
    # links_list = [read_write_links(study=study, mode='read', i=i, OUTPUT_DIR=output_dir, method=method) for i in range(n)]

    #- pool the weights 
    ws_pool = np.array([ll['Weight'].values for ll in links_list]).T
    #- create average links df
    ws_mean = np.mean(ws_pool, axis=1)
    links = pd.DataFrame()
    links.loc[:,['Regulator','Target']] = links_list[0].loc[:,['Regulator','Target']]
    links['Weight'] = ws_mean
    links['WeightPool'] = list(ws_pool)
    return links
def _pool_GRN_oo(method, study, replica_n, output_dir):
    """ pools iteration results such links and scores
    """
    DIR_METHOD = os.path.join(output_dir, 'GRN', method)
    DIR_STUDY = os.path.join(DIR_METHOD, study)
    DIR_LINKS = os.path.join(DIR_STUDY, 'links')
    DIR_TRAINSCORES = os.path.join(DIR_STUDY, 'trainscores')
    DIR_TESTSCORES = os.path.join(DIR_STUDY, 'testscores')
    testscores_stack = []
    trainscores_stack = []
    links_stack = []
    for i in range(replica_n):
        #- retreive
        links = pd.read_csv(os.path.join(DIR_LINKS, f'data_{i}.csv'), index_col=False)
        testscores = np.genfromtxt(os.path.join(DIR_TESTSCORES, f'data_{i}.csv'))
        trainscores = np.genfromtxt(os.path.join(DIR_TRAINSCORES, f'data_{i}.csv'))
        #- stack them
        testscores_stack.append(testscores)
        trainscores_stack.append(trainscores)
        links_stack.append(links)
    # - pool the weights and  create average links df
    ws_pool = np.array([ll['Weight'].values for ll in links_stack]).T
    ws_mean = np.mean(ws_pool, axis=1)
    links = pd.DataFrame()
    links.loc[:, ['Regulator', 'Target']] = links_stack[0].loc[:, ['Regulator', 'Target']]
    links_mean = links.copy()
    links_pool = links.copy()
    links_mean['Weight'] = ws_mean
    links_pool['WeightPool'] = list(ws_pool)
    links_pool['Weight'] = ws_mean
    #- average pool scores
    testscores_mean = np.mean(testscores_stack, axis=1)
    trainscores_mean = np.mean(trainscores_stack, axis=1)
    #- save
    links_mean.to_csv(os.path.join(DIR_METHOD, f'links_{study}.csv'))
    links_pool.to_pickle(os.path.join(DIR_METHOD, f'links_pool_{study}.csv'))

    write_scores(method=method, study=study, trainscores=trainscores_mean,
                 testscores=testscores_mean, output_dir=DIR_METHOD)


def _write_scores(method, study, trainscores, testscores , output_dir):
    """
    Write train test scores to files
    """
    np.savetxt(os.path.join(output_dir, f'testscores_{study}.csv'), testscores)
    np.savetxt(os.path.join(output_dir, f'trainscores_{study}.csv'), trainscores)
def filter_ttest(links) ->pd.DataFrame:
    '''
        Conducts t test and select those significanly different than 0
    '''
    
    def t_test(vector):
        mode = statistics.mode(vector)
        m = len(vector)
        std = np.std(vector)
        t_value = mode/(std/np.sqrt(m))
        p_value = scipy.stats.t.ppf(q=1-.05,df=m-1)
        return t_value, p_value
    # tvalue = np.array([t_test(item)[0] for item in links['WeightPool']])
    # pvalue = np.array([t_test(item)[1] for item in links['WeightPool']])

    tvalue = np.array([scipy.stats.ttest_1samp(item,0)[0] for item in links['WeightPool']])
    pvalue = np.array([scipy.stats.ttest_1samp(item,0)[1] for item in links['WeightPool']])

    links['pvalue'] = pvalue
    links['tvalue'] = tvalue
    #- keep those with sig less than 0.05
    links = links.loc[pvalue<0.05,:]
    print(f'number of links after ttest {len(links)}')
    return links
def _filter_fitscore(links) -> pd.DataFrame:
    '''
        Filters the given df based on the fitscore
    '''
    #- keep those with FitScore over 0
    links = links.loc[links['FitScore']>0,:].reset_index(drop=True)
    # print(f'number of links after fitscore {len(links)}')
    return links
def choose_top_quantile(links: pd.DataFrame, quantile=0.75)->pd.DataFrame:
    '''
        Filters the given df based on the top quantile
    ''' 
    # links.reset_index(inplace=True, drop=True)
    # links.sort_values('Weight',ascending=True,inplace=True)
    cut_off = np.quantile(links['Weight'].values.tolist(), q=quantile)

    links_short = links.loc[links['Weight']>=cut_off,:].reset_index(drop=True)
    return links_short
def choose_top_count(links: pd.DataFrame, n=100)->pd.DataFrame:
    '''
        Filters the given df based on the top count
    ''' 
    links.reset_index(inplace=True, drop=True)
    links.sort_values('Weight',ascending=False,inplace=True)
    links_short = links.iloc[:n,:].reset_index(drop=True)
    return links_short
def plot_mean_weights(ax, links, name, studies):
    serif_font()
    colors = ['lightblue', 'pink']
    for i,study in enumerate(links):
        ax.hist(study['Weight'], bins=50, alpha=0.6,
                        histtype='bar', #'bar', 'barstacked', 'step', 'stepfilled'
                        color=colors[i],
                        # ec='black',
                        rwidth=1.5,
                        # density = True
                       )
    handles = []

    for i, color in enumerate(colors):
        handles.append(ax.scatter([],[], marker='o', label=studies[i],
         edgecolor='black', color=color, linewidth=.2))

        ax.legend(handles=handles,
            bbox_to_anchor=(1,1), prop={'size': 9}
            # title='Enriched Term'
            )
        # ax.add_artist(ll)
        ax.set_xlabel('Interaction strength')
        ax.set_ylabel('Density')

        ax.set_ymargin(.2)
#         ax.set_xmargin(.1)
        # ax.set_xlim([-.5,8])
        ax.set_title(name)

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
def compare_network_string_batch(DE_type, links, top_quantile, n_repeat, enrich_output_dir):
    """
    Compare the given links to model_selection for each weight set in weightpool
    """
    match_counts = []
    weightpool = np.array(links['WeightPool'].values.tolist()).T
    for weight in weightpool[0:n_repeat]:
        links['Weight'] = weight
        match_count = compare_network_string(DE_type=DE_type,links=links, top_quantile=top_quantile, enrich_output_dir=enrich_output_dir)
        match_counts.append(match_count)
    return np.array(match_counts)
def format_links_string(links, gene_names) -> pd.DataFrame:
    """Converts links to golden links style
    Links are given as gene 1 to gene 2 with a weight.
    Golden link is combination of all genes with 0 or 1 for the true links.
    We assume that any link given in links is 1 disgarding quantity of Weight.
    """
    n_genes = len(gene_names)
    regs = np.repeat(gene_names, n_genes-1)
    targs = []
    for i in range(n_genes):
        targs.extend(gene_names[:i] + gene_names[i+1:])
    golden_links = pd.DataFrame({'Regulator':regs, 'Target': targs})
    golden_links['Weight'] = (golden_links['Regulator'].isin(links['Regulator']) & golden_links['Target'].isin(links['Target'])).astype(int)
    #- control
    if True:
        gl = golden_links
        assert(len(golden_links)==n_genes*(n_genes-1))
        for i in range(10): # for 10 random selection, golden link should match the links
            idx = random.randint(0, len(links)-1)
            reg, targ,_ = links.iloc[idx,:]
            should_be_one = gl.loc[(gl['Regulator'] == reg) & (gl['Target'] == targ), 'Weight'].iloc[0]
            assert isinstance(should_be_one, np.int32)
            assert (should_be_one==1)


        print('Golden links creation control successful')


    return golden_links
def remove_mg(links):
    return links.loc[(links['Regulator']!='mg') & (links['Target']!='mg'),:]

def plot_match_counts_series(match_counts_list, links_names, top_quantile_list, ax=None):
    """
    Plots match counts for different top selected links as a line plot. Each line is for different methods such as RF and Ridge
    """
    top_quantile_list = np.round(top_quantile_list, 2)
    match_counts_sum = np.round(np.sum(match_counts_list, axis=1),2) #final scores
    # add scores to the names
    links_names = [f'{name}: (score={score})' for name, score in zip(links_names, match_counts_sum)]
    serif_font()
    if ax is None:
        fig, ax = plt.subplots(1, 1, tight_layout=True, figsize=(4.7, 3.5))
    colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey', 'blue']
    linestyles = ['-', '--', '-.', ':', '-', '--']
    assert (len(match_counts_list) == len(links_names))
    for i, data in enumerate(match_counts_list):
        ax.plot(top_quantile_list, data, label=links_names[i], color=colors[i], alpha=1, linewidth=2,
                linestyle=linestyles[i],
                marker='o')
    ax.set_ylabel('Fraction of network match')
    ax.set_xlabel('Selected top quantile')
    xticks = np.arange(min(top_quantile_list), max(top_quantile_list), .05)
    ax.set_xticks(xticks)
    # ax.set_xmargin(.15)
    # ax.set_ymargin(.25)
    ax.set_ylim([-.01,.25])

    ax.legend(frameon=False)
    # - annotate score on the right side of each line
    
    # for i, (score, y) in enumerate(zip(match_counts_sum, match_counts_list)):
    #     ax.annotate(score, xy=(top_quantile_list[-1] + .01, y[-1] - 0.5),
    #                 ha='center',
    #                 va='bottom',
    #                 size=10,
    #                 alpha=.5
    #
    #                 )

    try:
        return fig
    except:
        pass

def normalize_links(links):
    """
        Nornalize the links based on the std
    """
    links_n = links.copy()
    links_n.loc[:,'Weight'] = links['Weight']/np.std(links['Weight'])
    return links_n

def _convert_links_to_nodes_edges(links, protnames):
    '''
        Takes links as pd and seperates nodes and edges
    '''
    edges = links.copy()
    edges.rename(columns={'Regulator':'source', 'Target':'target'}, inplace=True)
    sum_ws = [sum(links.query(f"Regulator == '{gene}'")['Weight']) for gene in protnames]
    nodes = pd.DataFrame(data={'Protein':protnames,'Weight':sum_ws})
    print(edges)
    aa
    return nodes, edges

    
def read_write_nodes_edges(nodes=None, edges=None, study='ctr', mode='read', OUTPUT_DIR=''):
    '''
        Manages links for network visualization (nodes and edges)
    '''
    DIR = os.path.join(OUTPUT_DIR, 'GRN')
    edges_FILE = os.path.join(DIR, f'edges_{study}.csv')
    nodes_FILE = os.path.join(DIR, f'nodes_{study}.csv')
    if mode == 'write':
        edges.to_csv(edges_FILE,index=False)
        nodes.to_csv(nodes_FILE,index=False)
        print('successfully wrote nodes and edges to ', DIR)
    elif mode =='read':
        edges = pd.read_csv(edges_FILE, index_col=False)
        nodes = pd.read_csv(nodes_FILE, index_col=False)
        print('successfully read nodes and edges from ', DIR)
        return nodes, edges
# def plot_scores(data_ctr, data_sample, xlabel=''):
#     """plots oob adn train scores as a box plot for ctr and mg side by side"""
#     serif_font()
#     fig, axes = plt.subplots(1, 1, tight_layout=True, figsize=(2.5, 3),
#                              # gridspec_kw={'width_ratios': [2, 2]}
#                              )
#     data_s = [data_ctr, data_sample]
#     labels = ['ctr', 'mg']
#     ax = axes
#     for data in data_s:
#         bplot = ax.boxplot(data, notch=True, widths=[.5, .5], patch_artist=True, meanline=True)
#     # bplot = ax.violinplot(data_s, showmeans=True, showextrema=True, bootstrap=True
#     #     )
#     ax.set_ylabel('Score')
#     ax.set_xticks(range(1, len(labels) + 1))
#     ax.set_xticklabels(labels, rotation=0)
#     ax.axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=.5)
#     ax.set_ymargin(.1)
#     ax.set_xmargin(.15)
#     # # - face colors
#     # colors = ['pink', 'lightblue', 'lightgreen', 'cyan']
#     # for patch, color in zip(bplot['boxes'], colors):
#     #     patch.set_facecolor(color)
#     #     patch.set_edgecolor('black')
#     #     patch.set_alpha(1)
#
#     # colors = ['black' for i in range(len(labels))]
#     # tags = ['cbars']
#     # for tag in tags:
#     #     bplot[tag].set_color(colors)
#     #     bplot[tag].set_decay_coeff(.5)
#     return fig