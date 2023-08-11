"""
    Sets of functions useful for processing network inference
"""
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import portia as pt
from typing import List, Tuple, Dict, Optional
from matplotlib.figure import Figure


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from common_tools import create_check_dir, serif_font

from geneRNI.core import network_inference
from geneRNI.data import Data
from geneRNI.links import format_links

def batch_run_generni(study:str, method:str, DE_type:str, i_start:int, i_end:int, output_dir:str, **specs) -> None:
    """Runs GRN inference using geneRNI with multiple repeatitions
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
        ests, train_scores, links_df, oob_scores, test_scores = run_generni(**specs)
        if method=='RF':
            test_scores = oob_scores
        #- save
        np.savetxt(os.path.join(DIR_TESTSCORES, f'data_{i}.csv'), test_scores)
        np.savetxt(os.path.join(DIR_TRAINSCORES, f'data_{i}.csv'), train_scores)
        links_df.to_csv(os.path.join(DIR_LINKS, f'data_{i}.csv'), index=False)
    #------ pools iteration results, i.e. links and scores ---
    testscores_stack = []
    links_stack = []
    for i in range(i_end):
        # - retreive
        links = pd.read_csv(os.path.join(DIR_LINKS, f'data_{i}.csv'), index_col=False)
        testscores: np.adarray = np.genfromtxt(os.path.join(DIR_TESTSCORES, f'data_{i}.csv'))
        # - stack them
        testscores_stack.append(testscores)
        links_stack.append(links)
    # - pool the weights and  create average links df
    ws_pool = np.array([ll['Weight'].values for ll in links_stack]).T
    ws_mean = np.mean(ws_pool, axis=1)
    links = pd.DataFrame()
    links.loc[:, ['Regulator', 'Target']] = links_stack[0][['Regulator', 'Target']]
    links_mean = links.copy()
    links_pool = links.copy()
    links_mean['Weight'] = ws_mean
    links_pool['WeightPool'] = list(ws_pool)
    links_pool['Weight'] = ws_mean
    # - save links
    to_save = os.path.join(DIR_METHOD, f'links_{DE_type}_{study}.csv')
    links_mean.to_csv(to_save)
    print(f'output -> {to_save}')
    links_pool.to_pickle(os.path.join(DIR_METHOD, f'links_pool_{DE_type}_{study}.csv'))
    # - average test scores for multiple runs
    testscores = np.mean(np.asarray(testscores_stack), axis=0)
    np.savetxt(os.path.join(DIR_DE_type, f'testscores_{study}.csv'), testscores)
def retrieve_scores(study: str, method:str, DE_type:str, output_dir:str) -> Tuple[list[float], list[float]]:
    """Retrieves train/test scores
    """
    testscores = list(np.genfromtxt(os.path.join(output_dir, method, DE_type, f'testscores_{study}.csv')))
    trainscores = list(np.genfromtxt(os.path.join(output_dir, method, DE_type, f'trainscores_{study}.csv')))

    return trainscores, testscores
def run_portia(data, gene_names, **kwargs):
    """Runs GRN inference using Portia

    """
    data = np.asarray(data)
    # - create portia dataset
    portia_dataset = pt.GeneExpressionDataset()
    for exp_id, data_i in enumerate(data):
        portia_dataset.add(pt.Experiment(exp_id, data_i))
    # - GRN inference
    M_bar = pt.run(portia_dataset, method='fast', verbose=False, normalize=True)
    links_df = format_links(M_bar, gene_names)
    return links_df

def run_generni(data, time_points, gene_names, **specs):
    """Runs GRN using geneRNI
    """
    dataset = Data(ts_data=[data], ss_data=None, time_points=[time_points], gene_names=gene_names)
    return network_inference(dataset, gene_names=gene_names, **specs)

def choose_top_quantile(links: pd.DataFrame, quantile=0.75, col_name='Weight') -> pd.DataFrame:
    '''
        Filters the given df based on the top quantile
    '''
    cut_off = np.quantile(links[col_name].values.tolist(), q=quantile)
    links_short = links.loc[links[col_name].tolist()>=cut_off,:].reset_index(drop=True)
    return links_short
def choose_top_count(links: pd.DataFrame, n=100, col_name: str='Weight') -> pd.DataFrame:
    '''
        Filters the given df based on the top count
    ''' 
    links.reset_index(inplace=True, drop=True)
    links.sort_values(col_name,ascending=False,inplace=True)
    links_short = links.iloc[:n,:].reset_index(drop=True)
    return links_short

def plot_grn_weights_distributions_per_top_links(links:pd.DataFrame, color:str, name:str, dist_key:str='WeightPool') -> Figure:
    """ Plots weights distributions for the top g-g interactions

    """
    serif_font()
    links = links.sort_values('Weight', ascending=False).reset_index(drop=True)
    top_n = 5
    nrows = 1
    ncols = 5
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(ncols * 2.2, nrows * 2))
    for idx in range(top_n):
        ax = axes[idx]
        ax.hist(links[dist_key][idx],
                bins=30,
                alpha=1,
                histtype='stepfilled',  # 'bar', 'barstacked', 'step', 'stepfilled'
                color=color,
                ec='lightgreen',
                rwidth=1,
                )
        ax.set_xlabel('Interaction strength')
        if idx==0:
            ax.set_ylabel('Model count')
        else:
            ax.set_ylabel('')
            ax.set_yticks([])
        ax.set_title(links['Regulator'][idx]+'-'+links['Target'][idx])
        ax.set_xlim([0.1, 0.45])
        if idx == top_n-1:
            handles = []
            for i, color in enumerate([color]):
                handles.append(ax.scatter([], [], marker='o', label=name, s=50,
                                          edgecolor='black', color=color, linewidth=.2))
            ax.legend(handles=handles,
                      bbox_to_anchor=(1, .8), prop={'size': 12}, loc='right',
                      frameon = False
                      )
    return fig
def plot_grn_weights_distributions_per_network(ax, links, name, studies, study_colors) -> None:
    """ Plots the distribution of weights obtained for a network

    """
    serif_font()
    for i, study in enumerate(links):
        ax.hist(study['Weight'], bins=50, alpha=0.6,
                        histtype='bar', #'bar', 'barstacked', 'step', 'stepfilled'
                        color=study_colors[i],
                        rwidth=1.5,
                       )
    handles = []

    for i, color in enumerate(study_colors):
        handles.append(ax.scatter([],[], marker='o', label=studies[i],
         edgecolor='black', color=color, linewidth=.2))

        ax.legend(handles=handles,
            bbox_to_anchor=(1,1), prop={'size': 9}
            )
        ax.set_xlabel('Interaction strength')
        ax.set_ylabel('Density')

        ax.set_ymargin(.2)
        ax.set_title(name)
def format_links_string(links: pd.DataFrame, gene_names: List[str]) -> pd.DataFrame:
    """Converts links to golden links style
    Links are given as gene 1 to gene 2 with a weight.
    Golden link is combination of all genes with 0 or 1 for the true links.
    We assume that any link given in links is 1 disgarding quantity of Weight.
    """
    n_genes = len(gene_names)
    regs = np.repeat(np.asarray(gene_names), n_genes-1)
    targs = []
    for i in range(n_genes):
        targs.extend(gene_names[:i] + gene_names[i+1:])
    golden_links = pd.DataFrame({'Regulator':regs, 'Target': targs})
    golden_links['Weight'] = (golden_links['Regulator'].isin(links['Regulator']) & golden_links['Target'].isin(links['Target'])).astype(int)

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
def normalize_links(links):
    """
        Nornalize the links based on the std
    """
    links_n = links.copy()
    links_n.loc[:,'Weight'] = links['Weight']/np.std(links['Weight'])
    return links_n

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
