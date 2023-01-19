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


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from imports import *
from utils import create_check_dir, serif_font

from geneRNI.geneRNI import network_inference
from geneRNI.data import Data

def batch_GRN(study, method, i_start, i_end, output_dir, **specs):
    """
        Runs network inference for a number of times
    """
    #- create/check dirs for method, study, (links and scores)
    DIR_METHOD = create_check_dir(output_dir, method)
    DIR_STUDY = create_check_dir(DIR_METHOD, study)
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

def retreive_scores(study, method, output_dir):
    """
        Retreiev train and test scores
    """
    testscores = np.genfromtxt(os.path.join(output_dir, 'GRN', method, f'testscores_{study}.csv'))
    trainscores = np.genfromtxt(os.path.join(output_dir, 'GRN', method, f'trainscores_{study}.csv'))

    return trainscores, testscores
def grn(data, time_points, gene_names, **specs):
    dataset = Data(ts_data=[data], ss_data=None, time_points=[time_points], gene_names=gene_names)
    return network_inference(dataset, gene_names=gene_names, **specs)

def compare_network_string(links, OUTPUT_DIR, verbose=True) -> int:
    '''
        Compare extracted links by GRN to those suggested by vs_string. Returns number of match links.
    '''

    STR_LINKS_FILE = os.path.join(OUTPUT_DIR,'enrichment_analysis','string_interactions.tsv')
    links_string = pd.read_csv(STR_LINKS_FILE,sep='\t',index_col=False)
    links_string = links_string.rename(columns={'#node1':'Regulator','node2':'Target','combined_score':'Weight'})
    links_string = links_string.loc[:,['Regulator','Target','Weight']]  
    if verbose:
        print(f'Number of vs_string links: {len(links_string)}')
        print(f'Number of extracted links: {len(links)}')
    if len(links_string)>len(links):
        raise ValueError('Extracted links cannot be lesser than golden links')

    with open(os.path.join(OUTPUT_DIR, 'postprocess/map_genename_protname.json')) as f:
        map_genename_protname = json.load(f)['map']
    #- convert genenames to protnames in vs_string
    links_string.loc[:,'Regulator'] = [map_genename_protname[name] for name in links_string['Regulator']]
    links_string.loc[:,'Target'] = [map_genename_protname[name] for name in links_string['Target']]
    #- find vs_string links in the extracted links, label them and add stringweight
    links['inString'] = False
    for reg, target, weight in zip(links_string['Regulator'].values, links_string['Target'].values, links_string['Weight'].values):
        links.loc[(links['Regulator']==reg) & (links['Target']==target),'inString'] = True 
        links.loc[(links['Regulator']==reg) & (links['Target']==target),'StringWeight'] = weight
    n= links['inString'].values.tolist().count(True)
    
    return n
def read_write_links(method, study, mode:str, links:pd.DataFrame=None, output_dir='') -> pd.DataFrame:
    '''
        Read write links extracted from GRN 
    '''
    assert(study in ['ctr','mg','combined','random'])
    assert(mode in ['read', 'write'])
    #- determine file location
    DIR = os.path.join(output_dir, 'GRN', method)
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
def pool_GRN_oo(method, study, replica_n, output_dir):
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

def write_scores(method, study, trainscores, testscores , output_dir):
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
def filter_fitscore(links) -> pd.DataFrame:
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
    cut_off = np.quantile(links['Weight'].values.tolist(),q=quantile)
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
def plot_mean_weights(links_s, labels, colors):
    serif_font()
    nrows = 1
    ncols = 3
    fig, axes = plt.subplots(nrows, ncols, tight_layout=True, figsize=(ncols*3, nrows*2.5))
    for idx in range(len(labels)):
        # i = int(idx/(nrows-1))
        # j = idx%ncols
        # ax = axes[i][j]

        ax = axes[idx]
        for i,study in enumerate(links_s[idx]):
            ax.hist(study['Weight'], bins=100, alpha=0.5,
                            histtype='bar', #'bar', 'barstacked', 'step', 'stepfilled'
                            color=colors[i],
                            ec='black',
                            rwidth=1.1,
                            # density = True
                           )
        handles = []
        tags = ['ctr','mg']
        for i, color in enumerate(colors):
            handles.append(ax.scatter([],[],marker='o', label=tags[i],
             edgecolor='black', color=color, linewidth=.2))
        ll = plt.legend(handles=handles, 
            bbox_to_anchor=(1,1), prop={'size': 9}
            # title='Enriched Term'
            )
        # ax.add_artist(ll)
        ax.set_xlabel('Normalized interaction strength')
        ax.set_ylabel('Density')
        ax.set_ymargin(.2)
#         ax.set_xmargin(.1)
        # ax.set_xlim([-.5,8])
        ax.set_title(labels[idx])
def plot_match_counts(datas, labels, sig_signs):
    matplotlib.rcParams.update({'font.size': 12})

    fig, axes = plt.subplots(1, 1, tight_layout=True, figsize=(4.7,3.5), 
        )

    ax = axes
    bplot = ax.violinplot(datas, showmeans=True, showextrema=False)

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
    ys = np.max(datas, axis=1) 
    for i, sign in enumerate(sig_signs):
        if sign != '':
            ax.annotate(sign, xy=(xs[i],ys[i]),
                ha='center', 
                va='bottom',
                )
    return fig


def plot_match_counts_series(match_counts_list, links_names, top_quantile_list):
    """
    Plots match counts for different top selected links as a line plot. Each line is for different methods such as RF and Ridge
    """
    top_quantile_list = np.round(top_quantile_list, 2)
    match_counts_sum = np.sum(match_counts_list, axis=1).astype(int) #final scores
    # add scores to the names
    links_names = [f'{name}: {score}' for name, score in zip(links_names,match_counts_sum)]
    serif_font()
    fig, ax = plt.subplots(1, 1, tight_layout=True, figsize=(4.7, 3.5))
    # x_range = [min(xticks_labels), max(xticks_labels)]
    x = np.linspace(min(top_quantile_list), max(top_quantile_list), len(top_quantile_list))
    colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey']
    linestyles = ['-', '--', '-.', ':']
    assert (len(match_counts_list) == len(links_names))
    for i, data in enumerate(match_counts_list):
        ax.plot(top_quantile_list, data, label=links_names[i], color=colors[i], alpha=1, linewidth=2,
                linestyle=linestyles[i],
                marker='o')
        # ax.fill_between(x, data, alpha=.1)

    ax.set_ylabel('Number of matched interactions')
    ax.set_xlabel('Selected top quantile')
    xticks = np.arange(min(top_quantile_list), max(top_quantile_list), .05)
    ax.set_xticks(xticks)
    ax.set_xmargin(.15)
    ax.set_ymargin(.15)

    plt.legend(frameon=False)
    # - annotate score on the right side of each line
    
    # for i, (score, y) in enumerate(zip(match_counts_sum, match_counts_list)):
    #     ax.annotate(score, xy=(top_quantile_list[-1] + .01, y[-1] - 0.5),
    #                 ha='center',
    #                 va='bottom',
    #                 size=10,
    #                 alpha=.5
    #                 )

    return fig
def create_random_links(links_assembly, n=1000):
    #- TODO: create matrix (using protname)
    links_assembly = [normalize_links(links) for links in links_assembly]
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j] #flatten
    random_links = links_assembly[0].copy()
    weightpoolvector = []
    for i in range(len(links_assembly[0])):
        weightpoolvector.append(random.sample(weights, n))
    random_links['WeightPool'] = weightpoolvector
    random_links['Weight']= np.mean(weightpoolvector, axis=1)
    return random_links
def normalize_links(links):
    """
        Nornalize the links based on the std
    """
    links_n = links.copy()
    links_n.loc[:,'Weight'] = links['Weight']/np.std(links['Weight'])
    return links_n

def convert_links_to_nodes_edges(links, protnames, scores):
    '''
        Takes links as pd and seperates nodes and edges
    '''
    edges = links.copy()
    edges.rename(columns={'Regulator':'source', 'Target':'target'}, inplace=True)
    sum_ws = [sum(links.loc[links['Regulator']==gene,:]['Weight']) for gene in protnames]
    nodes = pd.DataFrame(data={'Entry':protnames,'SumWeight':sum_ws, 'FitScore': scores})
    return nodes, edges

def visualize_network(nodes, edges, study, protnames, preferred_names=None, OUTPUT_DIR=''):
    import igraph as ig
    DIR = os.path.join(OUTPUT_DIR, 'GRN')
    shapes = ['circle', 'circle','circle', 'rectangle']
    #- replace the names with indices: index of the genes in the protnames
    edges_index = copy.deepcopy(edges)
    edges_index.loc[:,'source'] = [protnames.index(name) for name in list(edges['source'].values)]
    edges_index.loc[:,'target'] = [protnames.index(name) for name in list(edges['target'].values)]
    edges_list = edges_index.loc[:,['source','target']].values
    #- Construct a graph 
    n_vertices = len(protnames)
    g = ig.Graph(n_vertices, edges=edges_list, directed=False)
    g["title"] = "Regulatory network"
    #- gene names instead of protnames
    if preferred_names is None:
        g.vs["name"] =  protnames
    else:
        g.vs["name"] = preferred_names
    #- create the layout
    coords = []
    for i in range(len(g.vs)):
        coords.append([i,i])
    layout_obj = ig.Layout(coords)
    #- plot
    weights = edges['Weight'].values
    vertex_sizes = nodes['FitScore']
    grad_colors = ig.GradientPalette("red", "green", len(protnames))
    fig, ax = plt.subplots(figsize=(10,10))
    vertex_size_correction = .8/np.max(vertex_sizes)
    edge_size_correction = 1/np.max(weights)
    # print(f'mean : {np.mean(weights)}, max: {np.max(weights)}, min: {np.min(weights)}')
    ig.plot(
        g,
        target=ax,
        layout="kk",
    #     layout= layout_obj,
        vertex_size=[x*vertex_size_correction for x in vertex_sizes],
        vertex_color=[grad_colors.get(i) for i in range(len(protnames))],
    #     vertex_frame_width=8,
    #     vertex_frame_color="white",
        vertex_label=g.vs["name"],
        vertex_label_size=12,
        edge_width=[x*edge_size_correction for x in weights],
        edge_arrow_width = .5,
        vertex_shape = [shapes[i] for i in nodes['Role'].values]
        # vertex_shape = ig.drawing.shapes.DiamondDrawer()
    #     edge_color=["#7142cf" if x>np.mean(g.es['weights']) else "#AAA" for married in g.es["married"]]
    
    )
    plt.show()
    fig.savefig(os.path.join(DIR,f'GRN_{study}.png'), dpi=300, transparent=True)


    
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