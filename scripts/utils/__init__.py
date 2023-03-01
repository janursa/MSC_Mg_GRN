# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:08:06 2022

@author: nourisa
"""
import os
import pathlib
import sys

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

import typing
import copy

from sklearn import preprocessing
from sklearn import inspection

import matplotlib
import matplotlib.pyplot as plt

# MAIN_DIR = os.path.join(pathlib.Path(__file__).parent.resolve(), '../..')
# sys.path.insert(0, MAIN_DIR)
#
# -- import from geneRNI
# geneRNI_dir = os.path.join(MAIN_DIR,'..','geneRNI')
# sys.path.insert(0, geneRNI_dir)
# from geneRNI import geneRNI, search_param
# from geneRNI.data import Data
# from geneRNI.models import get_estimator_wrapper
# from . import links
# from utils import links
# from utils import calibration
# from utils import enrichplot
# from utils import sensitivity_analysis
# from utils import VSA
#

def MG_noise_F(links, n_relica=100, std=.3) -> typing.Tuple[pd.DataFrame]:
    """ Multiplicitive noise
         Creates n_relica noised links
    """
    noised_link = links.copy()
    original_weights = links['Weight'].to_numpy(float)
    noised_links = [noised_link.assign(Weight= original_weights*np.random.normal(loc=1, scale=std, size=len(links)))
                                       for i in range(n_relica)]
    return noised_links
def AG_noise_F(links, n_relica=100, rel_std=.1)-> typing.Tuple[pd.DataFrame]:
    """ Additive noise
             Creates n_relica noised links
    """
    noised_link = links.copy()
    original_weights = links['Weight'].to_numpy(float)
    std = rel_std * np.std(original_weights)
    noised_links = [noised_link.assign(Weight= original_weights+np.random.normal(loc=0, scale=std, size=len(links)))
                                       for i in range(n_relica)]
    return noised_links
def make_title_pretty(name):
    # parts = name.split('_')
    # parts[0] = parts[0][0:1].upper()+parts[0][1:] + ' phase'# upper case for first letter
    # if 'Combined' in parts[0]:
    #     parts[0] = 'Combined'
    # parts[1] = '('+parts[1]+' proteins)'
    # return '\n'.join(parts)
    name = name.replace('day1_11', 'Early')
    name = name.replace('day1_21', 'Late')
    name = name.replace('portia', 'Portia')
    name = name.replace('ridge', 'Ridge')
    name = name.replace('_','-')

    # parts = name.split('-')
    # parts[-1] = f'({ parts[-1]})'
    # name = ' '.join(parts)
    return name
def comic_font():
    matplotlib.rc('font', family='Comic Sans MS')
    matplotlib.rc('text', usetex='false')
    matplotlib.rcParams.update({'font.size': 10})
def serif_font():
    matplotlib.rc('font', family='serif')
    matplotlib.rc('text', usetex='false')
    matplotlib.rcParams.update({'font.size': 10})

# def plot_time_series(df, prots, c_tag='ctr_', s_tag='mg_', p_name='Protein', time=time, ee=0.5, **kywrds):
#     """ plots ctr and sample in time series indicating the sig margin """
#     n_prots = len(prots)
#     if n_prots == 1:
#         n_cols = 1
#         n_rows = 1
#         fig = plt.figure(tight_layout=True, figsize=(5*n_prots,3))
#         axs = [[fig.add_subplot(1, 1, 1)]]
#     else:
#         n_cols = min([3,n_prots])
#         n_rows = 1+int((n_prots-n_cols)/3) + (1+n_prots-n_cols)%3
#
#         fig, axs = plt.subplots(n_rows, n_cols,  sharey=True, tight_layout=True, figsize=(4*n_cols,2*n_rows+2))
#
#     linewidth = 2
#     ctrs = [c_tag + str(ii) for ii in time]
#     samples = [s_tag + str(ii) for ii in time]
#     def plot_single(ax, prot, ctrs_data, samples_data, ee_u, ee_d):
#         ax.plot(time, ctrs_data, '-o', linewidth=linewidth, color ='b', label ='Ctr')
#         ax.plot(time, samples_data, '-x', linewidth=linewidth, color ='g', label ='Mg')
#         ax.plot(time, ee_u, color = 'r',linewidth=linewidth-1, linestyle = '--', label = '0.5 error bounds')
#         ax.plot(time, ee_d, linewidth=linewidth-1, linestyle = '--', color = 'r')
#         ax.set_xticks(time)
#         ax.set_title(prot)
#         ax.set_xlabel('Days')
#         ax.set_ylabel('Intensity (log2)')
#
#     count = 0
#
#     for i in range(n_rows):
#         for j in range(n_cols):
#             try:
#                 prot = prots[count]
#             except:
#                 return
#             df_tag = df.loc[df[p_name] == prot]
#             try:
#
#                 ctrs_data = df_tag[ctrs].iloc[0,:].to_list()
#             except:
#                 print(f'check if {prot} exist in df')
#                 raise ValueError()
#             ee_u = [ii + ee for ii in ctrs_data]
#             ee_d = [ii - ee for ii in ctrs_data]
#             samples_data = df_tag[samples].iloc[0,:].to_list()
#             # print(ctrs_data)
#             # print(samples_data)
#             plot_single(axs[i][j], prot, ctrs_data, samples_data, ee_u, ee_d)
#             if count == 0:
#                 axs[i][j].legend(loc='best', bbox_to_anchor=(1.1, 1.7),  ncol=3)
#             count+=1
#             if count == len(prots):
#                 break
            
# def plot_time_series_mutual(df1, df2, prots, c_tag='ctr_', s_tag='mg_', p_name='Entry', time=time, ee=0.5, **kywrds):
#     """ plots ctr and sample in time series indicating the sig margin """
#     n_prots = len(prots)
#
#     n_cols = 2
#     n_rows = len(prots)
#
#     fig, axs = plt.subplots(n_rows, n_cols,  sharey=True, tight_layout=True, figsize=(8,2*n_rows+2))
#
#     linewidth = 2
#     ctrs = [c_tag + str(ii) for ii in time]
#     samples = [s_tag + str(ii) for ii in time]
#     def plot_single(ax, prot, ctrs_data, samples_data, ee_u, ee_d):
#         ax.plot(time, ctrs_data, '-o', linewidth=linewidth, color ='b', label ='Ctr')
#         ax.plot(time, samples_data, '-x', linewidth=linewidth, color ='g', label ='Mg')
#         ax.plot(time, ee_u, color = 'r',linewidth=linewidth-1, linestyle = '--', label = '0.5 error bounds')
#         ax.plot(time, ee_d, linewidth=linewidth-1, linestyle = '--', color = 'r')
#         ax.set_xticks(time)
#         ax.set_title(prot)
#         ax.set_xlabel('Days')
#         ax.set_ylabel('Intensity (log2)')
#     def func(df_tag, ax): #auxillary function
#         try:
#             ctrs_data = df_tag[ctrs].iloc[0,:].to_list()
#         except:
#             print(f'check if {prot} exist in df')
#             raise ValueError()
#         ee_u = [ii + ee for ii in ctrs_data]
#         ee_d = [ii - ee for ii in ctrs_data]
#         samples_data = df_tag[samples].iloc[0,:].to_list()
#
#         plot_single(ax, prot, ctrs_data, samples_data, ee_u, ee_d)
#
#     count = 0
#
#     for i in range(n_rows):
#
#         try:
#             prot = prots[count]
#         except:
#             return
#
#         func(df1.loc[df1[p_name] == prot],axs[i][0])
#         func(df2.loc[df2[p_name] == prot],axs[i][1])
#
#         if count == 0:
#             axs[i][0].legend(loc='best', bbox_to_anchor=(1.1, 1.7),  ncol=3)
#         count+=1
#         if count == len(prots):
#             break

def listwise_deletion(df): 
    """ removes rows with a single zero """
    df_copy = copy.deepcopy(df)
    for col in df_copy.columns:
        drop_line = df_copy.index[df_copy[col].isnull()].tolist()
        df_copy = df_copy.drop(drop_line)
    df_copy.reset_index(inplace=True,drop=True)
    return df_copy
def rename_missing_symbols(df,p_name='Entry', OUTPUT_DIR=''):
    """ Name the missing symbols in protein names """
    df_c = copy.deepcopy(df)
    nulls = df_c[p_name].isnull()
    unique_names = ['p_{}'.format(ii) for ii in range(sum(nulls))]
    map_ = {i:name for i,name in zip([i for i, x in enumerate(nulls) if x],unique_names)}
    df_c.loc[nulls, p_name] = unique_names
    print('Remaining missing names: ',[x for x in df_c[p_name] if not isinstance(x,str)])
    with open(os.path.join(OUTPUT_DIR,'missing_names.txt'),'w') as f:
        for key,value in map_.items():
            f.write('original index: '+str(key)+' -> '+value+'\n')
    return df_c

def tailor_names(o_df, time, p_name, c_func, s_func, o_c_func, o_s_func, **kywrds):
    """ changes names from the original df to make the reading easier """
    df = pd.DataFrame()
    df[p_name] = o_df[p_name]
    df['Gene'] = o_df['Gene names  (primary)']
    missing_names = df['Gene'].isnull()
    df['Gene'].loc[missing_names] = df[p_name].loc[missing_names].values
    # Rename the data columns for ctr and sample 
    for i in time: 
        df[c_func(i)] = o_df[o_c_func(i)]
    for i in time:
        df[s_func(i)] = o_df[o_s_func(i)]
    print('Data size, original: {}'.format(len(df)))
    return df
    

def remove_zeros(df,c_tag, s_tag, time, **kywrds):
    """ Drop those rows with all zeros """
    df.replace(0,np.nan,inplace=True) 
    cols = [c_tag + str(i) for i in time] + [s_tag + str(i) for i in time]
    df = df.loc[~(df[cols].isnull()).all(axis=1)] 
    df.reset_index(inplace = True)
    df.drop(['index'],axis=1,inplace = True)
    print('Data size, rows with all zeros were removed: {}'.format(len(df)))
    return df

def feature_importance(reg, X_test, y_test,feature_names, params):
    result = inspection.permutation_importance(
        reg, X_test, y_test, **params
    )
    print('*** importance specs****')
    for i,feature_name in enumerate(feature_names):
        print('%s: mean: %.3f std: %.3f'%(feature_name,result['importances_mean'][i],result['importances_std'][i]))
    sorted_idx = result.importances_mean.argsort()
    fig,ax = plt.subplots(1,1,tight_layout=True)
    ax.boxplot(
        result.importances[sorted_idx].T,
        vert=False,
        labels=np.array(feature_names)[sorted_idx],
    )
    ax.set_title("Permutation Importance (test set)")
    return result,plt

def plot_hist(xs, names, **specs):
    '''
        Plots histograpm (distribution) for each x in a seperate window
    '''
    def plot_host_single(ax, x, name, **specs):
        ax.hist(x, **specs)
        ax.set_title(name)
        ax.set_xlabel('Value')
        ax.set_ylabel('Dist')
    fig, axs = plt.subplots(1 , len(names) ,tight_layout = True, figsize = (5*len(names),5))
    for (x, name, ax) in zip(xs, names, axs):
        plot_host_single(ax, x, name, **specs)


def create_check_dir(master_dir, name):
    DIR = os.path.join(master_dir, name)
    if not os.path.isdir(DIR):
        os.makedirs(DIR)
    else:
        if os.listdir(DIR):
            print(f'{name} directory is not empty')
    return DIR


def process_data(df_target, study, standardize=False) -> np.array :
    '''
        Extract training data from df and returns it in a from of array
    '''
    assert(study in ['ctr', 'mg', 'all-in'])
    # extract ctr and mg data
    if study == 'ctr' or study == 'mg':
        # df = df_target.loc[:, [study+'_'+str(day) for day in time_points]].T
        data = df_target.filter(like=study).values.T
        return data

    elif study == 'all-in':
        raise ValueError('not defined')
        # df_ctr = df_target.loc[:, ['ctr_'+str(day) for day in time_points]].T
        # df_mg = df_target.loc[:, ['mg_'+str(day) for day in time_points]].T
        # # add mg as a regulatory factor with 0 for ctr and 1 for mg
        # df_ctr['mg'] = np.zeros(len(time_points))
        # df_mg['mg'] = np.ones(len(time_points))
        #
        # df = pd.concat([df_ctr, df_mg], axis=0)


    if standardize:
        df.iloc[:, :] = preprocessing.scale(df.iloc[:, :])

    return df.values


def estimate_decay_rates(TS_data: typing.Tuple[np.array], time_points:  typing.Tuple[np.array]):
    
    """
    this function is exactly taken from dynGENIE3 code.

    For each gene, the degradation rate is estimated by assuming that the gene expression x(t) follows:
    x(t) =  A exp(-decay_coeff * t) + C_min,
    between the highest and lowest expression values.
    C_min is set to the minimum expression value over all genes and all samples.
    """
    
    ngenes = TS_data[0].shape[1]
    nexp = len(TS_data)
    
    C_min = TS_data[0].min()
    if nexp > 1:
        for current_timeseries in TS_data[1:]:
            C_min = min(C_min,current_timeseries.min())
    
    decay_coeffs = np.zeros((nexp,ngenes))
    
    for (i,current_timeseries) in enumerate(TS_data):
        current_time_points = time_points[i]
        
        for j in range(ngenes):
            
            idx_min = np.argmin(current_timeseries[:,j])
            idx_max = np.argmax(current_timeseries[:,j])
            
            xmin = current_timeseries[idx_min,j]
            xmax = current_timeseries[idx_max,j]
            
            tmin = current_time_points[idx_min]
            tmax = current_time_points[idx_max]
            
            xmin = max(xmin-C_min,1e-6)
            xmax = max(xmax-C_min,1e-6)
                
            xmin = np.log(xmin)
            xmax = np.log(xmax)
            
            decay_coeffs[i,j] = (xmax - xmin) / abs(tmin - tmax)
                
    decay_coeffs = decay_coeffs.max(axis=0)

    return decay_coeffs

