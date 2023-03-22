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

from typing import Tuple, List
import copy
from pathlib import Path

from sklearn import preprocessing
from sklearn import inspection

import matplotlib
import matplotlib.pyplot as plt






sys.path.insert(0, os.path.join(pathlib.Path(__file__).parent.resolve(), '../..'))
from scripts.imports import DATA_DIR



def make_title_pretty(name):
    name = name.replace('day1_11', 'Early')
    name = name.replace('day1_21', 'Late')
    name = name.replace('portia', 'Portia')
    name = name.replace('ridge', 'Ridge')
    name = name.replace('_','-')
    return name
def comic_font():
    matplotlib.rc('font', family='Comic Sans MS')
    matplotlib.rc('text', usetex='false')
    matplotlib.rcParams.update({'font.size': 10})
def serif_font():
    matplotlib.rc('font', family='serif')
    matplotlib.rc('text', usetex='false')
    matplotlib.rcParams.update({'font.size': 10})
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

def read_write_data(tag, mode, data=None):
    file_name = Path(DATA_DIR)/f"data_{tag}.csv"
    assert mode in ['read','write']
    if mode == 'read':
        return np.loadtxt(file_name,delimiter=",")
    else:
        assert data is not None
        np.savetxt(file_name, data, delimiter=",")
def process_data(df_target, study) -> np.array :
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




def estimate_decay_rates(TS_data: Tuple[np.array], time_points:  Tuple[np.array]):
    
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

