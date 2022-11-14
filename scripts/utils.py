# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:08:06 2022

@author: nourisa
"""
import copy
import pandas as pd
import numpy as np
from sklearn import ensemble
from sklearn import inspection
import matplotlib.pyplot as plt
from scipy import stats

time = [1,2,3,4,7,8,9,10,11,14,21]

def plot_time_series(df, prots, c_tag='ctr_', s_tag='mg_', p_name='Protein', time=time, ee=0.5, **kywrds):
    """ plots ctr and sample in time series indicating the sig margin """
    n_prots = len(prots)
    if n_prots == 1:
        n_cols = 1
        n_rows = 1
        fig = plt.figure(tight_layout=True, figsize=(5*n_prots,3))
        axs = [[fig.add_subplot(1, 1, 1)]]
    else:
        n_cols = min([3,n_prots])
        n_rows = 1+int((n_prots-n_cols)/3) + (1+n_prots-n_cols)%3

        fig, axs = plt.subplots(n_rows, n_cols,  sharey=True, tight_layout=True, figsize=(4*n_cols,2*n_rows+2))
    
    linewidth = 2
    ctrs = [c_tag + str(ii) for ii in time]
    samples = [s_tag + str(ii) for ii in time]
    def plot_single(ax, prot, ctrs_data, samples_data, ee_u, ee_d):
        ax.plot(time, ctrs_data, '-o', linewidth=linewidth, color ='b', label ='Ctr')
        ax.plot(time, samples_data, '-x', linewidth=linewidth, color ='g', label ='Mg')
        ax.plot(time, ee_u, color = 'r',linewidth=linewidth-1, linestyle = '--', label = '0.5 error bounds')
        ax.plot(time, ee_d, linewidth=linewidth-1, linestyle = '--', color = 'r')
        ax.set_xticks(time)
        ax.set_title(prot)
        ax.set_xlabel('Days')
        ax.set_ylabel('Intensity (log2)')

    count = 0

    for i in range(n_rows):
        for j in range(n_cols):
            try:
                prot = prots[count]
            except:
                return
            df_tag = df.loc[df[p_name] == prot]
            try:
                
                ctrs_data = df_tag[ctrs].iloc[0,:].to_list()
            except:
                print(f'check if {prot} exist in df')
                raise ValueError()
            ee_u = [ii + ee for ii in ctrs_data]
            ee_d = [ii - ee for ii in ctrs_data]
            samples_data = df_tag[samples].iloc[0,:].to_list()
            # print(ctrs_data)
            # print(samples_data)
            plot_single(axs[i][j], prot, ctrs_data, samples_data, ee_u, ee_d)
            if count == 0:
                axs[i][j].legend(loc='best', bbox_to_anchor=(1.1, 1.7),  ncol=3)
            count+=1
            if count == len(prots):
                break
            
def plot_time_series_mutual(df1, df2, prots, c_tag='ctr_', s_tag='mg_', p_name='Entry', time=time, ee=0.5, **kywrds):
    """ plots ctr and sample in time series indicating the sig margin """
    n_prots = len(prots)

    n_cols = 2
    n_rows = len(prots)

    fig, axs = plt.subplots(n_rows, n_cols,  sharey=True, tight_layout=True, figsize=(8,2*n_rows+2))
    
    linewidth = 2
    ctrs = [c_tag + str(ii) for ii in time]
    samples = [s_tag + str(ii) for ii in time]
    def plot_single(ax, prot, ctrs_data, samples_data, ee_u, ee_d):
        ax.plot(time, ctrs_data, '-o', linewidth=linewidth, color ='b', label ='Ctr')
        ax.plot(time, samples_data, '-x', linewidth=linewidth, color ='g', label ='Mg')
        ax.plot(time, ee_u, color = 'r',linewidth=linewidth-1, linestyle = '--', label = '0.5 error bounds')
        ax.plot(time, ee_d, linewidth=linewidth-1, linestyle = '--', color = 'r')
        ax.set_xticks(time)
        ax.set_title(prot)
        ax.set_xlabel('Days')
        ax.set_ylabel('Intensity (log2)')
    def func(df_tag, ax): #auxillary function
        try:
            ctrs_data = df_tag[ctrs].iloc[0,:].to_list()
        except:
            print(f'check if {prot} exist in df')
            raise ValueError()
        ee_u = [ii + ee for ii in ctrs_data]
        ee_d = [ii - ee for ii in ctrs_data]
        samples_data = df_tag[samples].iloc[0,:].to_list()

        plot_single(ax, prot, ctrs_data, samples_data, ee_u, ee_d)

    count = 0

    for i in range(n_rows):
        
        try:
            prot = prots[count]
        except:
            return
        
        func(df1.loc[df1[p_name] == prot],axs[i][0])
        func(df2.loc[df2[p_name] == prot],axs[i][1])

        if count == 0:
            axs[i][0].legend(loc='best', bbox_to_anchor=(1.1, 1.7),  ncol=3)
        count+=1
        if count == len(prots):
            break
            
    
    
    

def sig_test(df, day, c_tag, s_tag, ee, **kywrds): 
    """ Significant changes between ctr and sample are kept considering a single time point """
    df_sig = df.copy()
    diff = df[s_tag + str(day)] - df[c_tag + str(day)].values
    indices = diff>ee    
    df_sig = df_sig.loc[indices,:]
    df_sig.reset_index(inplace=True,drop=True)
    return df_sig

def test_overexpression(df,time, c_tag, s_tag, ee ,min_s, **kywrds):
    """ Over-expression between ctr and sample exceeding ee threshold, for min_s occurance"""
    df_sig = df.copy()
    cols_ctr = [x for x in df.columns if c_tag in x]
    cols_mg = [x for x in df.columns if s_tag in x]
    cols_diff = ['d_'+ str(i) for i in time]
    diff = df[cols_mg]-df[cols_ctr].values
    df_sig[cols_diff] = diff
    bools = diff > ee
    indices = [False for i in range(len(df))]
    for ii in range(len(df)):
        try: 
            indices[ii] = bools.iloc[ii,:].value_counts()[True] >= min_s
        except:
            continue  
    df_sig = df_sig.loc[indices,:]
    df_sig.reset_index(inplace=True,drop=True)
    print('Sig proteins: ',len(df_sig))
    return df_sig

def test_sig_gene_ttest(s1, s2):
    aa, pvalue = stats.ttest_ind(s1, s2)
    if pvalue < 0.05:
        return True
    else:
        return False
def test_sig_gene(s1, s2):
    return test_sig_gene_ttest(s1, s2)
    
        
def test_sig(df,time, c_tag, s_tag, **kywrds):
    """ Significant changes between ctr and sample are kept for time series """

    cols_ctr = [x for x in df.columns if c_tag in x]
    cols_mg = [x for x in df.columns if s_tag in x]
    # we take all samples for now. TODO: change it for early and late phases
    df_ctr = df.loc[:,cols_ctr]
    df_mg = df.loc[:,cols_mg]
    indices = []
    for i, row in df.iterrows():
        gene_ctr = row[cols_ctr].to_list()
        gene_mg = row[cols_mg].to_list()
        if test_sig_gene(gene_ctr, gene_mg):
            indices.append(i)

    return df.iloc[indices,:]
    # indices = [False for i in range(len(df))]
    # for ii in range(len(df)):
    #     try: 
    #         indices[ii] = bools.iloc[ii,:].value_counts()[True] >= min_s
    #     except:
    #         continue  
    # df_sig = df_sig.loc[indices,:]
    # df_sig.reset_index(inplace=True,drop=True)
    # print('Sig proteins: ',len(df_sig))

def listwise_deletion(df): 
    """ removes rows with a single zero """
    df_copy = copy.deepcopy(df)
    for col in df_copy.columns:
        drop_line = df_copy.index[df_copy[col].isnull()].tolist()
        df_copy = df_copy.drop(drop_line)
    df_copy.reset_index(inplace=True,drop=True)
    return df_copy
def rename_missing_symbols(df,p_name,**kywrds): 
    """ Name the missing symbols in protein names """
    df_c = copy.deepcopy(df)
    nulls = df_c[p_name].isnull() 
    unique_names = ['p_{}'.format(ii) for ii in range(sum(nulls))]
    map_ = {i:name for i,name in zip([i for i, x in enumerate(nulls) if x],unique_names)}
    df_c.loc[nulls, p_name] = unique_names
    print('Remaining missing names: ',[x for x in df_c[p_name] if not isinstance(x,str)])
    with open('results/missing_names.txt','w') as f:
        for key,value in map_.items():
            f.write('original index: '+str(key)+' -> '+value+'\n')
    return df_c

def tailor_names(o_df, time, p_name, c_func, s_func, o_c_func, o_s_func, **kywrds):
    """ """
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

def plot_hist(xs = [], names = [], **specs):
    if not xs or not names:
        raise ValueError('xs and names are not given')
    def plot_host_single(ax, x, name, **specs):
        ax.hist(x, **specs)
        ax.set_title(name)
        ax.set_xlabel('Value')
        ax.set_ylabel('Dist')
    fig, axs = plt.subplots(1 , len(names) ,tight_layout = True, figsize = (5*len(names),5))
    for (x, name, ax) in zip(xs, names, axs):
        plot_host_single(ax, x, name, **specs)