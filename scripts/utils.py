# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:08:06 2022

@author: nourisa
"""
import os
import copy
import pandas as pd
import numpy as np
from sklearn import ensemble
from sklearn import inspection
import matplotlib.pyplot as plt
from scipy import stats
import typing
import random
from sklearn import preprocessing
import json
import matplotlib
import shutil 
import matplotlib.font_manager as font_manager
import matplotlib.markers as mmarkers
import scipy
import statistics
# shutil.rmtree(matplotlib.get_cachedir())

time = [1,2,3,4,7,8,9,10,11,14,21]


def comic_font():
    matplotlib.rc('font', family='Comic Sans MS') 
    matplotlib.rc('text', usetex='false') 
    matplotlib.rcParams.update({'font.size': 10})
def serif_font():
    matplotlib.rc('font', family='serif') 
    matplotlib.rc('text', usetex='false') 
    matplotlib.rcParams.update({'font.size': 10})


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
    with open('results/data/missing_names.txt','w') as f:
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

def read_write_oo(study='ctr', mode='read', bestparams=None, bestscores=None, i=None, OUTPUT_DIR=''):
    '''
        Read and writes calibration results (best params and best scores) to file
    '''
    assert(mode in ['read','write'])
    assert(study in ['ctr','mg','combined'])
    DIR = os.path.join(OUTPUT_DIR,'calibration')
    if i is None:
        FILE = os.path.join(DIR, f'oo_{study}.txt')
    else:
        FILE = os.path.join(DIR, 'pool', f'oo_{study}_{i}.txt')
    if mode == 'write':
        with open(FILE,'w') as f:
            print({'best_params':bestparams, 'best_scores':bestscores},file=f)
    elif mode == 'read':
        with open(FILE,'r') as f:
            oo = eval(f.read())
        return oo['best_params'], oo['best_scores']
def process_data(df_target, study='ctr', standardize=False) -> np.array :
    '''
        Extract training data from df and returns it in a from of array
    '''
    assert(study in ['ctr','mg','combined'])
    # extract ctr and mg data
    if study=='ctr' or study=='mg':
        df = df_target.loc[:,[study+'_'+str(day) for day in time]].T
    elif study=='combined':
        df_ctr = df_target.loc[:,['ctr_'+str(day) for day in time]].T
        df_mg = df_target.loc[:,['mg_'+str(day) for day in time]].T
        # add mg as a regulatory factor with 0 for ctr and 1 for mg
        df_ctr['ctr_mg'] = np.zeros(len(time)) 
        df_mg['mg_mg'] = np.ones(len(time)) 
        df = None #TODO: how to concatenate
        protnames.append('mg') #TODO: needs evaluation

    if standardize:
        df.iloc[:,:] = preprocessing.scale(df.iloc[:,:])

    return df.values
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


class VSA:
    '''
        Scatter plot for Vester's sensitivity analysis
    '''
    linewidth = 2
    markers = ['o','o','o','*']
    colors = ['lightgreen','lightblue','yellow','r']
    sizes = [40*i for i in [1,1.5,2,5]]
    roles = ['Buffering','Passive','Active','Critical']
    @staticmethod
    def thresholds(AS, PS) -> typing.Tuple[float,float]:
        '''
            Calculate the threshold to mark fold change and pvalue 
        '''
        # def pool(df_1, df_2, key):
        #     return list(df_1[key])+list(df_2[key])+
        # PS = pool(df_ctr, df_sample, 'PS')
        # AS = pool(df_ctr, df_sample, 'AS')
        # PS_t, AS_t = np.mean(PS), np.mean(AS)
        PS_t, AS_t = np.quantile(PS,.75), np.quantile(AS,.75)
        # AS_t = 0
        # PS_t = 0
        return AS_t, PS_t

    @staticmethod
    def analyse(links, protnames) -> pd.DataFrame:
        '''
            Vester's sensitivity analysis. 
                - AS: active sum. Sum along rows of the influence matrix and it indicates how much does a variable influence all the others.
                - PS: passive sum. Its is the sum along columns of the influence matrix and it indicates how sensitive a variable is, how does it react to the influence of others
                - Q: AS/PS -> how dominant
                - P: AS.PS -> how participative a variable is
                - Active: +Q
                - Passive: -Q, -P
                - Critical: +Q, +P
                - Buffering: -Q, -P
        ------------------------------------------------------
        inputs: links (DataFrame)
        output: VSA (DataFrame) -> protname: AS, PS, Role
        '''
        

        AS = [sum(links.loc[links['Regulator']==gene,:]['Weight']) for gene in protnames]
        PS = [sum(links.loc[links['Target']==gene,:]['Weight']) for gene in protnames] 
        #- normalize links
        AS, PS = AS/np.std(AS), PS/np.std(PS)
        vsa_results = pd.DataFrame(data={'Entry':protnames,'AS':AS, 'PS':PS})
        #- define the roles: ['Buffering', 'Passive', 'Active', 'Critical'] -> [0, 1, 2, 3] 
        AS_t, PS_t = VSA.thresholds(vsa_results['AS'].values, vsa_results['PS'].values)

        vsa_results.loc[:,'AS'] = vsa_results['AS'] - AS_t 
        vsa_results.loc[:,'PS'] = vsa_results['PS'] - PS_t 
        AS_flags = vsa_results['AS'].values>0
        PS_flags = vsa_results['PS'].values>0
        roles = [2*AS_flag+PS_flag for AS_flag, PS_flag in zip(AS_flags,PS_flags)]
        vsa_results.loc[:,'Role'] = roles
        return vsa_results
    @staticmethod
    def mark_y(ax, y_t):
        """
        Mark threshold line for the sig
        """
        line_color = 'grey'
        dx = ax.get_xlim()
        ax.axhline(y_t,color='grey', linestyle='-.')
    @staticmethod
    def mark_x(ax, x_t):
        """
        Mark threshold line for fold change
        """
        line_color = 'grey'
        ax.axvline(x_t,color=line_color, linestyle='--',linewidth=VSA.linewidth)
    @staticmethod
    def scatter(ax,PS, AS, roles):
        
        sc = ax.scatter(PS,AS,
                   color=[VSA.colors[i] for i in roles],
#                    alpha=[1 if flag else .5 for flag in flags]
#                    alpha = .7,
                   linewidths=.1,
                   edgecolors='black',
                   s=50
                  # s = [VSA.sizes[i] for i in roles],
                  )
        # markers = [VSA.markers[i] for i in roles]
        # paths = []
        # for marker in markers:
        #     marker_obj = mmarkers.MarkerStyle(marker)
        #     path = marker_obj.get_path().transformed(
        #                 marker_obj.get_transform())
        #     paths.append(path)
        # sc.set_paths(paths)
    @staticmethod
    def plot_roles(ax):
        '''
            Plots the roles on the corners of the scatter plot
        '''
        # fontsize = 9

        ax.text(.01, .99, 'Active', ha='left', va='top', transform=ax.transAxes, fontweight='bold')
        ax.text(.01, .01, 'Buffering', ha='left', va='bottom', transform=ax.transAxes, fontweight='bold')
        ax.text(.99, .01, 'Passive', ha='right', va='bottom', transform=ax.transAxes, fontweight='bold')
        ax.text(.99, .99, 'Critical', ha='right', va='top', transform=ax.transAxes, fontweight='bold')
        
    @staticmethod
    def plot_prot_names(ax, AS, PS, gene_names, roles):
        '''
            Prints names of critical genes on the scatter plot
        '''
        xlim = ax.get_xlim()[1]-ax.get_xlim()[0]
        ylim = ax.get_ylim()[1]-ax.get_ylim()[0]
        
        arrow_specs = {'arrow_type':'arc3', 'rad':.2}
        for i, gene_name in enumerate(gene_names):
            offset_x = 0.2
            offset_y = 0.2
            rad = .2
            role = roles[i]
            if role != 3: #just critical 
                continue
            x = PS[i]
            y = AS[i]
            if gene_name == 'Q02218' or gene_name=='P26447':
                offset_y = -.05
                offset_x = .3
                rad=.2
            if gene_name == 'P62854':
                offset_y = -.25
                offset_x = .15
                rad=.2
            if gene_name == 'Q07065':
                offset_y = .2
                offset_x = -.25
                rad=-.2
            if gene_name == 'Q9UJZ1':
                offset_y = .15
                offset_x = 0.1
                rad=-.2
            if gene_name == 'P21281':
                offset_y = -.05
                offset_x = 0.22
                rad=.2
            # if gene_name == 'P21281', 'P26447', 'P62854','Q02218','Q07065','Q9UJZ1'


            ax.annotate(gene_name,xy=(x,y), xytext=(x+offset_x*xlim,y+offset_y*ylim), fontsize=7,
             arrowprops={'arrowstyle':'->','lw':0.9, 'connectionstyle':f'arc3, rad={rad}','color':'black', 'alpha':.7}
             ,horizontalalignment='center')
            # ax.annotate(gene_name, (x+offset_x,y+offset_y), fontsize=fontsize)
        

        
        # fontsize = 7
        
        

    @staticmethod
    def postprocess(ax, title, show_axis_names=True):
        if show_axis_names:
            ax.set_xlabel('PS')
            ax.set_ylabel('AS')

        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        xlen, ylen = xlim[1]-xlim[0], ylim[1]-ylim[0]
        ax.set_title(title)
        ax.margins(.3)
        ax.set_xticks([])
        ax.set_yticks([])
    @staticmethod
    def role_change(df_ctr, df_sample, target_roles=[3]):    
        '''
            Finds those genes with a change in the role from ctr to sample.

        Inputs: target_role -> one of 0,1,2,3
        '''
        flags = []
        for target_role in target_roles:
            flags.append(list(df_ctr['Role'] == target_role))
            flags.append(list(df_sample['Role'] == target_role))
        flags = np.array(flags).T  
        flags = [np.any(row) for row in flags]
        
        df_ctr_c = df_ctr.loc[flags,:].reset_index(drop=True)
        df_sample_c = df_sample.loc[flags,:].reset_index(drop=True)

        df_role_changed = pd.DataFrame()
        df_role_changed['Entry'] = df_ctr_c['Entry']
        df_role_changed[['AS_1','PS_1','Role_1']]=df_ctr_c[['AS','PS','Role']]
        df_role_changed[['AS_2','PS_2','Role_2']]=df_sample_c[['AS','PS','Role']]
        return df_role_changed
    def plot_ctr_vs_sample(df_ctr, df_sample, preferred_names, OUTPUT_DIR=''):
        '''
            Plots AS/PS for ctr and mg conditions in a 1*2 subplot
        '''
        comic_font()
        DIR = os.path.join(OUTPUT_DIR, 'VSA')
        rows = 1
        cols = 2
        fig, axes = plt.subplots(rows, cols, tight_layout=True, figsize=(cols*3.5, rows*3))
        dfs = [df_ctr, df_sample]
        titles = ['Ctr','Mg']
        for j in range(cols):
            df = dfs[j]
            if df is None:
                continue
            ax = axes[j]
            VSA.scatter(ax, PS=df['PS'], AS=df['AS'], roles=df['Role'])
            # AS_t, PS_t = VSA.thresholds(df['AS'], df['PS'])
            VSA.mark_x(ax, 0)
            VSA.mark_y(ax, 0)
            VSA.plot_prot_names(ax=ax, AS=df['AS'], PS=df['PS'], gene_names=preferred_names, roles=df['Role'])
            VSA.postprocess(ax, title=titles[j])
            VSA.plot_roles(ax)

        handles = []
        for i, color in enumerate(VSA.colors):
            handles.append(ax.scatter([],[],marker='o', label=VSA.roles[i],
             edgecolor='black', color=color, linewidth=.2))
        ll = plt.legend(handles=handles, 
            bbox_to_anchor=(1,1), prop={'size': 8}
            # title='Enriched Term'
            )
        # ax.add_artist(ll)

        fig.savefig(os.path.join(DIR, 'roles_ctr_mg.pdf'),bbox_extra_artists=(ll,), bbox_inches='tight')
        fig.savefig(os.path.join(DIR, 'roles_ctr_mg.png'), dpi=300, transparent=True, bbox_extra_artists=(ll,), bbox_inches='tight')
    def plot_SA(VSAs, name, OUTPUT_DIR):
        '''
            Plots AS/PS for ctr and mg conditions, for n different runs of sensitivity analysis, one window for each prot
        '''
        comic_font() 
        matplotlib.rcParams.update({'font.size': 8})
        arrow_specs = {'arrow_type':'arc3', 'rad':.2}
        ncols = 3
        nrows = 2
        fig, axes = plt.subplots(ncols=ncols, nrows=nrows, layout='tight', figsize=(1.6 * ncols, 1.6 * nrows))
        
        colors = ['lightgreen','purple']
        for idx,(prot, data) in enumerate(VSAs.items()):
            i = int(idx/(ncols))
            j = idx%ncols
            # print(i,j)
            ax = axes[i][j]

            # ax = axes[idx]
            #- plot role change 
            
            for idxx, study in enumerate(data.keys()):
                sc = ax.scatter(data[study]['PS'], data[study]['AS'],
                   color=colors[idxx],
                   alpha = .7,
                   linewidths=.2,
                   edgecolors='black',
                   s = 25
                  # s = [VSA.sizes[i] for i in roles],
                  )
                if np.mean(data[study]['PS'])>0 and np.mean(data[study]['AS'])>0:
                    x = np.mean(data[study]['PS'])
                    y = np.mean(data[study]['AS'])+.8

                    ax.annotate(r'$*$', xy=(x,y), fontsize = 12, color='red')

            VSA.postprocess(ax, title=prot, show_axis_names=False)
            VSA.mark_x(ax, 0)
            VSA.mark_y(ax, 0)
            VSA.plot_roles(ax)
            ax.set_xlim([-3,3])
            ax.set_ylim([-3,3])
            #- arrow
            fontsize = 8
            x0, y0, x1, y1  = np.mean(data['ctr']['PS']), np.mean(data['ctr']['AS']), np.mean(data['mg']['PS']), np.mean(data['mg']['AS'])
            arrow_t = arrow_specs['arrow_type']
            rad = arrow_specs['rad']
            ax.annotate('',xy=(x1,y1), xytext=(x0,y0),
             arrowprops={'arrowstyle':'->','lw':1.5, 'connectionstyle':f'{arrow_t}, rad={rad}','color':'black', 'alpha':.7}
             ,horizontalalignment='center')

        # tags = ['Ctr','Mg']
        # handles = []
        # for i, color in enumerate(colors):
        #     handles.append(ax.scatter([],[],marker='o', label=tags[i], edgecolor='black', color=color, linewidth=.2))
        # ll = plt.legend(handles=handles, 
        #     bbox_to_anchor=(1.3,3.4), 
        #     )
        # fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.png'), dpi=300, transparent=True, bbox_extra_artists=(ll,), bbox_inches='tight')
        # fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.pdf'),bbox_extra_artists=(ll,), bbox_inches='tight')
        fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.pdf'), bbox_inches='tight')
        fig.savefig(os.path.join(OUTPUT_DIR, f'VSA/{name}.png'), bbox_inches='tight',dpi=300)

    def plot_role_change(df, OUTPUT_DIR=''):
        '''
            Plots AS/PS for role change betwen ctr and sample
        '''
        DIR = os.path.join(OUTPUT_DIR, 'VSA')
        #- define default specs for arrows
        arrow_specs = {}
        for prot in df['Entry'].values:
            arrow_specs[prot] = {'arrow_type':'arc3', 'rad':.3}
        #- manual changes to arrow specs
        # arrow_specs['P46379']['rad'] = -.1
        # gene names shown on ctr and not mg
        # an arrow from ctr to mg only for any change from or to critical role
        # names of ctr or mg on the symbols
        rows = 1
        cols = 1
        fig, axes = plt.subplots(rows, cols, tight_layout=True, figsize=(cols*4, rows*4))
        ax = axes 
        #- first ctr
        VSA.scatter(ax, PS=df['PS_1'], AS=df['AS_1'], roles=df['Role_1'])
        #- second sample
        VSA.scatter(ax, PS=df['PS_2'], AS=df['AS_2'], roles=df['Role_2'])
        #- arrow from ctr to sample
        fontsize = 8
        for prot, x0, y0, x1, y1 in zip(df['Entry'].values, df['PS_1'].values, df['AS_1'].values, df['PS_2'].values, df['AS_2'].values): 
            print(prot)
            x, y = x0-.15, y0-.012
            # if prot == 'Q07954' or prot=='Q14444':
            #     x, y = x0-.14, y0+.01
            # if prot == 'P12956':
            #     x, y = x1-.1, y1-.012
            ax.annotate(f'{prot}', xy=(x,y), fontsize = 7)
            arrow_t = arrow_specs[prot]['arrow_type']
            rad = arrow_specs[prot]['rad']
            ax.annotate('',xy=(x1,y1), xytext=(x0,y0),
             arrowprops={'arrowstyle':'->', 'connectionstyle':f'{arrow_t}, rad={rad}'}
             ,horizontalalignment='center')
        
        VSA.mark_x(ax, 0)
        VSA.mark_y(ax, 0)
        VSA.plot_roles(ax)
        VSA.postprocess(ax, title='')

        fig.savefig(os.path.join(DIR, 'roles_ctr_to_mg.png'), dpi=300, transparent=True)
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


    
def convert_links_to_nodes_edges(links, protnames, scores):
    '''
        Takes links as pd and seperates nodes and edges
    '''
    edges = links.copy()
    edges.rename(columns={'Regulator':'source', 'Target':'target'}, inplace=True)
    sum_ws = [sum(links.loc[links['Regulator']==gene,:]['Weight']) for gene in protnames]
    nodes = pd.DataFrame(data={'Entry':protnames,'SumWeight':sum_ws, 'FitScore': scores})
    return nodes, edges


class EnrichPlot:
    
    def mscatter(x,y,ax=None, m=None, **kw):
        """ a custom plot built on scatter plot to enable multi marker visualization
        """
        import matplotlib.markers as mmarkers
        if not ax: ax=plt.gca()
        sc = ax.scatter(x,y,**kw)
        if (m is not None) and (len(m)==len(x)):
            paths = []
            for marker in m:
                if isinstance(marker, mmarkers.MarkerStyle):
                    marker_obj = marker
                else:
                    marker_obj = mmarkers.MarkerStyle(marker)
                path = marker_obj.get_path().transformed(
                            marker_obj.get_transform())
                paths.append(path)
            sc.set_paths(paths)
        return sc
    def plot(datas, tags, size_tag, color_tag, xlabel, marker_types, figsize,legend_color=True, 
            legend_size=True, legend_marker=True, title=''):
        #-----------define prop---------------
        comic_font()

        scale_scatter_size = 60
        
        #-----------calculations---------------
        x = [j for sub in [data[xlabel].values.tolist() for data in datas] for j in sub]
        terms =[j for sub in [data['Description'].values.tolist() for data in datas] for j in sub]

        y = list(range(len(terms), 0, -1)) # numbers for each term
        sizes = np.array([j for sub in [data[size_tag].values.tolist() for data in datas] for j in sub])
        colors = [j for sub in [data[color_tag].values.tolist() for data in datas] for j in sub]
        counts = [len(data['Description'].values.tolist()) for data in datas]
        markers = [marker_types[i] for sub in [np.full(shape=n, fill_value=i, dtype=int) for i,n in enumerate(counts)] for i in sub]

        #-----------plot-----------------------
        fig, ax = plt.subplots(1,1, figsize=figsize, tight_layout=True)
        sc = EnrichPlot.mscatter(x, y, c=colors, s=scale_scatter_size*sizes, 
            m=markers, ax=ax, 
            cmap='autumn'
            # cmap='Spectral'
            )
        ax.set_yticks(y)
        ax.set_yticklabels(terms)
        ax.set_xticks([min(x), int((max(x)+min(x))/2), max(x)])
        ax.set_xlabel('Protein count')
        ax.set_xmargin(0.2)
        ax.set_title(title)
        #- marker legend
        if legend_marker:
            handles = []
            for i, marker in enumerate(marker_types):
                handles.append(ax.scatter([],[],marker=marker, label=tags[i], color='black'))
            l2 = plt.legend(handles=handles, 
                bbox_to_anchor=(2.6,1), 
                # title='Enriched Term'
                )
            ax.add_artist(l2)
        #- size legend 
        if legend_size:
            handles = []
            labels = []
            n_classes = 4
            sizes_classes = np.arange(min(sizes), max(sizes), (max(sizes)-min(sizes))/n_classes)
            for i, size in enumerate(sizes_classes):
                handles.append(ax.scatter([],[], marker='o', label=round(size,1),
                    color='black', 
                    s= size*scale_scatter_size,
                    alpha=1
                    ))
            l1 = plt.legend(bbox_to_anchor=(1, 1), handles=handles, 
                            title=size_tag, fancybox=False, frameon=False)

            ax.add_artist(l1)
        #- color legend
        if legend_color:
            PCM=ax.get_children()[0]
            CB = plt.colorbar(PCM, ax=ax, 
                        aspect=3, 
                        location='right',
                        anchor=(1, 0), 
                        # extend='both', 
                        ticks=[min(colors), max(colors)], 
                        # format="%4.2e",
                        shrink=1
                        )
            CB.ax.set_title(color_tag)
            # CB.ax.tick_params(labelsize=fontsize['fontsize'])
        return fig
class Links:
    def compare_network_string(links, OUTPUT_DIR, verbose=True) -> int:
        '''
            Compare extracted links by GRN to those suggested by string. Returns number of match links.
        '''

        STR_LINKS_FILE = os.path.join(OUTPUT_DIR,'enrichment_analysis','string_interactions.tsv')
        links_string = pd.read_csv(STR_LINKS_FILE,sep='\t',index_col=False)
        links_string = links_string.rename(columns={'#node1':'Regulator','node2':'Target','combined_score':'Weight'})
        links_string = links_string.loc[:,['Regulator','Target','Weight']]  
        if verbose:
            print(f'Number of string links: {len(links_string)}')
            print(f'Number of extracted links: {len(links)}')
        if len(links_string)>len(links):
            raise ValueError('Extracted links cannot be lesser than golden links')

        with open(os.path.join(OUTPUT_DIR, 'postprocess/map_genename_protname.json')) as f:
            map_genename_protname = json.load(f)['map']
        #- convert genenames to protnames in string 
        links_string.loc[:,'Regulator'] = [map_genename_protname[name] for name in links_string['Regulator']]
        links_string.loc[:,'Target'] = [map_genename_protname[name] for name in links_string['Target']]
        #- find string links in the extracted links, label them and add stringweight
        links['inString'] = False
        for reg, target, weight in zip(links_string['Regulator'].values, links_string['Target'].values, links_string['Weight'].values):
            links.loc[(links['Regulator']==reg) & (links['Target']==target),'inString'] = True 
            links.loc[(links['Regulator']==reg) & (links['Target']==target),'StringWeight'] = weight
        n= links['inString'].values.tolist().count(True)
        
        return n
    def read_write_links(links=None, study='ctr', mode='read', i=None, OUTPUT_DIR='', method=None) -> pd.DataFrame or None:
        '''
            Read write links extracted from GRN 
        '''
        assert(study in ['ctr','mg','combined'])
        assert(mode in ['read','write'])
        if i is None:
            FILE=os.path.join(OUTPUT_DIR, f'GRN/links_{study}_{method}.csv')
        else:
            FILE=os.path.join(OUTPUT_DIR, f'GRN/pool/links_{study}_{i}.csv')
        
        if mode=='read':
            links_df = pd.read_csv(FILE, index_col=False)
            return links_df
        elif mode=='write':
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
        
        links_list = [Links.read_write_links(study=study, mode='read', i=i, OUTPUT_DIR=output_dir, method=method) for i in range(n)]

        #- pool the weights 
        ws_pool = np.array([ll['Weight'].values for ll in links_list]).T
        #- create average links df
        ws_mean = np.mean(ws_pool, axis=1)
        links = pd.DataFrame()
        links.loc[:,['Regulator','Target']] = links_list[0].loc[:,['Regulator','Target']]
        links['Weight'] = ws_mean
        links['WeightPool'] = list(ws_pool)
        return links
    
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
        cut_off = np.quantile(links['Weight'],q=quantile)
        return links.loc[links['Weight']>cut_off,:].reset_index(drop=True)

    def choose_top_count(links: pd.DataFrame, n=100)->pd.DataFrame:
        '''
            Filters the given df based on the top count
        ''' 
        links.reset_index(inplace=True, drop=True)
        links.sort_values('Weight',ascending=True,inplace=True)
        return links.iloc[:n,:].reset_index(drop=True)

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
            # gridspec_kw={'width_ratios': [2, 2]}
            )

        ax = axes
        # bplot = ax.boxplot(datas, notch=False, widths =.2, patch_artist=False, 
        #                 meanline=False, showfliers=False)
        bplot = ax.violinplot(datas, showmeans=True, showextrema=False)

        ax.set_ylabel('Number of matched interactions')
        # ax.set_title('(A)')
        ax.set_xticks(list(range(1,len(labels)+1)))
        ax.set_xticklabels(labels,rotation=0)
        ax.set_ymargin(.25)
        # print(bplot)
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
            # ax.annotate(sign, (xs[i],ys[i]), fontsize=14, color='red')
            if sign != '':
                x,y = (xs[i-1]+xs[i])/2,ys[i-1]+.05*max(ys)
                # print(i,x,y)
                # print(ys[i]-1)
                # x,y=1,5

                ax.annotate(sign, xy=(x,y), xytext=(x,y+.04*max(ys)), 
                    ha='center', 
                    va='bottom',
                    arrowprops=dict(arrowstyle='-[, widthB=2.3, lengthB=0.2', lw=1.2))
        return fig
class SensitivityAnalysis:
    def multiGaussNoise(links_ctr, links_sample, protnames, target_prots, OUTPUT_DIR):
        """
            Multi step function to conduct multiplicative gaussian noise analysis
        """
        #- create noised dfs
        links_ctr_noised = SensitivityAnalysis.MG_noise_F(links_ctr)
        links_sample_noised = SensitivityAnalysis.MG_noise_F(links_sample)
        #- run VSA for each df
        series_1 = SensitivityAnalysis.VSA_F(links_ctr_noised, protnames, target_prots)
        series_2 = SensitivityAnalysis.VSA_F(links_sample_noised, protnames, target_prots)
        #- filter series for target proteins and sort them as a dict
        rr_sorted = SensitivityAnalysis.filter_and_sort(series_1, series_2, target_prots)
        #- plot
        VSA.plot_SA(rr_sorted, name='multiGaussNoise', OUTPUT_DIR=OUTPUT_DIR)
        
    def addGaussNoise(links_ctr, links_sample, protnames, target_prots, OUTPUT_DIR):
        """
            Multi step function to conduct additive gaussian noise analysis
        """
        #- create noised dfs
        links_ctr_noised = SensitivityAnalysis.AG_noise_F(links_ctr)
        links_sample_noised = SensitivityAnalysis.AG_noise_F(links_sample)
        #- run VSA for each df
        series_1 = SensitivityAnalysis.VSA_F(links_ctr_noised, protnames, target_prots)
        series_2 = SensitivityAnalysis.VSA_F(links_sample_noised, protnames, target_prots)
        #- filter series for target proteins and sort them as a dict
        rr_sorted = SensitivityAnalysis.filter_and_sort(series_1, series_2, target_prots)
        #- plot
        VSA.plot_SA(rr_sorted, name='addGaussNoise', OUTPUT_DIR=OUTPUT_DIR)

    def VSA_F(links_series, protnames, target_prots):
        """ 
            Runs VSA for the entire noised dfs 
        """
        rr_series = []
        for i in range(len(links_series)):
            links = links_series[i]
            rr = VSA.analyse(links, protnames)
            rr_t = rr.loc[rr['Entry'].isin(target_prots),:]
            rr_series.append(rr)
        return rr_series
    def filter_and_sort(series_1, series_2, target_prots):
        """
            Filters series based on target proteins and puts the data in a form of dict
        """
        sorted_VSA = {}
        extract = lambda VSA_s, prot, tag='AS': [VSA.loc[VSA['Entry']==prot,:][tag].values.tolist()[0] for VSA in VSA_s]
        for prot in target_prots:
            ASs_ctr = extract(series_1, prot=prot, tag='AS')
            PSs_ctr = extract(series_1, prot=prot, tag='PS')
            Roles_ctr = extract(series_1, prot=prot, tag='Role')
            
            ASs_mg = extract(series_2, prot=prot, tag='AS')
            PSs_mg = extract(series_2, prot=prot, tag='PS')
            Roles_mg = extract(series_2, prot=prot, tag='Role')
            
            sorted_VSA[prot] = {'ctr':{'AS':ASs_ctr,'PS':PSs_ctr, 'Role':Roles_ctr},'mg':{'AS':ASs_mg,'PS':PSs_mg, 'Role':Roles_mg}}
        return sorted_VSA
    def MG_noise_F(links):
        n = 100 # noise replica 
        std = 0.15 # standard deviation
        noised_link = links.copy()
        noised_links = [noised_link.assign(Weight= noised_link['Weight']*np.random.normal(loc=1, scale=std, size=len(links)))
                                           for i in range(n)]
        return noised_links
    def AG_noise_F(links):
        n = 100 # noise replica 
        std = 0.5 # standard deviation
        noised_link = links.copy()
        noised_links = [noised_link.assign(Weight= noised_link['Weight']+np.random.normal(loc=1, scale=std, size=len(links)))
                                           for i in range(n)]
        return noised_links

def estimate_decay_rates(TS_data, time_points):
    
    """
    this function is exactly taken from dynGENIE3 code.

    For each gene, the degradation rate is estimated by assuming that the gene expression x(t) follows:
    x(t) =  A exp(-alpha * t) + C_min,
    between the highest and lowest expression values.
    C_min is set to the minimum expression value over all genes and all samples.
    """
    
    ngenes = TS_data[0].shape[1]
    nexp = len(TS_data)
    
    C_min = TS_data[0].min()
    if nexp > 1:
        for current_timeseries in TS_data[1:]:
            C_min = min(C_min,current_timeseries.min())
    
    alphas = np.zeros((nexp,ngenes))
    
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
            
            alphas[i,j] = (xmax - xmin) / abs(tmin - tmax)
                
    alphas = alphas.max(axis=0)


 
    return alphas