def sig_test(df, day, ee = 0.5): # receives the df and the specs, and return a df based on the sig change in proteins at the given time
    diff = df['Mg_'+str(day)]-df['ctr_'+str(day)].values
    df_sig = df.copy()
    df_sig['d_'+str(day)] = diff
    indices = diff>ee    
    df_sig = df_sig.loc[indices,:]
    df_sig.reset_index(inplace=True,drop=True)
    return df_sig
def edit_GO(df,time,**kywrds): 
    """ Keep the GO terms only """
    df_e = copy.deepcopy(df)
    for GO_term in ['Gene ontology (biological process)']:
        new_col = []
        for item in df[GO_term]:
            try :
                item_e1 = item.split(';')
            except:
                new_col.append(item)
                continue
            item_e2 = [x.split('[')[0] for x in item_e1]
            item_e2 = [x.strip() for x in item_e2]
            new_col.append(item_e2)
        df_e[GO_term] = new_col
    return df_e
def _sig_test(df, ee = 0.5, _min = 1, time = []): # receives the df and the specs, and return the sig proteins and their log2fc values
        cols = value_cols(df)
        def ff(day):
            ctr, Mg = df['ctr_{}'.format(day)], df['Mg_{}'.format(day)]
            log2FC = Mg - ctr
            sigs = abs(log2FC) > ee
            return sigs,log2FC
    
    sigs_bool_time = pd.DataFrame()
    sigs = pd.DataFrame()
    for day in time:
        sigs_bool_time[day],sigs[day] = ff(day)  

    indices = []
    for ii in range(len(sigs_bool_time)):
        try:
            count_true = sigs_bool_time.iloc[ii,:].value_counts()[True]
        except:
            continue
        if count_true >= _min:
            indices.append(ii)
    sigs = sigs.iloc[indices,:] # filter
    symbols = df.loc[indices,'Symbols'] # filter
    sigs['Symbols'] = symbols.values
    return symbols,sigs

## dynGENIE3
import importlib
from dynGENIE3 import dynGENIE3 
importlib.reload(dynGENIE3)
import _pickle
f = open('dynGENIE3/TS_data.pkl','rb')
(TS_data, time_points, decay_rates, gene_names) = _pickle.load(f)
            #TS_data: 3*21*10
            #time_points: 3*21
            #decay_rates: 10
            #gene_names: 10
f.close()
# (VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3.dynGENIE3(TS_data,time_points)
# VIM: putative regulatory links -> from i th gene to j th gene
# alphas: decay rates. 

def fill_single_missings(df): # fill single zeros by extrapolation from neighbors
    start_ii = 1
    if 'ID' in df:
        start_ii += 1
    for ii in range(1,len(df)): 
        for jj in range(start_ii,len(df.columns)):
            if df.iloc[ii,jj] == 0:
                if jj == start_ii: # only right hand side
                    rr = df.iloc[ii,jj+1]
                elif jj == len(df.columns)-1: #only left hand side
                    rr = df.iloc[ii,jj-1]
                else:
                    if (df.iloc[ii,jj-1]==0) or (df.iloc[ii,jj+1]==0):
                        rr = 0
                    else:
                        rr = (df.iloc[ii,jj-1]+df.iloc[ii,jj+1])/2
                df.iloc[ii,jj] = rr
    return df
##---extract gene ontology 
day = 1
GO_term = 'GO_BP'
indices = sig_time[day].index
df_sig = pd.concat([sig_time[day],df.loc[indices,[GO_term]]],axis = 1)
df_sig_f = df_sig.drop(df_sig.index[df_sig[GO_term].isnull()])
GO_terms_e = [x for xx in df_sig_f[GO_term].tolist() for x in xx]
GO_terms = list(set(GO_terms_e))
GO_terms_i = {i:GO_terms_e.count(i) for i in GO_terms}
GO_terms_i

df_sig_f.head()

## output files for casual path
df_o = df 
df_o['ID'] = df['Symbols']
df_o = listwise_deletion(df_o)
df_o['Sites'] = np.nan
df_o['Effect'] = np.nan
df_o.to_csv('data/sig_data.txt',sep='\t')
print('Data size, with at least one sig change: {}'.format(len(df_f)))
with open('data/parameters.txt') as ff:
    params = ff.read()
for tt in time:
    params += "control-value-column = ctr_{}\n test-value-column = Mg_{}\n".format(tt,tt)
with open('data/parameters.txt','w') as ff:
    ff.write(params)
## run time series network analysis using dynGENIE3
import importlib
from dynGENIE3 import dynGENIE3 
importlib.reload(dynGENIE3)
df_target = df_interp

params = dict(
    nthreads = 10,
    save_models = False,
    compute_quality_scores = True,
    ntrees = 100,
    max_depth = 20,
    K = 'all',
    tree_method = 'RF',
    regulators = 'all'
)
(VIM, alphas, oob_scores, stability_score, treeEstimators,scores_train) = dynGENIE3.dynGENIE3(TS_data,time_points,alpha='from_data',SS_data=None,gene_names=gene_names,**params)
links = dynGENIE3.get_link_list(VIM,gene_names=gene_names)

#------ interpolation
df_1 = df.interpolate(limit=1) # filling up one missing value
df_1 = listwise_deletion(df_1)
print('Data size, rows with more than 1 zero were removed: {}'.format(len(df_1)))
#------
df_2 = df.interpolate(limit=2) # filling up one missing value
df_2 = listwise_deletion(df_2)
print('Data size, rows with more than 2 zero were removed: {}'.format(len(df_2)))
df_ffill = df.ffill() 
df_ffill = listwise_deletion(df_ffill)
print('Data size, df_ffill: {}'.format(len(df_ffill)))