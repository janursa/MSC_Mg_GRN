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

## build sklearn random forest
n_estimators = 100
criterion='squared_error'
max_depth = None
min_samples_split = 2
min_samples_leaf = 1
min_weight_fraction_leaf = 0
max_features = sqrt(len(feature_names) - 1)/len(feature_names) # according to GENIE3 suggestion
max_leaf_nodes = None
min_impurity_decrease = 0
bootstrap = True
oob_score = True  # to use out-of-bag samples to estimate the generalization score (available for bootstrap=True)
n_jobs = None # number of jobs in parallel. fit, predict, decision_path, and apply can be done in parallel over the trees
random_state=None # controls randomness in bootstrapping as well as drawing features
verbose=0 
warm_start=False # reuse the slution of the previous call to fit and add more ensembles to the estimator. look up on Glosery
ccp_alpha=0 # complexity parameter used for minima cost-complexity pruning. by default, no prunning
max_samples = None # if bootstrap is True, the number of samples to draw from the samples. if none, draw X.shape[0]

sk_rf_reg = RandomForestRegressor(n_estimators=n_estimators, criterion=criterion, max_depth=max_depth, min_samples_split=min_samples_split,
                                  min_samples_leaf=min_samples_leaf, min_weight_fraction_leaf=min_weight_fraction_leaf, 
                                  max_features=max_features, max_leaf_nodes=max_leaf_nodes, min_impurity_decrease=min_impurity_decrease,
                                  bootstrap=bootstrap, oob_score=oob_score, n_jobs=n_jobs, random_state=random_state, verbose=verbose,
                                  warm_start=warm_start, ccp_alpha=ccp_alpha, max_samples=max_samples)


## feature importance
feature_importance = reg.feature_importances_
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + 0.5
fig = plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.barh(pos, feature_importance[sorted_idx], align="center")
plt.yticks(pos, np.array(diabetes.feature_names)[sorted_idx])
plt.title("Feature Importance (MDI)")

result = permutation_importance(
    reg, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)
sorted_idx = result.importances_mean.argsort()
plt.subplot(1, 2, 2)
plt.boxplot(
    result.importances[sorted_idx].T,
    vert=False,
    labels=np.array(diabetes.feature_names)[sorted_idx],
)
plt.title("Permutation Importance (test set)")
fig.tight_layout()
plt.show()

## load the data
import sys
import os
main_dir = "C:/Users/nourisa/Downloads/testProjs/omics"
sys.path.insert(0,main_dir)
import _pickle

data_dir = 'C:/Users/nourisa/Downloads/testProjs/omics/dynGENIE3/dynGENIE3_data/real_networks/data/yeast_data.pkl'
goldern_links_dir = ''
with open(data_dir,'rb') as f:
    (TS_data, time_points, genes, TFs, alphas) = _pickle.load(f)
print('exp n:',len(TS_data))
print('time*genes',TS_data[0].shape)
print('len of alphas:',len(alphas))
print('len of genes:',len(genes))
print('len of TFs:',len(TFs))
# print(genes)
# Xs,ys = newEstimator.process_data(TS_data, time_points)
# print('genes: ',len(Xs))
# print('time across all experiments*genes:',Xs[0].shape)
# check_estimator(TargetEstimator)



from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import numpy as np

# Define the model inputs
problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359]]
}

# Generate samples
param_values = saltelli.sample(problem, 1024)

# Run model (example)
Y = Ishigami.evaluate(param_values)

# Perform analysis
Si = sobol.analyze(problem, Y, print_to_console=True)

# Print the first-order sensitivity indices
print(Si['S1'])