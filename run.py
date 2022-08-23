## grid search
import newEstimator 
import importlib
import numpy as np
importlib.reload(newEstimator)
import sys
import os
import pandas as pd
import numpy as np
import copy
import math
import matplotlib.pyplot as plt
import json
import utils
main_dir = "C:/Users/nourisa/Downloads/testProjs/omics"
sys.path.insert(0,main_dir)
specs = dict(
    o_df_dir = os.path.join(main_dir,'data','original_omics.xlsx'),
    df_dir = os.path.join(main_dir,'data','omics.csv'),
    time = [1,2,3,4,7,8,9,10,11,14,21],
    p_ID = 'Entry',  # The column name in the original data to be used as protein ID
    c_tag = 'ctr_',
    s_tag = 'mg_',
    c_func = lambda t: 'ctr_' + str(t),  # The tag used in the tailored dataframe for control
    s_func = lambda t: 'mg_' + str(t),  # The tag used in the tailored dataframe for sample
    o_c_func = lambda t: 'LogLFQ intensity ' + str(t) + '_0',  # Func to extract the control data from the original database 
    o_s_func = lambda t: 'LogLFQ intensity ' + str(t) + '_1',  # Func to extract the sample data from the original database
    ee = 0.5,  # Threshold for log2FC to retain the data
    min_s = 2,  # Min occurance in time series in order to retain the data  
)
## read the original data
df = pd.read_csv('data/primary_omics.csv')
# extract sig proteins
df_sig = utils.sig_df(df, **specs)
print('Sig proteins: {}'.format(len(df_sig)))
# impute missing values: available techniques to try: interpolate, univariate and multivariate: https://scikit-learn.org/stable/modules/impute.html
df_interp = df_sig.interpolate() 
df_interp = utils.listwise_deletion(df_interp)
print('Sig interpolated proteins: {}'.format(len(df_interp)))

# reformat the data for time series
df_target = df_interp
time = specs['time']
p_ID = specs['p_ID']
data_ctr = np.array(df_target.loc[:,['ctr_'+str(day) for day in time]].values)
data_mg = np.array(df_target.loc[:,['mg_'+str(day) for day in time]].values)
(TS_data, time_points, decay_rates, gene_names) = \
    [data_ctr.transpose(),data_mg.transpose()], [time,time],[],df_target[p_ID].tolist()
print('n_exp n:',len(TS_data))
print('n_time*n_genes',TS_data[0].shape)

# convert data to n_samples*n_genes format
Xs,ys = newEstimator.process_data(TS_data, time_points, regulators = 'all')
print('n_genes: ',len(ys),len(Xs))
print('n_sample*n_genes : ',Xs[0].shape)
print('n_sample for y : ',len(ys[0]))

if __name__ == '__main__':
    param = dict(
        estimator_t = 'RF',
        # estimator_t = 'HGB',
        n_estimators = 100,
        # validation_fraction = None,
        # loss='squared_error'
    )
    param_grid = dict(
    #     max_depth = [10,30],
        n_estimators = np.arange(100, 400, 50),
        alpha = np.arange(0, 1.1, .01),
        # learning_rate = np.arange(0.01, 0.5, .05)
    #     max_leaf_nodes=31, max_depth=None
    #     min_samples_leaf
    #     max_bins = np.arange(50, 255, 50)
    )
    specs = dict(
        n_jobs = 50,
        # cv = 5,
    )
    best_scores, best_params, best_ests = newEstimator.grid_search(Xs, ys, param, param_grid, **specs)
    print('mean:', np.mean(best_scores), ' std:', np.std(best_scores))
    with open('results/grid_search.txt', 'w') as f:
        print(outs, file=f)

    