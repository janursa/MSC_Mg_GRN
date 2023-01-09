"""
Tune hyper params and decay coeffs using geneRNI
- Only train set was used
"""
from imports import *

data_ctr = utils.process_data(df_target, study='ctr')
data_mg = utils.process_data(df_target, study='mg')

specs = dict(
#     train_flag=True # maximizing training data
)

param = dict(estimator_t = 'RF')

def tune(study):
    """
        Read training data, tune the params, and output best params and best scores.
    """
    if study == 'ctr':
        data = data_ctr
    else:
        data = data_mg

    dataset = Data(
        gene_names=protnames,
        ss_data=None,
        ts_data=[data],
        time_points=[time]
        )

    
    best_scores, best_params, best_ests, sampled_permts_sorted = \
            search_param.rand_search(dataset, param=param, 
                                     param_grid=param_grid, 
                                     n_jobs=10, n_sample=60, **specs)
    print(f'mean best score {np.mean(best_scores)}')
    return best_scores, best_params

def batch_run(study, istart=None, iend=None):
    if istart is None:
        print(f'----Tuning for {study} -----')
        best_scores, best_params = tune(study)
        utils.read_write_oo(study, mode='write', best_params=best_params, best_scores=best_scores, OUTPUT_DIR=OUTPUT_DIR)
    for i in range(istart, iend):
        print(f'----Tuning for {study} iteration {i}-----')
        best_scores, best_params = tune(study)
        utils.read_write_oo(study, mode='write', best_params=best_params, best_scores=best_scores, OUTPUT_DIR=OUTPUT_DIR, i=i)

# batch_run('ctr',50,100)
batch_run('mg',50,100)
