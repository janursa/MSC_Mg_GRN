from imports import *
data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
data_mg = utils.process_data(df_target, study='mg',standardize=False)
print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')

def GRN(data, study='ctr', i=0):
    dataset = Data(gene_names=protnames,
        ss_data=None,
        ts_data=[data],
        time_points=[time])

    #- read the outputs of tunning
    best_params, _ = utils.read_write_oo(study, mode='read', OUTPUT_DIR=OUTPUT_DIR)
    #- run the network inference
    param = dict(estimator_t='RF')

    ests, train_scores, links_df, oob_scores, test_scores = \
        geneRNI.network_inference(dataset, gene_names=protnames, 
                                  param=param, param_unique=best_params, verbose=False)

    utils.Links.read_write_links(links=links_df, study=study, mode='write',i=i,OUTPUT_DIR=OUTPUT_DIR)
def run_batch(data, study, istart, iend):
    for i in range(istart, iend):
        GRN(data, study, i)
        print(f'{i} finished')

# run_batch(data_ctr,'ctr',0, 200)
run_batch(data_mg,'mg',100, 145) 

