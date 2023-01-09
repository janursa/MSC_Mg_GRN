from imports import *

data_ctr = utils.process_data(df_target, study='ctr', standardize=True)
data_mg = utils.process_data(df_target, study='mg',standardize=True)
print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
estimated_decay_coeffs = utils.estimate_decay_rates([data_ctr, data_mg], [time,time])

# np.savetxt(os.path.join(OUTPUT_DIR,'data', 'data_ctr.csv'), data_ctr, delimiter=",")
# np.savetxt(os.path.join(OUTPUT_DIR,'data', 'data_mg.csv'), data_mg, delimiter=",")
# np.savetxt(os.path.join(OUTPUT_DIR,'data', 'estimated_decay_coeffs.csv'), estimated_decay_coeffs, delimiter=",")

def GRN(data, study='ctr', i=0):
    dataset = Data(gene_names=protnames,
        ss_data=None,
        ts_data=[data],
        time_points=[time],
        test_size=0
    )
    # run the network inference
    param = dict(estimator_t='ridge', alpha=10)
    # param_unique = [{'alpha':alpha} for alpha in estimated_decay_coeffs]

    ests, train_scores, links_df, oob_scores, test_scores = \
        geneRNI.network_inference(dataset, gene_names=protnames, 
                                  param=param, 
                                  # param_unique=param_unique, 
                                  )
    utils.Links.read_write_links(links=links_df, study=study, mode='write',method='ridge',OUTPUT_DIR=OUTPUT_DIR)


GRN(data_ctr,'ctr')
GRN(data_mg,'mg')