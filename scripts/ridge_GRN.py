from imports import *
def GRN(data, study='ctr', i=0):
    dataset = Data(
        gene_names=protnames,
        ss_data=None,
        ts_data=[data],
        time_points=[time],
        test_size=0
    )
    # run the network inference
    param = dict(estimator_t='ridge')
    param_unique = [{'decay_coeff':decay_coeff} for decay_coeff in estimated_decay_coeffs]

    ests, train_scores, links_df, oob_scores, test_scores = \
        geneRNI.network_inference(dataset,
                                  gene_names=protnames,
                                  param=param, 
                                  # param_unique=param_unique,
                                  )
    print(test_scores)
    utils.Links.read_write_links(links=links_df, study=study, mode='write',method='ridge',OUTPUT_DIR=OUTPUT_DIR)
    return links_df

if __name__ == '__main__':
    #- read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
    data_mg = utils.process_data(df_target, study='mg', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    estimated_decay_coeffs = utils.estimate_decay_rates([data_ctr, data_mg], [time, time])
    #- grn
    links = GRN(data_ctr, 'ctr')
    GRN(data_mg, 'mg')
    #- compare to string
    top_n = 100
    links_short = utils.Links.choose_top_count(links, n=top_n)
    # links_short = utils.Links.choose_top_quantile(links, quantile=.75)
    match_count = utils.Links.compare_network_string(links_short.copy(), OUTPUT_DIR)
    print('match: ', match_count)


