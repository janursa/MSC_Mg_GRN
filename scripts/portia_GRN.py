from imports import *

# estimated_decay_coeffs = utils.estimate_decay_rates([data_ctr, data_mg], [time,time])

dir_portia = os.path.join(MAIN_DIR,'..','external/PORTIA-master')
sys.path.insert(0,dir_portia)
import portia as pt
def data_process(data, time):
    dataset = pt.GeneExpressionDataset()
    exp_id = 1
    # decay_coeffs = estimated_decay_coeffs
    for data_i, time_i in zip(data,time):
        dataset.add(pt.Experiment(exp_id, data_i))
        exp_id+=1
    return dataset


def GRN(data, study):
    dataset = data_process(data, time)
    M_bar = pt.run(dataset, method='fast')
    links_df = tools.Links.format(M_bar, protnames)
    
    utils.Links.read_write_links(links=links_df, study=study, mode='write',method='portia',OUTPUT_DIR=OUTPUT_DIR)
    return links_df




if __name__ == '__main__':
    # - read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
    data_mg = utils.process_data(df_target, study='mg', standardize=False)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')

    # - run
    links = GRN(data_ctr, 'ctr')
    GRN(data_mg, 'mg')

    #- compare to string
    top_n = 250

    # links_short = utils.Links.choose_top_count(links, n=top_n)
    links_short = utils.Links.choose_top_quantile(links, quantile=.75)
    match_count = utils.Links.compare_network_string(links_short.copy(), OUTPUT_DIR)
    print('match: ', match_count)
