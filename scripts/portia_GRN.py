from imports import *
data_ctr = utils.process_data(df_target, study='ctr', standardize=False)
data_mg = utils.process_data(df_target, study='mg',standardize=False)
print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
estimated_decay_coeffs = utils.estimate_decay_rates([data_ctr, data_mg], [time,time])

dir_portia = os.path.join(MAIN_DIR,'..','external/PORTIA-master')
sys.path.insert(0,dir_portia)
import portia as pt
def data_process(data, time):
    dataset = pt.GeneExpressionDataset()
    exp_id = 1
    alphas = estimated_decay_coeffs
    for data_i, time_i in zip(data,time):
        print(time_i)
        dataset.add(pt.Experiment(exp_id, data_i))
        exp_id+=1
    return dataset


def GRN(data, study):
    dataset = data_process(data, time)
    M_bar = pt.run(dataset, method='fast')
    links_df = tools.Links.format(M_bar, protnames)
    
    utils.Links.read_write_links(links=links_df, study=study, mode='write',method='portia',OUTPUT_DIR=OUTPUT_DIR)

#- read the data
data_ctr = np.genfromtxt(os.path.join(OUTPUT_DIR,'data', 'data_ctr.csv'),  delimiter=',')
data_mg = np.genfromtxt(os.path.join(OUTPUT_DIR,'data', 'data_mg.csv'),  delimiter=',')

#- run 
GRN(data_ctr, 'ctr')
GRN(data_mg, 'mg')