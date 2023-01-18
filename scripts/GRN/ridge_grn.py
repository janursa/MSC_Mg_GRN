import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *
if __name__ == '__main__':
    #- read the data
    data_ctr = utils.process_data(df_target, study='ctr', standardize=True)
    data_mg = utils.process_data(df_target, study='mg', standardize=True)
    print('Data shape:', np.array(data_ctr).shape, '(n_samples_time_series*n_genes)')
    #- params
    method='ridge'
    param = dict(estimator_t=method)
    # estimated_decay_coeffs = utils.estimate_decay_rates([data_ctr, data_mg], [time, time])

    #- network inference
    test_size = 0
    _, param_unique_ctr = utils.calibration.retreive_data(study='ctr', method=method, OUTPUT_DIR=OUTPUT_DIR)
    _, param_unique_mg = utils.calibration.retreive_data(study='mg', method=method, OUTPUT_DIR=OUTPUT_DIR)

    _, trainscores_ctr, links_ctr, _, testscores_ctr = utils.links.grn(data=data_ctr, gene_names=protnames, time_points=time,
                                    test_size=test_size, param=param, param_unique=param_unique_ctr)
    _, trainscores_mg, links_mg, _, testscores_mg = utils.links.grn(data=data_mg, gene_names=protnames, time_points=time,
                                            test_size=test_size, param=param, param_unique=param_unique_mg)
    #- write to files
    utils.links.read_write_links(links=links_ctr, study='ctr', mode='write', method=method, output_dir=OUTPUT_DIR)
    utils.links.read_write_links(links=links_mg, study='mg', mode='write', method=method, output_dir=OUTPUT_DIR)

    output_dir = os.path.join(OUTPUT_DIR, 'GRN', method)
    utils.links.write_scores(method=method, study='ctr', trainscores=trainscores_ctr, testscores=testscores_mg, output_dir=output_dir)
    utils.links.write_scores(method=method, study='mg', trainscores=trainscores_mg, testscores=testscores_mg, output_dir=output_dir)

    #- plot
    fig = utils.calibration.plot_scores(data_ctr=trainscores_ctr, data_sample=trainscores_mg, ylabel='Train score')
    fig.savefig(os.path.join(OUTPUT_DIR, 'GRN', method, 'trainscores.pdf'))
    #
    plt.show()
    #- compare to string
    def compare(links):
        # top_n = 100
        # links_short = utils.links.choose_top_count(links_ctr, n=top_n)
        links_short = utils.links.choose_top_quantile(links, quantile=.75)
        match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
        print('match: ', match_count)
    compare(links_ctr)


