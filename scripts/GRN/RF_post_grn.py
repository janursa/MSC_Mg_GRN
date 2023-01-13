"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links.
    Reads oob and train scores and plot them
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *

def add_score_to_links(links, scores, protnames):
    links['FitScore'] = None
    for target, score in zip(protnames, scores):
        links.loc[links['Target'] == target, 'FitScore'] = score
    return links

if __name__ == '__main__':
    method = 'RF'
    replica_n = 100
    # - read links extracted during different runs.
    links_ctr = utils.links.pool_links('ctr', protnames, output_dir=OUTPUT_DIR, n=replica_n, method=method)
    links_sample = utils.links.pool_links('mg', protnames, output_dir=OUTPUT_DIR, n=replica_n, method=method)
    # - retreive scores and add them to the links
    _, scores_ctr = utils.calibration.read_write_oo(study='ctr', mode='read', OUTPUT_DIR=OUTPUT_DIR)
    _, scores_sample = utils.calibration.read_write_oo(study='mg', mode='read', OUTPUT_DIR=OUTPUT_DIR)

    links_ctr = add_score_to_links(links_ctr, scores_ctr, protnames)
    links_sample = add_score_to_links(links_sample, scores_sample, protnames)

    #- save pooled links
    links_ctr.to_pickle(os.path.join(OUTPUT_DIR, 'GRN', method,
                                     f'links_ctr.csv'))  # save with pickle because there is a list of items in data
    links_sample.to_pickle(os.path.join(OUTPUT_DIR, 'GRN', method, f'links_mg.csv'))

    #- pool oob and train scores, and plot them
    def pool_scores(method, study, replica_n, OUTPUT_DIR):
        oob_scoress = []
        train_scoress = []
        for i in range(replica_n):
            oob_scores = utils.links.read_write_links(study=study, mode='read_oob_scores', i=i, OUTPUT_DIR=OUTPUT_DIR, method=method)
            train_scores = utils.links.read_write_links(study=study, mode='read_train_scores', i=i, OUTPUT_DIR=OUTPUT_DIR, method=method)

            oob_scoress.append(oob_scores)
            train_scoress.append(train_scores)
        oob_scores_means = np.mean(oob_scoress, axis=1)
        train_scores_means = np.mean(train_scoress, axis=1)
        return oob_scores_means, train_scores_means

    oob_scores_ctr, train_scores_ctr = pool_scores(method, 'ctr', replica_n, OUTPUT_DIR)
    oob_scores_mg, train_scores_mg = pool_scores(method, 'mg', replica_n, OUTPUT_DIR)

    # utils.links.plot_scores(data_ctr=[train_scores_ctr, train_scores_mg], data_sample=[oob_scores_ctr, oob_scores_mg])
    fig = utils.calibration.plot_scores(data_ctr=train_scores_ctr, data_sample=train_scores_mg, ylabel='Training score')
    fig.savefig(os.path.join(OUTPUT_DIR, 'GRN/RF', 'oobscores.pdf'))
    fig = utils.calibration.plot_scores(data_ctr=oob_scores_ctr, data_sample=oob_scores_mg, ylabel='OOB score')
    fig.savefig(os.path.join(OUTPUT_DIR, 'GRN/RF', 'trainscores.pdf'))

    plt.show()
