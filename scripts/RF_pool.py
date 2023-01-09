"""
    Process the results of GRN using RF by pooling them and adding oob scores to the links. 
"""
from imports import *
method = 'RF'
#- read links extracted during different runs. 
links_ctr = utils.Links.pool_links('ctr', protnames, output_dir=OUTPUT_DIR, n=200, method=method)
links_sample = utils.Links.pool_links('mg', protnames, output_dir=OUTPUT_DIR, n=200, method=method)
#- retreive scores
_, scores_ctr = utils.read_write_oo(study='ctr', mode='read', OUTPUT_DIR=OUTPUT_DIR)
_, scores_sample = utils.read_write_oo(study='mg', mode='read', OUTPUT_DIR=OUTPUT_DIR)

def add_score_to_links(links, scores, protnames):
    links['FitScore'] = None
    for target, score in zip(protnames, scores):
        links.loc[links['Target'] == target, 'FitScore'] = score
    return links
links_ctr = add_score_to_links(links_ctr, scores_ctr, protnames)
links_sample = add_score_to_links(links_sample, scores_sample, protnames)
# #- save
links_ctr.to_pickle(os.path.join(OUTPUT_DIR,'GRN', f'links_ctr_{method}.csv')) # save with pickle because there is a list of items in data
links_sample.to_pickle(os.path.join(OUTPUT_DIR,'GRN', f'links_mg_{method}.csv'))