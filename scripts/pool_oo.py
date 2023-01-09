"""
Reads the results of hyper param tuning (RF), 
pool them (because there are multiple runs), average them, and output them to files.
"""
from imports import *

def process_oo(study, n=10):
    '''
        Read the results of tunning and pool them. 

    '''
    best_paramss = []
    best_scoress = []
    for i in range(n):
        best_params, best_scores = utils.read_write_oo(study, mode='read', OUTPUT_DIR=OUTPUT_DIR, i=i)
        best_paramss.append(best_params)
        best_scoress.append(best_scores)

    #- process best scores
    bestscores_pool = np.array(best_scoress).T.tolist() # n_genes*n
    bestscores = np.mean(bestscores_pool, axis=1).tolist() # n_genes
    
    #- process best params
    bestparams_pool = [] # n_gene*best_params (n values)
    for prot_i in range(len(best_paramss[0])):
        bestparams = {}
        for key in best_paramss[0][0].keys():
            vector = []
            for item in best_paramss:
                vector.append(item[prot_i][key])
            bestparams[key] = vector
        bestparams_pool.append(bestparams)
    #- average best params for each gene
    bestparams = [{key:np.mean(vector) for key, vector in bestparams.items()} for bestparams in bestparams_pool]
    #- make some param int
    bestparams = [{key:value if key=='decay_coeff' else int(value) for key,value in bestparam.items()} for bestparam in bestparams]
    #- pool for all genes
    # pool_of_pool = {}
    # for key in best_paramss[0][0].keys():
    #     vector_of_vectors = []
    #     for item in bestparams_pool:
    #         vector_of_vectors+=list(item[key])

    #     pool_of_pool[key] = vector_of_vectors
    # bestparams_pool_allgenes = np.array([list(item.values()) for item in bestparams_pool])
    
    with open(os.path.join(OUTPUT_DIR,'calibration', f'oo_{study}.txt'),'w') as f:
        print({'best_params':bestparams, 'best_scores':bestscores},file=f)
    
    with open(os.path.join(OUTPUT_DIR,'calibration', f'oo_pool_{study}.txt'),'w') as f:
        print({'bestparams_pool':bestparams_pool, 'bestscores_pool':bestscores_pool},file=f)
        
process_oo('ctr', n=100)
process_oo('mg', n=100)