from imports import *
import scipy
import statistics
#- create random links
def create_random_links(links_assembly, n=1000):
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j]
    random_links = links_assembly[0].copy()
    WeightPool = [random.sample(weights, len(weights)) for i in range(n)]
    weightpoolpool = np.array([i for j in WeightPool for i in j])
    weightpoolvector = []
    for i in range(len(links_assembly[0])):
        weightpoolvector.append(np.random.choice(weightpoolpool, n))
    random_links['WeightPool'] = weightpoolvector
    random_links['Weight']= np.mean(weightpoolvector, axis=1)
    return random_links

def pool_comparision(links_target):
    """
    Extract match counts for different runs of the model, weightpool
    """
    match_count = []
    weightpool = np.array(links_target['WeightPool'].values.tolist()).T
    for weight in weightpool:
        aa = links_target.copy()
        aa['Weight'] = weight
        aa_short = utils.Links.choose_top_count(aa, n=top_n)
        n = utils.Links.compare_network_string(aa_short, OUTPUT_DIR, verbose=False)
        match_count.append(n)
    return match_count



if __name__ == '__main__':
    replica_n = 200
    #- retreive the data
    #-- rf
    links_rf = retreive_grn('ctr','rf')
    #-- ridge
    links_ridge_decay_coeff_estimated = pd.read_csv(os.path.join(OUTPUT_DIR, 'GRN/links_ctr_ridge_estimated.csv'), index_col=False)
    links_ridge = pd.read_csv(os.path.join(OUTPUT_DIR, 'GRN/links_ctr_ridge.csv'), index_col=False)

    #-- portia
    links_portia = pd.read_csv(os.path.join(OUTPUT_DIR, 'GRN/links_ctr_portia.csv'), index_col=False)
    #-- random
    links_all = [links_rf,links_ridge_decay_coeff_estimated,links_portia]
    links_all_n = [utils.Links.nomalize(item) for item in links_all]
    links_random = create_random_links(links_all_n, n=replica_n)

    #- calculate match count. choose top links for different methods
    top_n = 100
    #-- RF: 
    links_rf_short = utils.Links.filter_fitscore(links_rf)
    match_count_rf = pool_comparision(links_rf_short)
    #-- ridge
    def process(links):
        links_short = utils.Links.choose_top_count(links, n=top_n)
        match_count = utils.Links.compare_network_string(links_short.copy(), OUTPUT_DIR)
        match_count = [match_count for i in range(replica_n)]
        return match_count
    match_count_ridge_estimated = process(links_ridge_decay_coeff_estimated)
    match_count_ridge = process(links_ridge)
    #-- portia
    match_count_portia = process(links_portia)
    #-- random
    match_count_random = pool_comparision(links_random)


    #- visualize it

    # datas = [match_count_random_grni, match_count_grni, match_count_random_portia, match_count_portia, match_count_ensemble]
    # labels = ['Random\n geneRNI', 'geneRNI', 'Random\n Portia', 'Portia','Ensemble\n methods']

    datas = [match_count_random, match_count_rf, match_count_ridge_estimated,match_count_ridge, match_count_portia]
    labels = ['Random\n', 'RF', 'Ridge\n e', 'Ridge', 'Portia']
    print('Mean match counts:')
    for label, data in zip(labels,datas):
        print(f'{label}: {np.mean(data)}')



    # pvalues = []
    # s, p = scipy.stats.ttest_ind(match_count_grni, match_count_random_grni)
    # pvalues.append(1) #match_count_random_grni compared to itself
    # pvalues.append(p<0.05)
    # s, p = scipy.stats.ttest_ind(match_count_ensemble, match_count_random_grni)
    # pvalues.append(p<0.05)
    # s, p = scipy.stats.ttest_ind(match_count_portia, match_count_random_portia)
    # pvalues.append(1) #match_count_random_portia compared to itself
    # pvalues.append(p<0.05)

    # value_flags = np.array((np.mean(datas[1:], axis=1) - np.mean(datas[0]))>0)


    # def define_sign(p):
    #     if p:
    #         sign = r'$*$'
    #     else:
    #         sign=''
    #     return sign
    # flags = flags_sigs*value_flags
    # sig_signs = ['',r'$*$','',r'$*$','']
    sig_signs = ['' for i in range(len(datas))]
    #- plot violin for match scores
    importlib.reload(utils)
    fig=utils.Links.plot_match_counts(datas=datas, labels=labels, sig_signs=sig_signs)
    fig.savefig(os.path.join(OUTPUT_DIR, 'postprocess/match_count.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(OUTPUT_DIR, 'postprocess/match_count.pdf'))
