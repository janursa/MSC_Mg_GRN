from imports import *
#- create random links
def create_random_links(links_assembly, n=1000):
    links_assembly = [utils.links.nomalize(links) for links in links_assembly]
    weights = [links['Weight'].values.tolist() for links in links_assembly]
    weights = [i for j in weights for i in j]
    random_links = links_assembly[0].copy()
    weightpoolvector = []
    for i in range(len(links_assembly[0])):
        weightpoolvector.append(random.sample(weights, n))
    random_links['WeightPool'] = weightpoolvector
    random_links['Weight']= np.mean(weightpoolvector, axis=1)
    return random_links

def batch_comparision(links_target):
    """
    Extract match counts for different runs of the model, weightpool
    """
    match_count = []
    weightpool = np.array(links_target['WeightPool'].values.tolist()).T
    for weight in weightpool:
        aa = links_target.copy()
        aa['Weight'] = weight
        aa_short = utils.links.choose_top_count(aa, n=top_n)
        n = utils.links.compare_network_string(aa_short, OUTPUT_DIR, verbose=False)
        match_count.append(n)
    return match_count

if __name__ == '__main__':
    replica_n = 100
    #- retreive the data
    #-- rf
    links_rf = pd.read_pickle(os.path.join(OUTPUT_DIR, 'GRN', 'RF',
                                     f'links_ctr.csv'))  # save with pickle because there is a list of items in data
    #-- ridge
    links_ridge = utils.links.read_write_links('ctr',method='ridge',mode='read_links',OUTPUT_DIR=OUTPUT_DIR)

    #-- portia
    links_portia =utils.links.read_write_links('ctr',method='portia',mode='read_links',OUTPUT_DIR=OUTPUT_DIR)
    #-- random
    links_all = [links_rf,links_ridge,links_portia]
    links_all_n = [utils.links.nomalize(item) for item in links_all]
    links_random = create_random_links(links_all_n, n=replica_n)

    #- calculate match count. choose top links for different methods
    top_n = 100
    # -- random
    match_count_random = batch_comparision(links_random)
    #-- RF: 
    # links_rf_short = utils.links.filter_fitscore(links_rf)
    match_count_rf = batch_comparision(links_rf)
    #-- ridge
    def process(links):
        links_short = utils.links.choose_top_count(links, n=top_n)
        match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR)
        match_count = [match_count for i in range(replica_n)]
        return match_count
    match_count_ridge = process(links_ridge)
    #-- portia
    match_count_portia = process(links_portia)

    #- visualize it
    # datas = [match_count_random_grni, match_count_grni, match_count_random_portia, match_count_portia, match_count_ensemble]
    # labels = ['Random\n geneRNI', 'geneRNI', 'Random\n Portia', 'Portia','Ensemble\n methods']

    datas = [match_count_random, match_count_rf, match_count_ridge, match_count_portia]
    labels = ['Random\n', 'RF', 'Ridge', 'Portia']
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
    fig=utils.links.plot_match_counts(datas=datas, labels=labels, sig_signs=sig_signs)
    fig.savefig(os.path.join(OUTPUT_DIR, f'postprocess/match_count_{top_n}.png'), dpi=300, transparent=True)
    fig.savefig(os.path.join(OUTPUT_DIR, f'postprocess/match_count_{top_n}.pdf'))
    plt.show()