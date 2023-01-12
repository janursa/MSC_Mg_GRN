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

    random_links.to_pickle(os.path.join(OUTPUT_DIR, 'postprocess',
                                        f'links_random.csv'))
    return random_links
def compare(links, top_quantile):
    """
    Compare the given links to string
    """
    links_short = utils.links.choose_top_quantile(links,quantile=top_quantile)
    match_count = utils.links.compare_network_string(links_short.copy(), OUTPUT_DIR, verbose=False)
    return match_count
def batch_compare(links, top_quantile):
    """
    Compare the given links to string for each weight set in weightpool
    """
    match_counts = []
    weightpool = np.array(links['WeightPool'].values.tolist()).T
    for weight in weightpool:
        links['Weight'] = weight
        match_count = compare(links, top_quantile)
        match_counts.append(match_count)
    return match_counts
def determine_sig_signes(datas):
    """
    Conducts t test to determine whether datas[1:] are significantly different than datas[0], which is ctr
    Datas: Tuple(DataFrame), e.g. [ctr, RF, Ridge, Portia]
    """
    ctr = datas[0] #random
    #- determine p values: compared to ctr
    pvalues = np.array([])
    for data in datas[1:]:
        s, p = scipy.stats.ttest_ind(data, ctr)
        pvalues = np.append(pvalues, p)
    #- determine whether mean distribution is higher than ctr: we only plot higher ones
    increase_flags = np.array((np.mean(datas[1:], axis=1) - np.mean(ctr))>0)
    #- use p values with value flags
    def define_sign(p):
        if p:
            sign = r'$*$'
        else:
            sign=''
        return sign
    flags = (pvalues<0.05)*increase_flags
    sig_signs = ['']+[define_sign(flag) for flag in flags]
    return sig_signs
def retrieve_data():
    """
    Read links
    """
    # -- rf
    links_rf = pd.read_pickle(os.path.join(OUTPUT_DIR, 'GRN', 'RF',
                                           f'links_ctr.csv'))  # save with pickle because there is a list of items in data
    # -- ridge
    links_ridge = utils.links.read_write_links('ctr', method='ridge', mode='read_links', OUTPUT_DIR=OUTPUT_DIR)
    # -- portia
    links_portia = utils.links.read_write_links('ctr', method='portia', mode='read_links', OUTPUT_DIR=OUTPUT_DIR)

    return links_rf, links_ridge, links_portia
def main(links_random, links_rf, links_ridge, links_portia, top_quantile=0.75, pool_flag=False, plot_flag=False):
    """
    Main function to run the comparision for all datasets
    """
    # -- ridge
    match_count_ridge = compare(links_ridge, top_quantile)
    # -- portia
    match_count_portia = compare(links_portia, top_quantile)
    if pool_flag:
        match_count_random = batch_compare(links_random.copy(), top_quantile)
        match_count_rf = batch_compare(links_rf.copy(), top_quantile)

        match_count_ridge = [match_count_ridge for i in range(len(match_count_random))]
        match_count_portia = [match_count_portia for i in range(len(match_count_random))]
    else:
        match_count_random = compare(links_random.copy(), top_quantile)
        # match_count_random = int(np.mean(match_count_random))
        match_count_rf = compare(links_rf.copy(), top_quantile)
        # match_count_rf = int(np.mean(match_count_rf))

        # match_count_rf = np.mean(match_count_random)
    match_counts = match_count_random, match_count_rf, match_count_ridge, match_count_portia
    if plot_flag:
        labels = ['Random\n', 'RF', 'Ridge', 'Portia']
        for label, data in zip(labels, match_counts):
            print(f'{label}: {np.mean(data)}')
        sig_signs = determine_sig_signes(match_counts)
        fig = utils.links.plot_match_counts(datas=match_counts, labels=labels, sig_signs=sig_signs)

        fig.savefig(os.path.join(OUTPUT_DIR, f'postprocess/match_count_{top_quantile}.png'), dpi=300, transparent=True)
        fig.savefig(os.path.join(OUTPUT_DIR, f'postprocess/match_count_{top_quantile}.pdf'))
    return match_counts

if __name__ == '__main__':

    links_rf, links_ridge, links_portia = retrieve_data()
    # links_random = create_random_links([links_rf, links_ridge, links_portia], n=10) #only once
    links_random = pd.read_pickle(os.path.join(OUTPUT_DIR, 'postprocess', f'links_random.csv'))
    def section1():
        """
        comapre for top quartile 0.9 and 0.75 with distribution plot
        """
        top_quantile = 0.9
        main(links_random, links_rf, links_ridge, links_portia, top_quantile, pool_flag=True, plot_flag=True)

        top_quantile = 0.75
        main(links_random, links_rf, links_ridge, links_portia, top_quantile, pool_flag=True, plot_flag=True)
        plt.show()
    # section1()
    def section2():
        """
        Obtains top matches for a range of top quantile values. Only mean values are considered
        in match calculation
        """
        top_quantile_list = np.arange(.75,.9, 0.016) # 10 values
        match_counts_pool = []
        for top_quantile in top_quantile_list:
            match_counts = main(links_random, links_rf, links_ridge, links_portia, top_quantile)
            match_counts_pool.append(match_counts)
        match_counts_list = np.array(match_counts_pool).T
        match_counts_sum = np.sum(match_counts_pool, axis=0)
        #- plot
        def plot_match_counts_multiple(match_counts_list):
            utils.serif_font()
            fig, ax = plt.subplots(1, 1, tight_layout=True, figsize=(4.7, 3.5),
                                     )
            for data in match_counts_list:
                plt.plot(data)

            # ax.set_ylabel('Number of matched interactions')
            # ax.set_xticks(list(range(1, len(labels) + 1)))
            # ax.set_xticklabels(labels, rotation=0)
            # ax.set_ymargin(.25)
            # # - face colors
            # colors = ['lightpink', 'lightblue', 'lightgreen', 'cyan', 'grey']
            # for patch, color in zip(bplot['bodies'], colors):
            #     patch.set_facecolor(color)
            #     patch.set_edgecolor('black')
            #     patch.set_alpha(1)

            return fig


        plot_match_counts_multiple(match_counts_list)
        plt.show()

    section2()



