"""
    Sets of functions useful for sensitivity analysis (robustness analysis)
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
import utils
def multiGaussNoise(links_ctr, links_sample, protnames, target_prots, OUTPUT_DIR):
    """
        Multi step function to conduct multiplicative gaussian noise analysis
    """
    #- create noised dfs
    links_ctr_noised = sensitivity_analysis.MG_noise_F(links_ctr)
    links_sample_noised = sensitivity_analysis.MG_noise_F(links_sample)
    #- run VSA for each df
    series_1 = sensitivity_analysis.VSA_F(links_ctr_noised, protnames, target_prots)
    series_2 = sensitivity_analysis.VSA_F(links_sample_noised, protnames, target_prots)
    #- filter series for target proteins and sort them as a dict
    rr_sorted = sensitivity_analysis.filter_and_sort(series_1, series_2, target_prots)
    #- plot
    VSA.plot_SA(rr_sorted, name='multiGaussNoise', OUTPUT_DIR=OUTPUT_DIR)
    
def addGaussNoise(links_ctr, links_sample, protnames, target_prots, OUTPUT_DIR):
    """
        Multi step function to conduct additive gaussian noise analysis
    """
    #- create noised dfs
    links_ctr_noised = sensitivity_analysis.AG_noise_F(links_ctr)
    links_sample_noised = sensitivity_analysis.AG_noise_F(links_sample)
    #- run VSA for each df
    series_1 = sensitivity_analysis.VSA_F(links_ctr_noised, protnames, target_prots)
    series_2 = sensitivity_analysis.VSA_F(links_sample_noised, protnames, target_prots)
    #- filter series for target proteins and sort them as a dict
    rr_sorted = sensitivity_analysis.filter_and_sort(series_1, series_2, target_prots)
    #- plot
    VSA.plot_SA(rr_sorted, name='addGaussNoise', OUTPUT_DIR=OUTPUT_DIR)

def VSA_F(links_series, protnames, target_prots):
    """ 
        Runs VSA for the entire noised dfs 
    """
    rr_series = []
    for i in range(len(links_series)):
        links = links_series[i]
        rr = VSA.analyse(links, protnames)
        rr_t = rr.loc[rr['Entry'].isin(target_prots),:]
        rr_series.append(rr)
    return rr_series
def filter_and_sort(series_1, series_2, target_prots):
    """
        Filters series based on target proteins and puts the data in a form of dict
    """
    sorted_VSA = {}
    extract = lambda VSA_s, prot, tag='AS': [VSA.loc[VSA['Entry']==prot,:][tag].values.tolist()[0] for VSA in VSA_s]
    for prot in target_prots:
        ASs_ctr = extract(series_1, prot=prot, tag='AS')
        PSs_ctr = extract(series_1, prot=prot, tag='PS')
        Roles_ctr = extract(series_1, prot=prot, tag='Role')
        
        ASs_mg = extract(series_2, prot=prot, tag='AS')
        PSs_mg = extract(series_2, prot=prot, tag='PS')
        Roles_mg = extract(series_2, prot=prot, tag='Role')
        
        sorted_VSA[prot] = {'ctr':{'AS':ASs_ctr,'PS':PSs_ctr, 'Role':Roles_ctr},'mg':{'AS':ASs_mg,'PS':PSs_mg, 'Role':Roles_mg}}
    return sorted_VSA
def MG_noise_F(links):
    n = 100 # noise replica 
    std = 0.15 # standard deviation
    noised_link = links.copy()
    noised_links = [noised_link.assign(Weight= noised_link['Weight']*np.random.normal(loc=1, scale=std, size=len(links)))
                                       for i in range(n)]
    return noised_links
def AG_noise_F(links):
    n = 100 # noise replica 
    std = 0.5 # standard deviation
    noised_link = links.copy()
    noised_links = [noised_link.assign(Weight= noised_link['Weight']+np.random.normal(loc=1, scale=std, size=len(links)))
                                       for i in range(n)]
    return noised_links