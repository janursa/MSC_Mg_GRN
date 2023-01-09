from imports import *

def execute(method):
    """
    reads the links from the file, filters and outputs them.
    """

    links_ctr, links_mg = retreive_grn('ctr',method), retreive_grn('mg',method)
    print(f'number of original links; ctr: {len(links_ctr)}, mg:{len(links_mg)}')
    if method=='rf': # oob score
        links_ctr = utils.Links.filter_fitscore(links_ctr)
        links_mg = utils.Links.filter_fitscore(links_mg)
        print(f'number of links after filtering based on fit score; ctr: {len(links_ctr)}, mg:{len(links_mg)}')
    links_short_ctr = utils.Links.choose_top_count(links_ctr, count=100)
    links_short_mg = utils.Links.choose_top_count(links_mg, count=100)
    print(f'number of links after filter; ctr: {len(links_short_ctr)}, mg:{len(links_short_mg)}')
    links_short_ctr.to_pickle(os.path.join(OUTPUT_DIR,'GRN',f'links_short_ctr_{method}.csv'))
    links_short_mg.to_pickle(os.path.join(OUTPUT_DIR,'GRN',f'links_short_mg_{method}.csv'))
#-- portia 
execute('portia')
#-- rf 
execute('rf')
#-- ridge 
execute('ridge')