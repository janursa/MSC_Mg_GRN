#- retreive the links

# links_ctr_grni, links_sample_grni = retreive('ctr'), retreive_grni('mg')
# links_ctr_portia, links_sample_portia = retreive_portia('ctr'), retreive_portia('mg')
# #- normalize
# def normalize(links, study):
# #     links['Weight'] = links['Weight']/np.std(links['Weight'])
#     #- normalize pool of links
#     if study=='grni':
#         links['WeightPool'] = [row/np.std(links['Weight']) for row in links['WeightPool'].values]
#         links['Weight'] = [np.mean(row) for row in links['WeightPool'].values]        
#     if study =='portia':
#         links['Weight'] = links['Weight']/np.std(links['Weight'])
#     #- update the weight
#     return links
# links_ctr_grni, links_sample_grni = normalize(links_ctr_grni,'grni'), normalize(links_sample_grni,'grni')
# links_ctr_portia, links_sample_portia = normalize(links_ctr_portia,'portia'), normalize(links_sample_portia,'portia')

# #- combine
# def combine(links_1, links_2):
#     links = links_1.copy() #you bitch
#     keys = ['Weight']
#     for key in keys:
#         links[key] = (links_1[key] + links_2[key])/2
#     return links
# links_ctr_ensemble = combine(links_ctr_grni, links_ctr_portia)
# links_sample_ensemble = combine(links_sample_grni, links_sample_portia)
# #- save
# links_ctr_ensemble.to_pickle(os.path.join(OUTPUT_DIR,'GRN', 'links_ctr_ensemble.csv')) # save with pickle because there is a list of items in data
# links_sample_ensemble.to_pickle(os.path.join(OUTPUT_DIR,'GRN', 'links_mg_ensemble.csv'))