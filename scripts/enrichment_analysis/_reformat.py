"""
Map different annotations to each other.
String data should be available beforehand.
# geneID (entrenz gene id). https://www.uniprot.org/id-mapping/
"""
import sys
import os
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import *


def reireive_string_data_F(enrich_file: str, enrich_type: str, map_genename_protname: dict, term_length_cut_threshold: int) -> pd.DataFrame:
    '''
        Reads model_selection EA data such as enrichment.Function
    '''
    enrich_data = pd.read_csv(enrich_file, sep='\t', index_col=False)
    # keep the specified enrichment type
    filter_ = enrich_data.loc[:, '#category'] == enrich_type
    enrich_data = enrich_data.loc[filter_, :]
    enrich_data.drop('#category', axis=1, inplace=True)
    # remove some columns
    enrich_data.drop('matching proteins in your network (IDs)', axis=1, inplace=True)
    enrich_data.drop('background gene count', axis=1, inplace=True)
    # - redefine column names
    enrich_data.columns = ['ID', 'Description', 'ProteinCount', 'Strength', 'FDR', 'genenames']
    # replace model_selection gene name with uniprotname
    print(enrich_data.columns)
    unitprotnames = []
    for row in enrich_data['genenames']:
        a = ''
        for item in row.split(','):
            a += map_genename_protname[item] + ','
        unitprotnames.append(a[0:-1])
    enrich_data['protnames'] = unitprotnames
    enrich_data.drop('genenames', axis=1, inplace=True)

    #- fill the rest of the entries
    enrich_data.loc[:, 'ProteinRatio'] = enrich_data['ProteinCount'] / len(map_genename_protname)
    # not useful but required as input to the graph
    enrich_data['pvalue'] = enrich_data['FDR']
    enrich_data['qvalue'] = enrich_data['pvalue']
    #- shorten term names for visualization purposes
    short_descript = []
    for item in enrich_data.loc[:, 'Description']:
        if len(item) > term_length_cut_threshold:
            aa = item[0:term_length_cut_threshold] + '...'
        else:
            aa = item
        short_descript.append(aa)
    enrich_data['Description'] = short_descript
    return enrich_data

if __name__ == '__main__':
    #- retreiev map_genename_protname
    with open(os.path.join(OUTPUT_DIR,'data','map_genename_protname.json')) as file:
        map_genename_protname = json.load(file)['map']
    #- retreieve enrichmenetn analysis results
    FILE = os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'enrichment.all.tsv')
    #- setting
    term_length_cut_threshold = 10 # letters
    #- process
    data_process = reireive_string_data_F(FILE, 'GO Process', map_genename_protname, term_length_cut_threshold)
    data_process.to_csv(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'enrichment_process.csv'), index=False)
    print('number of enriched process tersm ', len(data_process))

    data_function = reireive_string_data_F(FILE, 'GO Function', map_genename_protname, term_length_cut_threshold)
    data_function.to_csv(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'enrichment_function.csv'), index=False)
    print('number of enriched function tersm ', len(data_function))

    data_component = reireive_string_data_F(FILE, 'GO Component', map_genename_protname, term_length_cut_threshold)
    data_component.to_csv(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'enrichment_component.csv'), index=False)
    print('number of enriched component tersm ', len(data_component))

    data_uniprot = reireive_string_data_F(FILE, 'UniProt Keywords', map_genename_protname, term_length_cut_threshold)
    data_uniprot.to_csv(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'enrichment_uniprot.csv'))
    print('number of enriched UniProt keywords terms ', len(data_uniprot))

    data_reactome = reireive_string_data_F(FILE, 'Reactome', map_genename_protname, term_length_cut_threshold)
    data_reactome.to_csv(os.path.join(OUTPUT_DIR, 'enrichment_analysis', 'enrichment_reactome.csv'))
    print('number of enriched Reactome tersm ', len(data_reactome))