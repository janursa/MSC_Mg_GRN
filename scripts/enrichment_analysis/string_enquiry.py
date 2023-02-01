import sys
import os
import typing
import json

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import DATA_DIR, ENRICH_DIR
import requests

string_api_url = "https://version-11-5.string-db.org/api"
def string_links(my_genes):

    output_format = "tsv-no-header"
    method = "network"

    request_url = "/".join([string_api_url, output_format, method])

    params = {

        "identifiers": "%0d".join(my_genes),  # your protein
        "species": 9606,  # species NCBI identifier
        "caller_identity": "www.awesome_app.org",  # your app name
        "show_query_node_labels":1

    }

    response = requests.post(request_url, data=params)

    for line in response.text.strip().split("\n"):

        l = line.strip().split("\t")
        p1, p2 = l[2], l[3]

        ## filter the interaction according to experimental score
        experimental_score = float(l[10])
        if experimental_score > 0.4:
            ## print
            print("\t".join([p1, p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))

def string_EA(my_genes):
    output_format = "tsv"
    method = "enrichment"
    request_url = "/".join([string_api_url, output_format, method])
    params = {
        "identifiers": "%0d".join(my_genes),  # your protein
        "species": 9606,  # species NCBI identifier
    }
    response = requests.post(request_url, data=params)
    for i, line in enumerate(response.text.strip().split("\n")):
        values = line.split("\t")
        if i == 0:
            df = pd.DataFrame(columns=values)
        else:
            df.loc[len(df.index)] = values
    # - columns to retain
    cols_to_stay = ['category', 'number_of_genes_in_background', 'number_of_genes', 'inputGenes',
                    'fdr', 'description']
    df = df.loc[:, cols_to_stay]
    #- rename columns
    df.rename(columns={'fdr':'FDR', 'description': 'Description', 'number_of_genes':'ProteinCount',
                       'inputGenes':'ProteinNames', 'number_of_genes_in_background':'BackgroundGeneCount'}, inplace=True)
    #- add new columns
    ProteinCount = np.asarray(list(map(int, df['ProteinCount'].values)))
    BackgroundGeneCount = np.asarray(list(map(int, df['BackgroundGeneCount'].values)))
    # df.loc[:,'ProteinRatio'] = ProteinCount/len(my_genes)
    df.loc[:, 'Strength'] = -1/np.log10(ProteinCount/BackgroundGeneCount)
    #- change category names
    df.replace('RCTM','Reactome Pathways', inplace=True)
    df.replace('Keyword', 'Uniprot Keywords', inplace=True)
    df.replace('Process', 'GO Biological Process', inplace=True)
    df.replace('Function', 'GO Biological Function', inplace=True)
    df.replace('Component', 'GO Cellular Component', inplace=True)
    return df


if __name__ == '__main__':
    with open(os.path.join(DATA_DIR, 'DE_protnames.txt'), 'r') as f:
        data = eval(f.read())
        DE_proteins_early = data['early']
        DE_proteins_late = data['late']
    for phase, DE_proteins in zip(['early', 'late'], [DE_proteins_early, DE_proteins_late]):
        df_enrich = string_EA(DE_proteins)
        df_enrich.to_csv(os.path.join(ENRICH_DIR, f'enrichment_all_{phase}.csv'), index=False)

