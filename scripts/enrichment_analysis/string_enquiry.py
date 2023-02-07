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
def string_enriched_network(my_genes):

    output_format = "tsv"
    method = "network"

    request_url = "/".join([string_api_url, output_format, method])

    params = {

        "identifiers": "%0d".join(my_genes),  # your protein
        "species": 9606,  # species NCBI identifier
        "network_type": "functional",
        "show_query_node_labels":1 #when available use submitted names in the preferredName column when (0 or 1) (default:0)

    }

    response = requests.post(request_url, data=params)
    for i, line in enumerate(response.text.strip().split("\n")):
        values = line.split("\t")
        if i == 0:
            df = pd.DataFrame(columns=values)
        else:
            df.loc[len(df.index)] = values
    # - columns to retain
    cols_to_stay = ['preferredName_A', 'preferredName_B', 'score']
    df = df.loc[:, cols_to_stay]
    # - rename columns
    df.rename(columns={'preferredName_A': 'Regulator', 'preferredName_B': 'Target', 'score':'Weight'},
              inplace=True)
    return df
def string_enriched_functions(my_genes):
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
    df.loc[:, 'Strength'] = -1/(np.log10(ProteinCount/BackgroundGeneCount)-0.000001)
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
    #- functional enrichment
    # for DE_type, DE_proteins in data.items():
    #     if DE_type in ['early_30', 'late_30']:
    #         df_enrich = string_enriched_functions(DE_proteins)
    #         df_enrich.to_csv(os.path.join(ENRICH_DIR, f'enrichment_all_{DE_type}.csv'), index=False)

    #- network enrichment
    for DE_type, DE_proteins in data.items():
        df_enrich = string_enriched_network(DE_proteins)
        df_enrich.to_csv(os.path.join(ENRICH_DIR, f'network_{DE_type}.csv'), index=False)

