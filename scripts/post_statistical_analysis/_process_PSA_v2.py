"""
Read the results of statistical analysis and extract DE proteins, imputated data, sig data, and output shortlisted data
This outputs DE data for early, late, combined phases for top 30 and 50 proteins based on ranking
"""
import json
import sys
import os
import typing

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import STATISTICAL_ANALYSIS_DIR, DATA_DIR
from edit_data.edit import correct_entry

def extract_DE_data():
    pass

if __name__ == '__main__':
    #- read ranking data
    early_tag = "short.term.prod.rank"
    late_tag = "long.term.prod.rank"
    RANKLIST_FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, 'BothMethods_rankedList.csv')
    df_ranklist = pd.read_csv(RANKLIST_FILE, index_col=False)
    #- extract early and late top proteins
    def F_rank_prots(df, tag, n):
        if tag == 'combined':
            early, late = [F_rank_prots(df_ranklist, tag, n) for tag in [early_tag, late_tag]]
            sum = list(set(early + late))
            df_ranklist_short = df_ranklist.loc[df_ranklist['Protein'].isin(sum), :]
            df_ranklist_short['rank_overall'] = df_ranklist_short[early_tag] * df_ranklist_short[late_tag]
            return df_ranklist_short.sort_values('rank_overall')['Protein'].values.tolist()[0:n]
        else:
            return df.sort_values(tag)['Protein'].values.tolist()[0:n]
    DE_proteins = {'early_30': F_rank_prots(df_ranklist, early_tag, 30),
               'early_50': F_rank_prots(df_ranklist, early_tag, 50),
               'late_30': F_rank_prots(df_ranklist, late_tag, 30),
               'late_50': F_rank_prots(df_ranklist, late_tag, 50),
               'combined_30': F_rank_prots(df_ranklist, 'combined', 30),
               'combined_50': F_rank_prots(df_ranklist, 'combined', 50),
               }
    #- output DE proteins
    with open(os.path.join(DATA_DIR, 'DE_protnames.txt'), 'w') as f:
        print(DE_proteins, file=f)
    #- extract DE data by averaging over differnt imputation methods.
    period = 'day1_21' # imputated data of this folder is considered
    imputation_methods = ['KNN','MinProb']
    def F_read_file(method):
        FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, period, f'ProteinAbundance_tables/ProtAbundance__Norm_n_Imp_{method}.csv')  # dir for the imputated df
        df = pd.read_csv(FILE)
        return df
    imputed_data_stack = [F_read_file(method) for method in imputation_methods]
    imputed_data_stack_sort = [df.sort_values('Protein') for df in imputed_data_stack]
    imputed_data_conct = pd.concat(imputed_data_stack_sort)
    imputed_data = imputed_data_conct.groupby('Protein').mean().reset_index()
    #- shortlist to DE proteins
    extract_DE_data = lambda vector: imputed_data.loc[imputed_data['Protein'].isin(vector), :]

    DE_data = {ky:extract_DE_data(value) for ky, value in DE_proteins.items()}

    # [df.to_csv(os.path.join(DATA_DIR, f'DE_data_{key}.csv'), index=False) for key, df in DE_data.items()]
    DE_data_serialized = {ky: df.to_json() for ky, df in DE_data.items()}
    # [df.to_csv(os.path.join(DATA_DIR, f'DE_data_{key}.csv'), index=False) for key, df in DE_data.items()]
    with open(os.path.join(DATA_DIR, 'DE_data.csv'), 'w') as f:
        json.dump(DE_data_serialized, f)
    print(f'Shortlisted data are outputted to {DATA_DIR}')