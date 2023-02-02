"""
Read the results of statistical analysis and extract DE proteins, imputated data, sig data, and output shortlisted data
We use KNN as the chosen imputation method.
"""
import sys
import os
import typing

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import STATISTICAL_ANALYSIS_DIR, DATA_DIR
from edit_data.edit import correct_entry

def main(IMPUT_FILE:str, SIG_FILE:str, sig_t:float) -> typing.Tuple[list,pd.DataFrame]:
    # load the imputed df
    df_imput = pd.read_csv(IMPUT_FILE)
    df_imput.rename(columns={'Unnamed: 0': 'Protein'}, inplace=True)
    df_imput = correct_entry(df_imput, 'Protein')
    # read the sig results
    df_sig = pd.read_csv(SIG_FILE, index_col=False)
    df_sig = correct_entry(df_sig, 'Protein')
    df_sig.drop(columns=['Unnamed: 0'], inplace=True)
    # sort the data to make them comparable
    df_imput.sort_values('Protein', inplace=True)
    df_imput.reset_index(drop=True, inplace=True)
    df_sig.sort_values('Protein', inplace=True)
    df_sig.reset_index(drop=True, inplace=True)
    # determine sign proteins

    sig_flags = (df_sig['P.Value'] < sig_t).tolist()
    print(f'Number of DE proteins {sig_flags.count(True)}')
    # output the DE proteins to a file
    DE_protnames = list(df_sig.loc[sig_flags, 'Protein'].values)
    # shortlisted DE dataframe
    DE_data = df_imput.loc[df_imput['Protein'].isin(DE_protnames), :].reset_index(drop=True)
    return DE_protnames, DE_data

if __name__ == '__main__':
    #- read ranking data
    early_rank_col_tag = "short.term.prod.rank"
    late_rank_col_tag = "long.term.prod.rank"
    RANKLIST_FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, 'BothMethods_rankedList.csv')
    df_ranklist = pd.read_csv(RANKLIST_FILE, index_col=False)
    #- extract early and late top proteins
    top_n = 25
    def F_rank_prots(df, tag, n):
        return df.sort_values(tag)['Protein'].values.tolist()[0:n]
    DE_proteins_early, DE_proteins_late = [F_rank_prots(df_ranklist, tag, top_n) for tag in [early_rank_col_tag, late_rank_col_tag]]
    #- from top 50 of early and late, choose top 50 of overall
    top_n = 50
    early, late = [F_rank_prots(df_ranklist, tag, top_n) for tag in [early_rank_col_tag, late_rank_col_tag]]
    sum = list(set(early + late))
    df_ranklist_short = df_ranklist.loc[df_ranklist['Protein'].isin(sum),:]
    df_ranklist_short['rank_overall'] = df_ranklist_short[late_rank_col_tag]*df_ranklist_short[early_rank_col_tag]
    DE_proteins_overall = F_rank_prots(df_ranklist_short, 'rank_overall', top_n)
    #- output DE proteins
    DE_proteins_early_str, DE_proteins_late_str, DE_proteins_overall_str = list(map(lambda vector: ','.join(vector), (DE_proteins_early, DE_proteins_late, DE_proteins_overall)))
    with open(os.path.join(DATA_DIR, 'DE_protnames.txt'), 'w') as f:
        print({'early': DE_proteins_early,
               'late': DE_proteins_late,
               'early_str': DE_proteins_early_str,
               'late_str': DE_proteins_late_str,
               'overall': DE_proteins_overall,
               'overall_str': DE_proteins_overall_str,
               }, file=f)
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
    DE_data = imputed_data.loc[imputed_data['Protein'].isin(DE_proteins_overall),:]
    DE_data.to_csv(os.path.join(DATA_DIR, 'DE_data.csv'), index=False)
    print(f'Shortlisted data are outputted to {DATA_DIR}')