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
    top_n = 25
    #TODO:
    # [1] read 'BothMethods_rankedList.csv', choose the top_n for early and late phase -> DE_proteins_early and DE_proteins_late
    # [2] extract ranking of top_n proteins from early and late. multiply them and to choose top_n final DE_proteins
    # [3] export DE_proteins
    # [4] read imputated results from two methods 'KNN' and 'prob'
    # [4] extract DE_proteins_early and DE_proteins_late from the imputated results and average them
    # [5] export DE_data.csv

    #- [1]
    early_rank_col_tag = "short.term.prod.rank"
    late_rank_col_tag = "long.term.prod.rank"
    RANKLIST_FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, 'BothMethods_rankedList.csv')
    df_ranklist = pd.read_csv(RANKLIST_FILE, index_col=False)
    def F_rank_prots(df, tag):
        return df.sort_values(tag)['Protein'].values.tolist()[0:top_n] # top 5 percent
    DE_proteins_early, DE_proteins_late = [F_rank_prots(df_ranklist, tag) for tag in [early_rank_col_tag, late_rank_col_tag]]
    #- [2]
    DE_proteins_sum = list(set(DE_proteins_early + DE_proteins_late))
    df_ranklist_short = df_ranklist.loc[df_ranklist['Protein'].isin(DE_proteins_sum),:]
    df_ranklist_short['rank_overall'] = df_ranklist_short[late_rank_col_tag]*df_ranklist_short[early_rank_col_tag]
    DE_proteins_overall = F_rank_prots(df_ranklist_short, 'rank_overall')
    #- [3]
    DE_proteins_early_str, DE_proteins_late_str, DE_proteins_overall_str = list(map(lambda vector: ','.join(vector), (DE_proteins_early, DE_proteins_late, DE_proteins_overall)))
    with open(os.path.join(DATA_DIR, 'DE_protnames.txt'), 'w') as f:
        print({'early': DE_proteins_early,
               'late': DE_proteins_late,
               'early_str': DE_proteins_early_str,
               'late_str': DE_proteins_late_str,
               'overall': DE_proteins_overall,
               'overall_str': DE_proteins_overall_str,
               }, file=f)
    #- [4]

    aa
    imputation_method = 'KNN'
    sig_t = .05
    IMPUT_FILE =  os.path.join(STATISTICAL_ANALYSIS_DIR, f'ProteinAbundance_tables/ProtAbundance__Norm_n_Imp_{imputation_method}.csv')  # dir for the imputated df
    SIG_FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, f'DiffExp_tables/TopTable_TimeCourse_{imputation_method}.csv')  # dir for the sig analysis

    DE_protnames, DE_data = main(IMPUT_FILE, SIG_FILE, sig_t)

    with open(os.path.join(DATA_DIR, 'DE_protnames.txt'), 'w') as f:
        print({'DE_protnames': DE_protnames, 'DE_protnames_str': DE_protnames_str}, file=f)
    print(f'DE proteins names are outputted to {DATA_DIR}')
    DE_data.to_csv(os.path.join(DATA_DIR, 'DE_data.csv'), index=False)
    print(f'Shortlisted data are outputted to {DATA_DIR}')