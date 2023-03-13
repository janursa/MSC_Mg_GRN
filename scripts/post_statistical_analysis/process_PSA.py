"""
Read the results of statistical analysis and extract DE proteins, imputated data, sig data, and output shortlisted data
This outputs DE data for early and late for top p values of 0.05 and 0.025
"""
import sys
import os
from typing import Tuple, List
import pandas as pd
import json
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import STATISTICAL_ANALYSIS_DIR, DATA_DIR
from edit_data.edit import correct_entry
from scripts.utils import read_write_data, process_data

def main(IMPUT_FILE:str, SIG_FILE:str, sig_t:float) -> Tuple[list,pd.DataFrame]:
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
    periods = ['day1_11', 'day1_21']
    imputs = ['MinProb', 'KNN']
    sig_t = 0.05

    DE_data = {}
    DE_proteins = {}
    for period in periods:
        for imput in imputs:
            IMPUT_FILE =  os.path.join(STATISTICAL_ANALYSIS_DIR, period, f'ProteinAbundance_tables/ProtAbundance__Norm_n_Imp_{imput}.csv')  # dir for the imputated df
            TOPTABLE_FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, period, f'DiffExp_tables/TopTable_TimeCourse_{imput}.csv')  # dir for the sig analysis

            DE_protnames_single, DE_data_single = main(IMPUT_FILE, TOPTABLE_FILE, sig_t)

            tag = '_'.join([period, imput])
            DE_proteins[tag] = DE_protnames_single
            DE_data[tag] = DE_data_single

            #- extract data in the format of arrays and store it for each model
            for study in ['ctr','mg']:
                data = process_data(DE_data_single, study=study)
                tag_study = '_'.join([tag,study])
                read_write_data(mode='write', tag=tag_study, data=data)

    with open(os.path.join(DATA_DIR, 'DE_protnames.txt'), 'w') as f:
        print(DE_proteins, file=f)
    DE_data_serialized = {ky: df.to_json() for ky, df in DE_data.items()}
    with open(os.path.join(DATA_DIR, 'DE_data.csv'), 'w') as f:
        json.dump(DE_data_serialized, f)
