"""
Read the results of statistical analysis and extract DE proteins, imputated data, sig data, and output shortlisted data
This outputs DE data for early and late for top p values of 0.05 and 0.025
"""
import sys
import os
from typing import Tuple, List
import pandas as pd
import json
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import DATA_DIR
from common_tools import read_write_data, process_data
from psa_aux import load_sig_df, load_imputated_df
def extract_de_data(imput_method, period_days, sig_t:float) -> Tuple[list,pd.DataFrame]:
    """ Extracts DE proteins and data
    """
    # load the imputed df
    df_imput = load_imputated_df(imput_method, period_days)
    # read the sig results
    df_sig = load_sig_df(imput_method, period_days)
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--periods', nargs='+', default=['early', 'late'])
    parser.add_argument('--periods_days', nargs='+', default=['day1_11', 'day1_21'])
    parser.add_argument('--imputs', nargs='+', default=['MinProb', 'KNN'])
    parser.add_argument('--studies', nargs='+', default=['ctr', 'mg'])
    parser.add_argument('--sig_t', nargs='+', default=0.05)
    args, remaining_args = parser.parse_known_args()

    periods = args.periods
    periods_days = args.periods_days
    imputs = args.imputs
    studies = args.studies
    sig_t = args.sig_t

    DE_data = {}
    DE_proteins = {}
    for period_phase, period in zip(periods, periods_days):
        for imput in imputs:
            DE_protnames_single, DE_data_single = extract_de_data(imput, period, sig_t)
            tag = '_'.join([period_phase, imput])
            DE_proteins[tag] = DE_protnames_single
            DE_data[tag] = DE_data_single
            #- extract data in the format of arrays and store it for each model
            for study in studies:
                data = process_data(DE_data_single, study=study)
                tag_study = '_'.join([tag,study])
                read_write_data(mode='write', tag=tag_study, data=data)
    # finding the proteins are that are mutually detected as DE across all datasets
    sets = [set(val) for val in DE_proteins.values()]
    # Find the intersection of all sets
    common_proteins = list(set.intersection(*sets))
    print(f'common DE proteins {common_proteins}')
    with open(os.path.join(DATA_DIR, 'DE_protnames.txt'), 'w') as f:
        print(DE_proteins, file=f)
    DE_data_serialized = {ky: df.to_json() for ky, df in DE_data.items()}
    with open(os.path.join(DATA_DIR, 'DE_data.csv'), 'w') as f:
        json.dump(DE_data_serialized, f)
