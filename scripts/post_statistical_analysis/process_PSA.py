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
                #- get the data as matrix
                data_matrix = process_data(DE_data_single, study=study)
                tag_study = '_'.join([tag,study])
                read_write_data(mode='write', tag=tag_study, data=data_matrix)
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


    #- to store the data as a df where gene names are given in the first column
    for period_phase, period in zip(periods, periods_days):
        for imput in imputs:
            tag = '_'.join([period_phase, imput])
            #- extract data in the format of arrays and store it for each model
            for study in studies:
                #- get the data for each study and add gene name to the df
                data_df = DE_data[tag].filter(like=study)
                data_df.insert(0, 'Protein', DE_proteins[tag])
                data_df.round(3).to_csv(f'{DATA_DIR}/df_data_{tag}_{study}.csv')

    # - check if certain proteins are up or down regulated for a given model
    def get_values(df, protein, day):
        protein_data = df[df['Protein'] == protein]

        if not protein_data.empty:
            ctr_value = protein_data[f'ctr_{day}'].values[0]
            mg_value = protein_data[f'mg_{day}'].values[0]
            return ctr_value, mg_value
        else:
            return "Protein not found"


    def get_mean_values(df, protein):
        protein_data = df[df['Protein'] == protein]

        if not protein_data.empty:
            ctr_cols = [col for col in df.columns if col.startswith('ctr')]
            mg_cols = [col for col in df.columns if col.startswith('mg')]

            ctr_mean = protein_data[ctr_cols].mean(axis=1).values[0]
            mg_mean = protein_data[mg_cols].mean(axis=1).values[0]

            return ctr_mean, mg_mean
        else:
            return "Protein not found"
    def print_values(df, protein):
        protein_data = df[df['Protein'] == protein]


        ctr_cols = [col for col in df.columns if col.startswith('ctr')]
        mg_cols = [col for col in df.columns if col.startswith('mg')]

        print('ctr', protein_data[ctr_cols].values[0])
        print('mg', protein_data[mg_cols].values[0])

    data = DE_data['late_KNN']
    proteins = {'P40926':'MDH2', "O94925": "GLS", "Q13263": "TRIM28"}
    day = '11'
    for protein, gene in proteins.items():
        # print(f'{gene}:{get_mean_values(data, protein)}')
        # print(f'{data[]}')
        print(gene)
        print_values(data, protein)

