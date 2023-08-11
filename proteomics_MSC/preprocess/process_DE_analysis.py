import sys
import os
from typing import Tuple, List
import pandas as pd
import json
from proteomics_MSC.common_tools import read_write_data, process_data
from proteomics_MSC.imports import DE_data_file, DE_protnames_file, DE_common_proteins_file, STATISTICAL_ANALYSIS_DIR
from proteomics_MSC.preprocess.edit_data import correct_entry
import yaml
from importlib.resources import open_text

with open_text('proteomics_MSC', 'config.yaml') as file:
    config = yaml.safe_load(file)

def load_imputated_df(imput_method, period_days):
    IMPUT_FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, period_days,
                              f'ProteinAbundance_tables/ProtAbundance__Norm_n_Imp_{imput_method}.csv')  # dir for the imputated df
    df_imput = pd.read_csv(IMPUT_FILE)
    df_imput.rename(columns={'Unnamed: 0': 'Protein'}, inplace=True)
    df_imput = correct_entry(df_imput, 'Protein')
    return df_imput
def load_sig_df(imput_method, period_days):
    TOPTABLE_FILE = os.path.join(STATISTICAL_ANALYSIS_DIR, period_days,
                                 f'DiffExp_tables/TopTable_TimeCourse_{imput_method}.csv')  # dir for the sig analysis
    df_sig = pd.read_csv(TOPTABLE_FILE, index_col=False)
    df_sig = correct_entry(df_sig, 'Protein')
    df_sig.drop(columns=['Unnamed: 0'], inplace=True)
    return df_sig


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
def post_process_DE_analysis() -> tuple[dict[str, list[str]], dict[str, pd.DataFrame]]:
    """Extract DE proteins and associated data from statistical analysis"""
    sig_t = 0.05
    DE_data: dict[str, pd.DataFrame] = {}
    DE_proteins: dict[str, list[str]] = {}
    for period_phase, period in zip(periods, periods_days):
        for imput in imputs:
            DE_protnames_single, DE_data_single = extract_de_data(imput, period, sig_t)
            tag = '_'.join([period_phase, imput])
            DE_proteins[tag] = DE_protnames_single
            DE_data[tag] = DE_data_single
            # - extract data in the format of arrays and store it for each model
            for study in studies:
                # - get the data as matrix
                data_matrix = process_data(DE_data_single, study=study)
                tag_study = '_'.join([tag, study])
                read_write_data(mode='write', tag=tag_study, data=data_matrix)
    return DE_proteins, DE_data
if __name__ == '__main__':
    periods_days = config['periods_days']
    periods = config['periods']
    imputs = config['imputs']
    studies = config['studies']


    DE_proteins, DE_data = post_process_DE_analysis()
    with open(DE_protnames_file, 'w') as f:
        print(DE_proteins, file=f)
    DE_data_serialized = {ky: df.to_json() for ky, df in DE_data.items()}
    with open(DE_data_file, 'w') as f:
        json.dump(DE_data_serialized, f)
    print(f'Output {DE_data_file}')
    # finding the proteins are that are mutually detected as DE across all datasets
    sets = [set(val) for val in DE_proteins.values()]
    # Find the intersection of all sets
    common_proteins = list(set.intersection(*sets))
    with open(DE_common_proteins_file, 'w') as f:
        print(common_proteins, file=f)
    print(f'Output {DE_common_proteins_file}')
    print(f'common DE proteins {common_proteins}')



    #- to store the data as a df where gene names are given in the first column
    if False:
        for period_phase, period in zip(periods, periods_days):
            for imput in imputs:
                tag = '_'.join([period_phase, imput])
                #- extract data in the format of arrays and store it for each model
                for study in studies:
                    #- get the data for each study and add gene name to the df
                    data_df = DE_data[tag].filter(like=study)
                    data_df.insert(0, 'Protein', DE_proteins[tag])
                    data_df.round(3).to_csv(f'{DATA_DIR}/df_data_{tag}_{study}.csv')


