"""
Read the results of statistical analysis and extract DE proteins, imputated data, sig data, and output shortlisted data
We use KNN as the chosen imputation method.
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *
from scripts.edit_data.edit import correct_entry

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
    # imputation_method = 'noImputation'
    imputation_method = 'KNN'
    sig_t = .05
    DIR_STATISTICAL_ANALYSIS =  os.path.join(MAIN_DIR,'statistical_analysis_2')
    IMPUT_FILE =  os.path.join(DIR_STATISTICAL_ANALYSIS, f'ProteinAbundance_tables/ProtAbundance__Norm_n_Imp_{imputation_method}.csv')  # dir for the imputated df
    SIG_FILE = os.path.join(DIR_STATISTICAL_ANALYSIS, f'DiffExp_tables/TopTable_TimeCourse_{imputation_method}.csv')  # dir for the sig analysis

    DE_protnames, DE_data = main(IMPUT_FILE, SIG_FILE, sig_t)
    DE_protnames_str = ''
    for gene in DE_protnames:
        DE_protnames_str += gene + ','
    DIR = os.path.join(OUTPUT_DIR, 'data')
    with open(os.path.join(DIR, 'DE_protnames.txt'), 'w') as f:
        print({'DE_protnames': DE_protnames, 'DE_protnames_str': DE_protnames_str}, file=f)
    print(f'DE proteins names are outputted to {DIR}')
    DE_data.to_csv(os.path.join(DIR, 'DE_data.csv'), index=False)
    print(f'Shortlisted data are outputted to {DIR}')