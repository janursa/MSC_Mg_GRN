import os
import pandas as pd
from edit_data.edit import correct_entry
from imports import STATISTICAL_ANALYSIS_DIR
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