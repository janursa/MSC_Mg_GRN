import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import *

specs = dict(
    o_df_dir = os.path.join(MAIN_DIR,'data','original_omics.xlsx'),
    df_dir = os.path.join(MAIN_DIR,'results/data','edited_data.csv'),
    time = time,
    p_name = 'Entry',  # The column name in the original data to be used as protein name
    c_tag = 'ctr_',
    s_tag = 'mg_',
    c_func = lambda t: 'ctr_' + str(t),  # The tag used in the tailored dataframe for control
    s_func = lambda t: 'mg_' + str(t),  # The tag used in the tailored dataframe for sample
    o_c_func = lambda t: 'LogLFQ intensity ' + str(t) + '_0',  # Func to extract the control data from the original database 
    o_s_func = lambda t: 'LogLFQ intensity ' + str(t) + '_1',  # Func to extract the sample data from the original database
)
def correct_entry(df, entry_name='Entry'):
    # - protnames needs correction: the original protnames are not homo sapiens
    entry_corrections = {'P81644': 'P02652'}
    for key in entry_corrections.keys():
        if key in df[entry_name].values:
            df.loc[df.loc[:,entry_name]==key,entry_name] = entry_corrections[key]
    return df

if __name__ == '__main__':
    #- read the original data 
    o_df = pd.read_excel(specs['o_df_dir'])
    #- process the data
    o_df_m = utils.rename_missing_symbols(o_df, OUTPUT_DIR=os.path.join(OUTPUT_DIR,'data'));
    df = utils.tailor_names(o_df_m, **specs)
    df = utils.remove_zeros(df, **specs)
    #- entry corrections
    df = correct_entry(df)
    df.to_csv(specs['df_dir'])
