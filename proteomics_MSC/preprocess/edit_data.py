import sys
import os
import pandas as pd
import numpy as np
# import importlib_resources as resources
from importlib.resources import open_text
from typing import List, Dict, Tuple, Optional
import yaml
from proteomics_MSC.imports import DATA_DIR, original_omics_file, processed_omics_file

with open_text('proteomics_MSC', 'config.yaml') as config_file:
    config = yaml.safe_load(config_file)

def correct_entry(df: pd.DataFrame, col_name: str='Entry') -> pd.DataFrame:
    """Replace incorrect entry names in the data with the correct versions."""
    entry_corrections: Dict[str, str] = {'P81644': 'P02652'}
    df[col_name].replace(entry_corrections, inplace=True)
    return df
def tailor_names(o_df: pd.DataFrame, days: List[int]) -> pd.DataFrame:
    """Change column names to make them easier to read."""
    df: pd.DataFrame = o_df[['Entry', 'Gene names  (primary)']].copy()
    df['Gene names  (primary)'].fillna(df['Entry'], inplace=True)
    # Rename the columns for ctr and sample
    for i in days:
        df[f'ctr_{i}'] = o_df[f'LogLFQ intensity {i}_0']
        df[f'mg_{i}'] = o_df[f'LogLFQ intensity {i}_1']
    print('Original data size: {}'.format(len(df)))
    return df
def remove_zeros(df: pd.DataFrame, days: List[int]) -> pd.DataFrame:
    """Drop rows with all zeros in the specified columns."""
    cols = [f'ctr_{i}' for i in days] + [f'mg_{i}' for i in days]
    df = df.replace(0, np.nan).dropna(subset=cols, how='all')
    df.reset_index(drop=True, inplace=True)
    print('Data size after removing rows with all zeros: {}'.format(len(df)))
    return df
def rename_missing_symbols(df: pd.DataFrame, p_name: str = 'Entry',
                           output_dir: Optional[str] = None) -> pd.DataFrame:
    """Rename missing symbols in protein names."""
    df_c = df.copy()
    nulls = df_c[p_name].isnull()
    unique_names = [f'p_{i}' for i in range(sum(nulls))]
    map_ = {i: name for i, name in zip(nulls.index[nulls], unique_names)}
    df_c.loc[nulls, p_name] = unique_names
    remaining_nulls = df_c[p_name].isnull()
    if remaining_nulls.any():
        print('Warning: {} missing names remain.'.format(remaining_nulls.sum()))
    if output_dir:
        with open(os.path.join(output_dir, 'missing_names.txt'), 'w') as f:
            for key, value in map_.items():
                f.write('Original index {}: {}\n'.format(key, value))
    return df_c
if __name__ == '__main__':
    if not os.path.isdir(DATA_DIR):
        os.makedirs(DATA_DIR)
    #- read the original data
    o_df = pd.read_excel(original_omics_file)
    #- process the data
    o_df_m = rename_missing_symbols(o_df, output_dir=DATA_DIR)
    df = tailor_names(o_df_m, days=config['time_points'])
    df = remove_zeros(df, days=config['time_points'])
    #- entry corrections
    df = correct_entry(df)
    print(f'number of proteins {len(df)}')
    print(df.columns)
    df.to_csv(processed_omics_file)

    print(f'ourputs are written to {DATA_DIR}')
