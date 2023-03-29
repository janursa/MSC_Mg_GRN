import sys
import os
import pandas as pd
import numpy as np
from typing import Dict, List, Optional
import argparse
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from imports import MAIN_DIR, DATA_DIR, time_points

def correct_entry(df: pd.DataFrame, col_name: str='Entry') -> pd.DataFrame:
    """
    Replace incorrect entry names in the data with the correct versions.
    """
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
    parser = argparse.ArgumentParser()

    parser.add_argument('--original_df_dir', type=str, default=Path(MAIN_DIR)/'data'/'original_omics.xlsx',
                        help='directory for original data')
    parser.add_argument('--df_dir', type=str, default=Path(DATA_DIR)/'edited_data.csv',
                        help='directory for processed data')
    args = parser.parse_args()

    original_df_dir = args.original_df_dir
    df_dir = args.df_dir

    if not os.path.isdir(DATA_DIR):
        os.makedirs(DATA_DIR)
    #- read the original data 
    o_df = pd.read_excel(original_df_dir)
    #- process the data
    o_df_m = rename_missing_symbols(o_df, output_dir=DATA_DIR)
    df = tailor_names(o_df_m, days=time_points())
    df = remove_zeros(df, days=time_points())
    #- entry corrections
    df = correct_entry(df)
    df.to_csv(df_dir)
