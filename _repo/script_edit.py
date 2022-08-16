#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:40:58 2022
Processes the omics data file, edits the missing information, relabels the columns and outputs the clean data
@author: matin
"""
import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import copy
import math
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file)
sys.path.insert(0,main_dir)

## read the data
original_df_dir = os.path.join(main_dir,'data','original_omics.xlsx')
df_dir = os.path.join(main_dir,'data','omics.csv')
if 'original_df' not in locals():
    print("data not found")
    original_df = pd.read_excel(original_df_dir)

## tailor the data  
time = [1,2,3,4,7,8,10,11] # in dayx
# time = [1,2,3] # in dayx
prefix = 'LogLFQ intensity '
df = pd.DataFrame()
df['Symbols'] = original_df['Gene names  (primary )']
df['ID'] = df['Symbols']

# db['Sites'] = ""
# db['Effect'] = ""
for ii in time: #rename the data for control
    tag = prefix + str(ii) + '_0'
    df['ctr_{}'.format(ii)] = original_df[tag]
# print(db.loc[[40]])
for ii in time: #rename the data for Mg
    tag = prefix + str(ii) + '_1'
    df['Mg_{}'.format(ii)] = original_df[tag]

## rename the missing symbols
df['Symbols'].replace('', np.nan, inplace=True) # to rename the empty lines
jj = 1
for ii in range(len(df['Symbols'])): # assign names to the missng names
    try:
        if math.isnan(df['Symbols'][ii]):
            df['Symbols'][ii] = 'tag_{}'.format(jj)
            jj+=1
    except:
        pass
# df.dropna(subset=['Symbols'], inplace=True) 

df.loc[~(df==0).all(axis=1)] # drop those rows with all zeros

start_ii = 1
if 'ID' in df:
    start_ii += 1

for ii in range(1,len(df)): # fill thoes single zeros with neighbor extrapolation
    for jj in range(start_ii,len(df.columns)):
        if df.iloc[ii,jj] == 0:
            if jj == start_ii: # only right hand side
                rr = df.iloc[ii,jj+1]
            elif jj == len(df.columns)-1: #only left hand side
                rr = df.iloc[ii,jj-1]
            else:
                if (df.iloc[ii,jj-1]==0) or (df.iloc[ii,jj+1]==0):
                    rr = 0
                else:
                    rr = (df.iloc[ii,jj-1]+df.iloc[ii,jj+1])/2
            df.iloc[ii,jj] = rr

for col in df.columns: # remove those with zeros
    drop_line = df.index[df[col] == 0].tolist()
    df = df.drop(drop_line)

def find_protein(tag):
    if df['Symbols'].str.contains(tag,case=False).any():
        print('{} found'.format(tag))
    else:
        print('{} not found'.format(tag))
find_protein('S100A10')
find_protein('BASP1')

df.to_csv(df_dir,index=False)
## visualize the data
