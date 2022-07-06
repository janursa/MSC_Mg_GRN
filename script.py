#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:40:58 2022

@author: matin
"""
import sys
import os
from pathlib import Path
import pandas as pd
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file)
sys.path.insert(0,main_dir)

## read the data
df_dir = os.path.join(main_dir,'data','omics.xlsx')
if 'df' not in locals():
    print("data is not found")
    original_df = pd.read_excel(df_dir)
## visualize the data
dates = [1,2] # in dayx
prefix = 'LogLFQ intensity '
db = {'Entry name':original_df['Entry name']}

for ii in dates:
    tag_ctr = prefix + str(dates[ii]) + '_0'
    tag_mg = prefix + str(dates[ii]) + '_1'
    db['ctr_{}'.format(ii)] = original_df[tag_ctr]
    db['Mg_{}'.format(ii)] = original_df[tag_mg]
    
# new_db = pd.DataFrame()
# print(df[tag])
