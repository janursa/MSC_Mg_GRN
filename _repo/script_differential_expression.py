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
import matplotlib.pyplot as plt
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file)
sys.path.insert(0,main_dir)

## read the data
df_dir =os.path.join(main_dir,'data','omics.csv')
df = pd.read_csv(df_dir)

cols = df.columns.to_list()

for opt in ['ID','Symbols']:
    cols.remove(opt)

time = []
for col in cols:
    nn = int(col.split('_')[-1])
    if nn not in time:
        time.append(nn)
## --- differential expression with respect to day 1 for both ctr and Mg groups----##

## --- differential expression of Mg with respect to ctr ----###
ee = 0.5
_min = 1
bool_diff = pd.DataFrame()
for ii in time:
    log2FC = df['Mg_{}'.format(ii)] - df['ctr_{}'.format(ii)]
    bool_diff[ii] = (abs(log2FC) > ee)

indices = []
for ii in range(len(bool_diff)):
    try:
        true_count = bool_diff.iloc[ii,:].value_counts()[True]
    except:
        continue
    if true_count >= _min:
        indices.append(ii)
##------------- filter df
df_s = df.iloc[indices,:]
def find_protein(df,tag):
    if df['Symbols'].str.contains(tag,case=False).any():
        print('{} found'.format(tag))
    else:
        print('{} not found'.format(tag))
find_protein(df_s,'S100A10')
find_protein(df_s,'BASP1') 

#------- plot S100A10 anyways
def plot_ctr_vs_mg(df,tag,time):
    linewidth = 2
    stich = lambda _str, ii: '{}_{}'.format(_str,ii)  
    ctrs = [stich('ctr',ii) for ii in time]
    Mgs = [stich('Mg',ii) for ii in time]
    df_tag = df.loc[df['Symbols'] == tag]
    ctrs_data = df_tag[ctrs].iloc[0,:].to_list()
    ee_u = [ii+0.5 for ii in ctrs_data]
    ee_d = [ii-0.5 for ii in ctrs_data]
    Mgs_data = df_tag[Mgs].iloc[0,:].to_list()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(time,ctrs_data,linewidth=linewidth, color = 'b', label = 'Ctr')
    ax.plot(time,Mgs_data,linewidth=linewidth,color = 'g',label = 'Mg')
    ax.plot(time,ee_u,color = 'r',linewidth=linewidth-1, linestyle = '--', label = '0.5 error bounds')
    ax.plot(time,ee_d,linewidth=linewidth-1,linestyle = '--',color = 'r')
    ax.set_ylim([6,10])
    ax.set_xticks(time)
    ax.set_title(tag)
    ax.set_xlabel('Days')
    ax.set_ylabel('Intensities (log2)')
    plt.legend()
    
    
plot_ctr_vs_mg(df=df,tag = 'S100A10',time=time)
    