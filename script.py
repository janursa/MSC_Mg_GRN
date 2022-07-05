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

data_dir = os.path.join(main_dir,'data','omics.xlsx')

data = pd.read_excel(data_dir)