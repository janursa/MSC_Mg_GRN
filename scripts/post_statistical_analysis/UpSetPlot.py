"""
Read the results of statistical analysis and extract DE proteins, imputated data, sig data, and output shortlisted data
This outputs DE data for early and late for top p values of 0.05 and 0.025
"""
import sys
import os
import typing
import pandas as pd
import json
import upsetplot
import matplotlib.pyplot as plt
from collections import namedtuple
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from scripts.imports import STATISTICAL_ANALYSIS_DIR, DATA_DIR, F_DE_protiens
from scripts.utils import make_title_pretty


if __name__ == '__main__':
    DE_proteins = F_DE_protiens()
    MySet = namedtuple('MySet', ['name', 'elements'])
    data = [MySet(name=k, elements=set(v)) for k, v in DE_proteins.items()]

    # Create UpSet plot
    upset_data = upsetplot.from_contents(data)
    # merged_list = []
    # [merged_list.extend(value) for value in DE_proteins.values()]
    # merged_list = list(set(merged_list))
    # names_pool = list(DE_proteins.values()
    upsetplot.plot(upset_data)
    # plt.show()
