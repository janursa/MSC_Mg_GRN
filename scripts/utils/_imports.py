import os
import pathlib
import sys

import pandas as pd
import numpy as np
pd.options.mode.chained_assignment

import typing
import random
import json
import shutil
import scipy
import statistics
import copy

from sklearn import preprocessing
from sklearn import ensemble
from sklearn import inspection
from scipy import stats

import matplotlib
import matplotlib.font_manager as font_manager
import matplotlib.markers as mmarkers
import matplotlib.pyplot as plt
# shutil.rmtree(matplotlib.get_cachedir())

MAIN_DIR = os.path.join(pathlib.Path(__file__).parent.resolve(), '../..')
sys.path.insert(0, MAIN_DIR)

#-- import from geneRNI
geneRNI_dir = os.path.join(MAIN_DIR,'..','geneRNI')
sys.path.insert(0, geneRNI_dir)
from geneRNI import geneRNI, tools, search_param
from geneRNI.data import Data
from geneRNI.models import get_estimator_wrapper
# print(pd)
