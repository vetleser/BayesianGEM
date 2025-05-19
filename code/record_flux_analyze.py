import os
import copy
import time
import math
import pickle
import random
import logging

import numpy as np
import pandas as pd
import numpy.typing as npt
from numpy.typing import NDArray
from typing import Dict, Iterable, List, Optional
from functools import partial
from itertools import starmap

from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import r2_score
from scipy.optimize import minimize

import evo_etc as evo
import evo_etc as CrowdingDE
import GEMS
from random_sampler import RV
from etcpy import etc
from etcpy.thermal_parameters import format_input
import multiprocessing


# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))


outdir = "../results/analysis/flux_analysis"
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")

df_r = load_pickle(f"{outdir}/combined_flux_data_119932.0.pkl")
#filtered_df = df_r[df_r.index.str.match(r"r_\d{4}No1")]

logging.info(f"Loaded df_r: \n {df_r}")
logging.info(f"df_r shape: {df_r.shape}")
logging.info(f"df_r columns: {df_r.columns[:10]}")

