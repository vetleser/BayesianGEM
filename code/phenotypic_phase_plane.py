#!/usr/bin/env python
# coding: utf-8

import os
import copy
import time
import pickle
import matplotlib.pyplot as plt
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

from reframed import CBModel
from reframed.solvers.solver import Solver
import reframed
import logging

import etcpy.reframed_mappers as reframed_mappers
import etcpy.thermal_parameters as thermal_parameters
#from .thermal_parameters import calculate_thermal_params
import cobra
from cobra.flux_analysis.phenotype_phase_plane import production_envelope

import gurobipy as gp


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
logging.info("BEGIN")
outdir = "../results/analysis"
start_full = time.time()


# task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
dfae_batch,dfan_batch = GEMS.load_exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))



# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))





# Load or define your reframed model
mae = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
man = pickle.load(open(os.path.join(path,'models/anaerobic.pkl'),'rb'))
logging.info("Models loaded")
logging.info(f"MAE reactions\n: {len(mae.reactions)}")

# Search for exchange reactions involving glucose or oxygen
for rxn_id, rxn in mae.reactions.items():
    if 'growth' in rxn.name.lower():
        logging.info(f"Reaction {rxn_id}: {rxn.name} (ID: {rxn.id})")
        logging.info(f"Stoichiometry: {rxn}")



oxygen_exchange = "r_1992"
glucose_exchange = "r_1714"
growth_reaction = "r_2111"  


model = mae.copy()

# Set growth reaction as objective
model.objective = "r_2111"

# Define ranges (note: uptake fluxes are negative)
glucose_range = np.linspace(-10, 0, 20)
oxygen_range = np.linspace(-20, 0, 20)

growth_rates = np.empty((len(oxygen_range), len(glucose_range)))

# Loop over oxygen and glucose uptake values
for i, o2 in enumerate(oxygen_range):
    for j, glc in enumerate(glucose_range):
        # Set bounds for oxygen and glucose uptake reactions
        model.reactions["r_1992"].lower_bound = o2
        model.reactions["r_1992"].upper_bound = 1000  # large upper bound
        
        model.reactions["r_1714"].lower_bound = glc
        model.reactions["r_1714"].upper_bound = 1000
        
        # Run FBA
        solution = model.solve()
        
        if solution.status == "optimal":
            growth_rates[i, j] = solution.objective_value
        else:
            growth_rates[i, j] = np.nan

# Plotting heatmap
plt.figure(figsize=(8,6))
plt.imshow(growth_rates, origin='lower', aspect='auto', 
           extent=[glucose_range.min(), glucose_range.max(), oxygen_range.min(), oxygen_range.max()],
           cmap='viridis')
plt.colorbar(label='Growth Rate (1/hr)')
plt.xlabel('Glucose uptake (mmol/gDW/hr)')
plt.ylabel('Oxygen uptake (mmol/gDW/hr)')
plt.title('Phenotype Phase Plane: Growth vs Glucose and Oxygen')
plt.show()

