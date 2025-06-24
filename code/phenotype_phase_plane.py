#!/usr/bin/env python
# coding: utf-8

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
import matplotlib.pyplot as plt

from reframed import CBModel
from reframed.solvers.solver import Solver
import reframed
import logging

import etcpy.reframed_mappers as reframed_mappers
import etcpy.thermal_parameters as thermal_parameters
#from .thermal_parameters import calculate_thermal_params

import gurobipy as gp


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")
outdir = "../results/analysis"
start_full = time.time()

simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]

# task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
dfae_batch,dfan_batch = GEMS.load_exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))



# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

#----------------------------------------------------------------------------------------------------
def simulate_growth_at_conditions(model, T, param_dict, max_attempts=10, sigma=0.5):
    solver = reframed.solver_instance(model)
    reframed_mappers.map_fNT(model, T, param_dict, solver_instance=solver)
    reframed_mappers.map_kcatT(model, T, param_dict, solver_instance=solver)
    reframed_mappers.set_NGAMT(solver, T)
    reframed_mappers.set_sigma(solver, sigma)
    solver.update()
    for attempt in range(max_attempts):
        try:
            solution = solver.solve(linear=model.get_objective(), minimize=False)
            if solution.status == reframed.solvers.solution.Status.OPTIMAL:
                return solution.fobj, solution.values
        except etc.OptimizationError:
            logging.warning(f"Attempt {attempt + 1} failed to solve the problem at temperature {T}. Retrying...")
            continue
    logging.error(f"Failed to solve the problem after {max_attempts} attempts at temperature {T}")
    return 0, None  # failure fallback

logging.info("Load particle and transform to dict")
n_particles = 1
file = load_pickle("../results/sa/sa_combined_df_R090_final.pkl")

# Columns representing parameters (exclude metadata)
param_columns = [col for col in file.columns if col not in {"particle_ID", "frame_ID", "r2"}]

# Sort by r2 descending, drop duplicates based on parameter columns, keep first (best)
file_unique = file.sort_values("r2", ascending=False).drop_duplicates(subset=param_columns, keep="first")

logging.info(f"Loaded unique particles: {len(file_unique)}")

# Get top n_particles best rows
best_rows = file_unique.nlargest(n_particles, "r2")

# Extract parameter dict from first best particle
thermalParams = best_rows.iloc[0][param_columns].to_dict()

# Then build the grid
shikimate_kinase = 'r_0997No1'
temps = np.linspace(280, 320, 5)
sigmas = np.linspace(0.1, 1.0, 5)
tm_values = np.linspace(312, 346, 25)  # Example Tm values in K
topt_values = np.linspace(272, 328, 25)  # Example Tm values in K
growth_grid = np.zeros((len(tm_values), len(topt_values)))
shikimate_grid = np.zeros((len(tm_values), len(topt_values)))
param_dict = format_input(params,thermalParams)
model = pickle.load(open(os.path.join(path, 'models/aerobic.pkl'), 'rb'))
enzyme_id = "P08566"
for i, tm in enumerate(tm_values):
    for j, topt in enumerate(topt_values):
        if tm < topt:
            logging.warning(f"Skipping invalid combination: Tm={tm}, Topt={topt}")
            continue
        T = temps[3]
        sigma = 0.5
        new_thermalParams = copy.deepcopy(thermalParams)
        new_thermalParams[f"{enzyme_id}_Tm"] = tm
        new_thermalParams[f"{enzyme_id}_Topt"] = topt
        param_dict = format_input(params,new_thermalParams)
        #logging.info(f"Simulating growth at Tm={tm}, Topt={topt}")
        growth, values = simulate_growth_at_conditions(model, T, param_dict=param_dict)
        shikimate_flux = values[shikimate_kinase]  # Example reaction ID for oxygen uptake
        logging.info(f"Tm={tm}, Topt={topt}:Shikimate flux {shikimate_flux}, growth {growth}")
        growth_grid[i, j] = growth
        shikimate_grid[i, j] = shikimate_flux

# Plotting as above


o2_uptake_range = np.linspace(0, 20, 5)  # in mmol/gDW/h
c_uptake_range = np.linspace(0, 20, 5)
o2_rxn_id = 'r_1992'  # example reaction ID for oxygen uptake
c_rxn_id = 'r_1714'  # example reaction ID for glucose (carbon source) uptake


dump_pickle(growth_grid, "../results/analysis/growth_grid.pkl")
dump_pickle(shikimate_grid, "../results/analysis/shikimate_grid.pkl")


plt.figure(figsize=(8,6))
X, Y = np.meshgrid(topt_values, tm_values)
cp = plt.contourf(X, Y, growth_grid, cmap='viridis')
plt.colorbar(cp, label='Growth rate')
plt.plot(X, X, color='white', linestyle='--', label='Tm = Topt')
plt.legend()
plt.xlabel('Enzyme optimum temperature (Topt, K)')
plt.ylabel('Enzyme melting temperature (Tm, K)')
plt.title('Phenotype Phase Plane: Growth vs Enzyme parameters')
plt.savefig(f"../figures/phase_plane_growth.png", dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(8,6))
cp = plt.contourf(X, Y, shikimate_grid, cmap='viridis')
plt.colorbar(cp, label='Shikimate kinase flux')
plt.plot(X, X, color='white', linestyle='--', label='Tm = Topt')
plt.legend()
plt.xlabel('Enzyme optimum temperature (Topt, K)')
plt.ylabel('Enzyme melting temperature (Tm, K)')
plt.title('Phenotype Phase Plane: Shikimate Flux vs Enzyme parameters')
plt.savefig(f"../figures/phase_plane_shikimate_flux.png", dpi=300, bbox_inches='tight')

logging.info("DONE")