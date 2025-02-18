#!/usr/bin/env python
# coding: utf-8

import copy
import pickle
import pandas as pd
import numpy as np
import logging
import evo_etc as CrowdingDE

from typing import Dict
import evo_etc as evo
import GEMS
import multiprocessing
from functools import partial
import numpy.typing as npt
from random_sampler import RV
import time

import random
from random_sampler import RV
from scipy.optimize import minimize
from numpy.typing import NDArray
import os

import gurobipy as gp


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")
outdir = "../results/analysis"
start_full = time.time()

simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]

task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

default_seed = 200
# Set NumPy global seed
np.random.seed(default_seed)
rng = np.random.Generator(np.random.PCG64(default_seed))

# # Set Python's built-in random module seed
# random.seed(default_seed)

gp.setParam('Seed', default_seed)
gp.setParam('Threads', 1)
gp.setParam('Method', 0)  # Default simplex method (may be more stable)
gp.setParam('Presolve', 0)  # Disable presolve for strict reproducibility
gp.setParam('Crossover', 0)  # Disable crossover for barrier method
np.set_printoptions(precision=15)


#Assert random seeds

# file = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")
# file = file.iloc[:, :-2] #Remove trailing columns of particle ID and simulation number
# model_particle: candidateType = file.loc[file["r2"].idxmax()].to_dict()
# r2_value = -model_particle.pop("r2")
# logging.info(f"r2 value of model particle is: {r2_value}")

# Load the data, create particle as dict
logging.info("Load particle and transform to dict")
file = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")
best_row = file.loc[file["particle_ID"] == 119932.0].iloc[0] #Particle ID of the particle with highest r2 score, found in previous simulations
model_particle: candidateType = best_row.drop(["r2", "particle_ID", "frame_ID"]).to_dict()
r2_value = -best_row["r2"]

logging.info(f"Selected particle ID: {best_row['particle_ID']}, r2 value: {r2_value}")


#Import necessary functions and data
logging.info("Import necessary functions and data")
simulator = GEMS.simulate_at_two_conditions_2
distance_function = GEMS.distance_2

Yobs_batch = GEMS.aerobic_exp_data()
dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}
Yobs = {'rae':Yobs_batch['data'],
            'ran':Yobs_batch_an['data']}

dump_pickle(Yobs, f"{outdir}/Yobs_{task_idx}.pkl")

def evaluate_candidate(param_values: NDArray[np.float64]):
    # Specifying timeout of 30 minutes
    # timeout = 30 * 60
    start = time.time()

    #Reconstructing particle
    candidate: candidateType = dict(zip(all_params, all_params_values))
    for i, param in enumerate(params_to_change):
        candidate[param] = param_values[i]
    
    dump_pickle(candidate, f"{outdir}/model_particle_{task_idx}.pkl")

    success = False
    simulated_data = None
    # No need for parallel processing
    try:
        simulated_data: simResultType = simulator(candidate)
        logging.info("Evaluation of candidate ran successfully")
        success = True
    except Exception as e:
        logging.error(f"Candidate evaluation failed: {e}")
    
    if success == True and simulated_data is not None:
        distance = distance_function(Yobs, simulated_data)
    dump_pickle(simulated_data, f"{outdir}/simulated_data_{task_idx}.pkl")
    #print(simulated_data)
        
        

    end = time.time()
    logging.debug(f'Completed evaluation of candidate in {end - start} seconds')
    logging.info(f"r2: {distance}")
    return distance

def extract_x0(nvar: int):
    cv_df = load_pickle("../results/analysis/cv_df.pkl")
    cv_df_sorted = cv_df.iloc[:, cv_df.iloc[0].argsort()]

    params_to_change = cv_df_sorted.columns[:nvar].tolist()
    logging.info(f"Parameters to change: {params_to_change}")
    x0 = np.array([model_particle.get(param) for param in params_to_change])
    return x0, params_to_change



nvar = 2
maxiter = 2
maxls = 10

x0, params_to_change = extract_x0(nvar=nvar)
all_params = np.array(list(model_particle.keys()))
all_params_values = np.array(list(model_particle.values()))

x0_fixed = copy.deepcopy(x0)

distances = []

# for i in range(3):
#     logging.info(f"Running evaluation {i+1}")
#     d = evaluate_candidate(x0)
#     distances.append(d)
#     logging.info(f"Simulation {i+1}: Distance is {d}, difference is {d-r2_value}")


calculated_r2 = evaluate_candidate(x0_fixed)
#logging.info(f"Calculated r2: {calculated_r2}")



# logging.info("Begin Gradient Search")
# result = minimize(fun=evaluate_candidate, 
#                   x0=x0, 
#                   method="L-BFGS-B", 
#                   options={"maxiter": maxiter,
#                            "disp": True,
#                            "maxls": maxls}
#                   )

# logging.info(result.message)
# logging.info(f"Original value: {x0_fixed}.  Optimal: {result.x}. Difference {result.x-x0_fixed}")
# logging.info(f"Original R2: {r2_value}. Calculated original R2: {calculated_r2}. Optimal R2: {result.fun}. Difference: {result.fun-calculated_r2}. (Negative=better)")

end_full = time.time()

logging.info(f"nvar: {nvar}, maxiter: {maxiter}, total time: {(end_full-start_full)/60} minutes")


#logging.info(f"Distances are {distances}")
# logging.info(f"Difference is {distance-r2_value}")
logging.info("DONE")

