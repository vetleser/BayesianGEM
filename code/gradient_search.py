#!/usr/bin/env python
# coding: utf-8

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

import gurobipy as gp


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")
outdir = "../results/analysis"

simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]


# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

default_seed = 200
# Set NumPy global seed
np.random.seed(default_seed)

# Set Python's built-in random module seed
random.seed(default_seed)

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

# Load the data
file = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")

# Select the row corresponding to particle ID 119932.0
best_row = file.loc[file["particle_ID"] == 119932.0].iloc[0] #Particle ID of the particle with highest r2 score, found in previous simulations

# Extract model parameters, removing unnecessary columns
model_particle: candidateType = best_row.drop(["r2", "particle_ID", "frame_ID"]).to_dict()

# Get the R2 value
r2_value = -best_row["r2"]

logging.info(f"Selected particle ID: {best_row['particle_ID']}, r2 value: {r2_value}")


simulator = GEMS.simulate_at_two_conditions_2
distance_function = GEMS.distance_2

Yobs_batch = GEMS.aerobic_exp_data()
dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}





Yobs = {'rae':Yobs_batch['data'],
            'ran':Yobs_batch_an['data']}

def evaluate_candidate(candidate: candidateType):
    # Specifying timeout of 30 minutes
    # timeout = 30 * 60
    start = time.time()

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
        
        

    end = time.time()
    logging.debug(f'Completed evaluation of candidate in {end - start} seconds')
    return distance

distances = []
for i in range(1):
    logging.info(f"Running evaluation {i+1}")
    d = evaluate_candidate(model_particle)
    distances.append(d)
    logging.info(f"Simulation {i+1}: Distance is {d}, difference is {d-r2_value}")

# logging.info(f"Distance is {distance}")
# logging.info(f"Difference is {distance-r2_value}")
logging.info("DONE")

