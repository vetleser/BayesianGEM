#!/usr/bin/env python
# coding: utf-8

import copy
import pickle
import numpy as np
import logging

from typing import Dict
import GEMS
from functools import partial
import numpy.typing as npt
from random_sampler import RV
import time
from scipy.optimize import minimize
from numpy.typing import NDArray



logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")
outdir = "../results/analysis"
start_full = time.time()

simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]


# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

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


def evaluate_candidate(param_values: NDArray[np.float64]):
    # Specifying timeout of 30 minutes
    # timeout = 30 * 60
    logging.info(f"Current param values to evaluate: {param_values}")
    start = time.time()
    iteration.append(1)
    logging.info(f"Current iteration is {len(iteration)-1}")
    #Reconstructing particle, takes approximately 0.001 seconds
    candidate: candidateType = dict(zip(all_params, all_params_values))
    for i, param in enumerate(params_to_change):
        candidate[param] = param_values[i]

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
        logging.info(f"Simulated data is {simulated_data}")
        
        

    end = time.time()
    logging.debug(f'Completed evaluation of candidate in {end - start} seconds')
    logging.info(f"r2: {distance}")
    return distance

def extract_x0_and_params(nvar: int):
    cv_df = load_pickle("../results/analysis/cv_df.pkl")
    cv_df_sorted = cv_df.iloc[:, cv_df.iloc[0].argsort()]

    params_to_change = cv_df_sorted.columns[:nvar].tolist()
    logging.info(f"Parameters to change: {params_to_change}")
    x0 = np.array([model_particle.get(param) for param in params_to_change])
    return x0, params_to_change


iteration = []

nvar = 100
maxiter = 10
maxls = 20

x0, params_to_change = extract_x0_and_params(nvar=nvar)
all_params = np.array(list(model_particle.keys()))
all_params_values = np.array(list(model_particle.values()))

x0_fixed = copy.deepcopy(x0)


calculated_r2 = evaluate_candidate(x0_fixed)
#logging.info(f"Calculated r2: {calculated_r2}")



logging.info("Begin Gradient Search")
result = minimize(fun=evaluate_candidate, 
                  x0=x0, 
                  method="L-BFGS-B", 
                  options={"maxiter": maxiter,
                           "disp": True,
                           "maxls": maxls}
                  )

logging.info(result.message)
logging.info(f"Original value: {x0_fixed}.  Optimal: {result.x}. Difference {result.x-x0_fixed}")
logging.info(f"Original R2: {r2_value}. Calculated original R2: {calculated_r2}. Optimal R2: {result.fun}. Difference: {result.fun-calculated_r2}. (Negative=better)")

end_full = time.time()

logging.info(f"nvar: {nvar}, maxiter: {maxiter}, maxls: {maxls}, total time: {(end_full-start_full)/60} minutes")


#logging.info(f"Distances are {distances}")
# logging.info(f"Difference is {distance-r2_value}")
logging.info("DONE")

