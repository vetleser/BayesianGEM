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
from collections import Counter


import gurobipy as gp

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

n_simulations = 10
tolerance = 1e-8  # Adjust tolerance if needed


model_particle_list = [load_pickle(f"{outdir}/model_particle_{ID}.pkl") for ID in range(n_simulations)]

logging.info("Evaluating particles")
particle_success = True
for key in model_particle_list[0]:
    #logging.info(f"Evaluating {key}")
    reference_value = model_particle_list[0][key]
    
    for i, model in enumerate(model_particle_list[1:], start=1):
        if model[key] != reference_value:
            logging.info(f"Value of parameter {key} different in 0 and {i}")
            particle_success = False
            break

if particle_success:
    logging.info("All particles are identical")
else:
    logging.info("Inconsistencies between particles")


# for key in model_particle1:
#     logging.info(f"Evaluating {key}")
#     if model_particle1[key] != model_particle2[key]:
#         logging.info(f"Value of parameter {key} different in 1 and 2")
#         break
#     if model_particle1[key] != model_particle3[key]:
#         logging.info(f"Value of parameter {key} different in 1 and 3")
#         break

logging.info("Evaluating Yobs")
Yobs_list = [load_pickle(f"{outdir}/Yobs_{ID}.pkl") for ID in range(n_simulations)]
Yobs_success = True
for key in Yobs_list[0]:
    #logging.info(f"Evaluating {key}")

    reference_value = Yobs_list[0].get(key)

    for i, Yobs in enumerate(Yobs_list[1:], start=1):
        current_value = Yobs.get(key)

        if isinstance(reference_value, np.ndarray) and isinstance(current_value, np.ndarray):
            # Find mismatches with tolerance
            diff_indices = np.where(~np.isclose(reference_value, current_value, atol=tolerance, rtol=0))[0]
            if diff_indices.size > 0:
                logging.info(f"Values of parameter {key} differ in 0 and {i} at indices {diff_indices}")
                logging.info(f"Values in 0: {reference_value[diff_indices]}")
                logging.info(f"Values in {i}: {current_value[diff_indices]}")
                Yobs_success = False
        else:
            if not np.isclose(reference_value, current_value, atol=tolerance):
                logging.info(f"Value of parameter {key} different in 0 and {i}: {reference_value} vs {current_value}")
                Yobs_success = False

if Yobs_success:
    logging.info("Yobs are identical")
else:
    logging.info("Yobs are not identical")


simulated_data_list = [load_pickle(f"{outdir}/simulated_data_{ID}.pkl") for ID in range(n_simulations)]

sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]

def evaluate_simulated_data1():
    logging.info("Evaluating simulated data")

    for key in simulated_data_list[0]:
        logging.info(f"Evaluating {key}")

        reference_value = simulated_data_list[0].get(key)

        for i, simulated_data in enumerate(simulated_data_list[1:], start=1):
            current_value = simulated_data.get(key)

            if isinstance(reference_value, np.ndarray) and isinstance(current_value, np.ndarray):
                # Find mismatches with tolerance
                diff_indices = np.where(~np.isclose(reference_value, current_value, atol=tolerance, rtol=0))[0]

                if diff_indices.size > 0:
                    differing_temps = [sel_temp[idx] for idx in diff_indices if idx < len(sel_temp)]
                    
                    logging.info(f"Values of parameter {key} differ in 0 and {i} at temperatures {differing_temps}")
                    logging.info(f"Values in 0: {reference_value[diff_indices]}")
                    logging.info(f"Values in {i}: {current_value[diff_indices]}")
            else:
                if not np.isclose(reference_value, current_value, atol=tolerance):
                    temp = sel_temp[i] if i < len(sel_temp) else f"Index {i}"
                    logging.info(f"Value of parameter {key} different at temperature {temp}: {reference_value} vs {current_value}")


def evaluate_simulated_data2():
    logging.info("Evaluating simulated data")
    temp_diff_summary = {}

    for key in simulated_data_list[0]:
        logging.info(f"Evaluating {key}")

        reference_value = simulated_data_list[0].get(key)

        if isinstance(reference_value, np.ndarray):
            for temp, ref_val in zip(sel_temp, reference_value):
                if temp not in temp_diff_summary:
                    temp_diff_summary[temp] = {}
                if key not in temp_diff_summary[temp]:
                    temp_diff_summary[temp][key] = set()
                temp_diff_summary[temp][key].add(ref_val)
        else:
            temp = sel_temp[0]  # Assign reference value to the first temperature
            if temp not in temp_diff_summary:
                temp_diff_summary[temp] = {}
            if key not in temp_diff_summary[temp]:
                temp_diff_summary[temp][key] = set()
            temp_diff_summary[temp][key].add(reference_value)

        for i, simulated_data in enumerate(simulated_data_list[1:], start=1):
            current_value = simulated_data.get(key)

            if isinstance(reference_value, np.ndarray) and isinstance(current_value, np.ndarray):
                # Find mismatches with tolerance
                diff_indices = np.where(~np.isclose(reference_value, current_value, atol=tolerance, rtol=0))[0]

                if diff_indices.size > 0:
                    differing_temps = [sel_temp[idx] for idx in diff_indices if idx < len(sel_temp)]
                    
                    for temp, ref_val, cur_val in zip(differing_temps, reference_value[diff_indices], current_value[diff_indices]):
                        if temp not in temp_diff_summary:
                            temp_diff_summary[temp] = {}
                        if key not in temp_diff_summary[temp]:
                            temp_diff_summary[temp][key] = []
                        temp_diff_summary[temp][key].add(cur_val)

            else:
                if not np.isclose(reference_value, current_value, atol=tolerance):
                    temp = sel_temp[i] if i < len(sel_temp) else f"Index {i}"
                    if temp not in temp_diff_summary:
                        temp_diff_summary[temp] = {}
                    if key not in temp_diff_summary[temp]:
                        temp_diff_summary[temp][key] = []
                    temp_diff_summary[temp][key].add(current_value)

    # Compute min-max differences
    temp_diff_results = {
        temp: {
            key: max(values) - min(values) if values else 0
            for key, values in key_values.items()
        }
        for temp, key_values in temp_diff_summary.items()
    }

    # Log the final dictionary
    filtered_results = {
        temp: {key: diff for key, diff in key_values.items() if diff > 0}
        for temp, key_values in temp_diff_results.items()
    }
    filtered_results = {temp: key_values for temp, key_values in filtered_results.items() if key_values}  # Remove empty entries

    logging.info(f"Differences by temperature (max-min): {filtered_results}")

def evaluate_simulated_data3():
    logging.info("Evaluating simulated data")
    temp_diff_summary = {}

    for key in simulated_data_list[0]:
        logging.info(f"Evaluating {key}")

        reference_value = simulated_data_list[0].get(key)

        if isinstance(reference_value, np.ndarray):
            for temp, ref_val in zip(sel_temp, reference_value):
                if temp not in temp_diff_summary:
                    temp_diff_summary[temp] = {}
                if key not in temp_diff_summary[temp]:
                    temp_diff_summary[temp][key] = set()
                temp_diff_summary[temp][key].add(ref_val)
        else:
            temp = sel_temp[0]  # Assign reference value to the first temperature
            if temp not in temp_diff_summary:
                temp_diff_summary[temp] = {}
            if key not in temp_diff_summary[temp]:
                temp_diff_summary[temp][key] = set()
            temp_diff_summary[temp][key].add(reference_value)

        for i, simulated_data in enumerate(simulated_data_list[1:], start=1):
            current_value = simulated_data.get(key)

            if isinstance(reference_value, np.ndarray) and isinstance(current_value, np.ndarray):
                # Find mismatches with tolerance
                diff_indices = np.where(~np.isclose(reference_value, current_value, atol=tolerance, rtol=0))[0]

                if diff_indices.size > 0:
                    differing_temps = [sel_temp[idx] for idx in diff_indices if idx < len(sel_temp)]
                    
                    for temp, cur_val in zip(differing_temps, current_value[diff_indices]):
                        if temp not in temp_diff_summary:
                            temp_diff_summary[temp] = {}
                        if key not in temp_diff_summary[temp]:
                            temp_diff_summary[temp][key] = set()
                        temp_diff_summary[temp][key].add(cur_val)

            else:
                if not np.isclose(reference_value, current_value, atol=tolerance):
                    temp = sel_temp[i] if i < len(sel_temp) else f"Index {i}"
                    if temp not in temp_diff_summary:
                        temp_diff_summary[temp] = {}
                    if key not in temp_diff_summary[temp]:
                        temp_diff_summary[temp][key] = set()
                    temp_diff_summary[temp][key].add(current_value)

    # Convert to number of unique values
    temp_diff_results = {
        temp: {key: len(values) for key, values in key_values.items()}
        for temp, key_values in temp_diff_summary.items()
    }

    filtered_results = {
        temp: {key: count for key, count in key_values.items() if count > 1}
        for temp, key_values in temp_diff_results.items()
    }
    filtered_results = {temp: key_values for temp, key_values in filtered_results.items() if key_values}  # Remove empty entries

    logging.info(f"Differences by temperature (unique value count): {filtered_results}")


                            


evaluate_simulated_data1()
evaluate_simulated_data3()
evaluate_simulated_data2()

Yobs = Yobs_list[0]
distance_function = GEMS.distance_2
distances =[]

logging.info("Measuring distances")
for i, simulated_data in enumerate(simulated_data_list):
    d = distance_function(Yobs, simulated_data)
    logging.info(f"R2 for simulation {i} is {d}")

logging.info("Evaluating Distance Function")
distances2 = []
for i in range(10):
    d = distance_function(Yobs_list[0], simulated_data_list[0])
    distances2.append(d)

logging.info(f"Distances are: {distances2}")

logging.info("DONE")