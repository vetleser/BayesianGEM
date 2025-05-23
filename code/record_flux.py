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

default_seed = 200
# Set NumPy global seed
np.random.seed(default_seed)
rng = np.random.Generator(np.random.PCG64(default_seed))




#Import necessary functions and data
logging.info("Import necessary functions and data")

def convert_to_dataframe(dictionary, filename: str):
    df = pd.DataFrame.from_dict(dictionary, orient='index', columns=['Importance'])
    dump_pickle(df, f"../results/analysis/flux_analysis/{filename}.pkl")
    


def simulate_at_two_conditions_2(args):
    ae_output = aerobic(args)
    data_batch = ae_output['data']
    reac_importance_ae = ae_output['reac_importance']



    an_output = anaerobic_reduced(args)
    data_batch_an= an_output['data']
    reac_importance_an = an_output['reac_importance']
    reac_importance_tot = {key: 0.0 for key in reac_importance_ae}

        # Combine flux dicts with metadata
    combined_data = []
    
    for temp, flux in ae_output['flux_dict_by_temp'].items():
        combined_data.append({**flux})

    for temp, flux in an_output['flux_dict_by_temp'].items():
        combined_data.append({**flux})
    
    # # Create a DataFrame from combined_data
    # df_flux = pd.DataFrame(combined_data)
    # # Save the DataFrame to a pickle file
    # dump_pickle(df_flux, f"../results/analysis/flux_analysis/combined_flux_data_{particle_id_str}.pkl")

    for reaction in reac_importance_ae:
        reac_importance_tot[reaction] = reac_importance_ae[reaction] + reac_importance_an[reaction]
        #if reac_importance_tot[reaction] == 0:
        #    reac_importance_tot.pop(reaction)
    
    # # Assuming reac_importance_ae and reac_importance_an are Dict[str, float]
    # reac_importance = {key: reac_importance_ae.get(key, 0) + reac_importance_an.get(key, 0) 
    #                    for key in set(reac_importance_ae) | set(reac_importance_an)}
    
    return combined_data#'reac_importance_tot': reac_importance_tot}# , reac_importance

def simulate_growth(model: CBModel, Ts,sigma,param_dict,Tadj=0, max_attempts = 10):
    '''
    # model, reframed model
    # Ts, a list of temperatures in K
    # sigma, enzyme saturation factor
    # param_dict, a dictionary containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # working_model, If provided, warm-start FBA will be used for accelerating computations. The working model will be modified during this
    # process
    # Ensure that Topt is in K. Other parameters are in standard units.
    # Tadj, as descrbed in map_fNT
    #
    '''
    rs = list()
    fluxes = list()
    solver: reframed.solvers.GurobiSolver = reframed.solver_instance(model)
    #solver: reframed.solvers.CplexSolver = reframed.solver_instance(model)
    #logging.info(f"Using solver: {type(solver).__name__}")
    #gp.setParam('Seed', 0)
    for T in Ts:
        # map temperature constraints
        mappers = reframed_mappers
        mappers.map_fNT(model,T,param_dict,solver_instance=solver)
        mappers.map_kcatT(model,T,param_dict,solver_instance=solver)
        mappers.set_NGAMT(solver,T)
        mappers.set_sigma(solver,sigma)
        solver.update()
        success = False
        for attempt in range(1, max_attempts+1):
            #logging.info(f"Attempt {attempt} to solve the problem at temperature {T}")
            # if(attempt>1):
            #     gp.setParam('Seed', attempt)
            try:
                solution = solver.solve(linear=model.get_objective(),minimize=False)
                if solution.status != reframed.solvers.solution.Status.OPTIMAL:
                    raise etc.OptimizationError(f"Solver status is {solution.status.value}")
                r = solution.fobj
                flux = solution.values
                logging.info(f"Model solved successfully at temperature {T} at attempt {attempt}")
                success = True
                break
            except etc.OptimizationError as err:
                logging.info(f'Attempt {attempt} failed to solve the problem, problem: {str(err)}. At Temperature {T}')
        if success:
            rs.append(r)
            fluxes.append(flux)
            #logging.info(flux)
        else:
            #logging.info(f"Failed to solve problem after {max_attempts} attempts at temperature {T}")
            
            rs.append(0) #Still returns 0 after failing to solve. Should fix later. For example: Return NaN, and stop checking if NaN is encountered in distance function
            #fluxes.append(None)
    return rs, fluxes



def aerobic(thermalParams):
    # thermalParams: a dictionary with ids like uniprotid_Topt 
    param_dict = format_input(params,thermalParams)
    mae = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
    rae, fluxes = simulate_growth(mae,dfae_batch.index+273.15,param_dict=param_dict,sigma=0.5)
    logging.info("Simulated aerobic growth finished")

    flux_dict_by_temp = {}
    temps = dfae_batch.index + 273.15
    for temp, flux in zip(temps, fluxes):
        flux_dict_by_temp[temp] = flux

    reac_importance : Dict[str, float] = {key: 0.0 for key in fluxes[0]}
    for flux in fluxes:
        for key in flux:
            if flux[key] > 0:
                reac_importance[key] += 1/len(fluxes)
    r_reactions_ae = {k: v for k, v in reac_importance.items() if k.startswith("r_") and not k.endswith("_REV")}

    # Log the sorted dictionary
    logging.info(f"Length of reac_importance_ae: {len(reac_importance)}")
    # logging.info(f"Length of r_reactions_ae: {len(r_reactions_ae)}")
    # logging.info(f"r_reactions_ae: \n {r_reactions_ae}")
    # logging.info(f"sorted_reac_importance: {sorted_reac_importance}")
    
    rae = [0 if x is None else x for x in rae]
    rae = [0 if x<1e-3 else x for x in rae]
    logging.info(f"rae: {rae}")

    if any(np.isnan(x) for x in rae):
        logging.info("NaN in rae GEMS.aerobic")
        return {'data':np.array(rae)}

    rexp = GEMS.aerobic_exp_data()['data']
    
    logging.info(f'r2_batch_ae: {r2_score(rexp,rae)}')
    logging.info(f'MSE_ae: {MSE(rexp,rae)}')
    
    return {'data':np.array(rae), 'reac_importance': reac_importance, 'flux_dict_by_temp': flux_dict_by_temp} #, reac_importance

def anaerobic_reduced(thermalParams):
    param_dict = format_input(params,thermalParams)
    man = pickle.load(open(os.path.join(path,'models/anaerobic.pkl'),'rb'))
    
    sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
    ran, fluxes = simulate_growth(man,np.array(sel_temp)+273.15,param_dict=param_dict,sigma=0.5)
    logging.info("Simulated anaerobic growth finished")

    flux_dict_by_temp = {}
    temps = np.array(sel_temp)+273.15
    for temp, flux in zip(temps, fluxes):
        flux_dict_by_temp[temp] = flux

    reac_importance : Dict[str, float] = {key: 0.0 for key in fluxes[0]}
    for flux in fluxes:
        for key in flux:
            if flux[key] > 0:
                reac_importance[key] += 1/len(fluxes)
    sorted_reac_importance = dict(sorted(reac_importance.items(), key=lambda item: item[1], reverse=True))

    # Log the sorted dictionary
    logging.info(f"Length of reac_importance_an: {len(reac_importance)}")
    # logging.info(f"sorted_reac_importance: {sorted_reac_importance}")
    ran = [0 if x is None else x for x in ran]
    logging.info(f"ran: {ran}")
    rexp = dfan_batch.loc[sel_temp,'r_an'].values
    #anaerobic_exp_data()['data']
    
    
    logging.info(f'r2_batch_an: {r2_score(rexp,ran)}')
    logging.info(f'MSE_an: {MSE(rexp,ran)}')
    logging.info(f'Model error: {len(rexp)} {len(ran)}')

    return  {'data':np.array(ran), 'reac_importance':reac_importance, 'flux_dict_by_temp': flux_dict_by_temp} #, reac_importance



Yobs_batch = GEMS.aerobic_exp_data()
dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}
Yobs = {'rae':Yobs_batch['data'],
            'ran':Yobs_batch_an['data']}

distance_function = GEMS.distance_2
simulator = simulate_at_two_conditions_2




def evaluate_candidate(candidate: candidateType):
    # Specifying timeout of 30 minutes
    # timeout = 30 * 60
    start = time.time()

    
    #dump_pickle(candidate, f"{outdir}/model_particle_{task_idx}.pkl")

    success = False
    #simulated_data = None
    # No need for parallel processing
    try:
        combined_data = simulator(candidate)
        logging.info("Evaluation of candidate ran successfully")
        success = True
    except Exception as e:
        logging.error(f"Candidate evaluation failed: {e}")
    
    #logging.info(f"Simulated data:\n {simulated_data}")

    #reac_importance_tot = simulated_data["reac_importance_tot"]
    #simulated_data= {'rae': simulated_data['rae'], 'ran': simulated_data['ran']}
    # logging.info(f"reac_importance_tot length: {len(reac_importance_tot)}")
    # #logging.info(reac_importance_tot)
    # r_REV_reactions = {k: v for k, v in reac_importance_tot.items() if k.startswith("r_")}
    # r_reactions = {k: v for k, v in reac_importance_tot.items() if k.startswith("r_") and not k.endswith("_REV")}


    
    #distance = distance_function(Yobs, simulated_data)

    #dump_pickle(simulated_data, f"{outdir}/simulated_data_{task_idx}.pkl")
    #print(simulated_data)
        
        

    end = time.time()
    logging.debug(f'Completed evaluation of candidate in {end - start} seconds')
    #logging.info(f"r2: {distance}")
    return combined_data

# Load the data, create particle as dict
# logging.info("Load particle and transform to dict")
# file = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")
# best_row = file.loc[file["particle_ID"] == 119932.0].iloc[0] #Particle ID of the particle with highest r2 score, found in previous simulations
# model_particle: candidateType = best_row.drop(["r2", "particle_ID", "frame_ID"]).to_dict()
# r2_value = -best_row["r2"]

logging.info("Load particle and transform to dict")
n_particles = 10
file = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")
param_columns = [col for col in file.columns if col not in ["particle_ID", "frame_ID", "r2"]]

# Sort by r2, drop duplicates based on parameter values
file_unique_particles = file.sort_values("r2", ascending=False).drop_duplicates(subset=param_columns, keep="first")


logging.info(f"Loaded file: {file_unique_particles}")
best_rows = file_unique_particles.nlargest(n_particles, "r2")
#model_particle: candidateType = best_row.to_dict()
particles_to_evaluate = []
for _, best_row in best_rows.iterrows():
    particle_dict = best_row.to_dict()
    particles_to_evaluate.append(particle_dict)

#logging.info(f"Selected particle ID: {best_row['particle_ID']}")
p_str_to_evaluate = []
counter = 0
for p in particles_to_evaluate:
    counter += 1
    particle_ID = p.pop("particle_ID")
    p_str_to_evaluate.append(particle_ID)
    frame_ID = p.pop("frame_ID")
    r2_value = p.pop("r2")

    logging.info(f"Evaluating candidate: {particle_ID}, {counter} of {n_particles}")
    combined_data = evaluate_candidate(p)
    df_flux = pd.DataFrame(combined_data)
    
    dump_pickle(df_flux, f"../results/analysis/flux_analysis/combined_flux_data_{particle_ID}.pkl")
        

        
    # Evaluate the candidate

dump_pickle(p_str_to_evaluate, f"../results/analysis/flux_analysis/particle_IDs_to_evaluate.pkl")
logging.info(f"Flux data:\n {df_flux}")


logging.info("DONE")