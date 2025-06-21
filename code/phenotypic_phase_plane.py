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
from typing import Iterable, List, Optional
import numpy as np
import pandas as pd
import time
from reframed import CBModel
from reframed.solvers.solver import Solver
import reframed
import logging

import etcpy.reframed_mappers as reframed_mappers
import etcpy.thermal_parameters as thermal_parameters
from .thermal_parameters import calculate_thermal_params

from sympy import Float
import gurobipy as gp

import gurobipy as gp

class OptimizationError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


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

def simulate_growth(model: CBModel, Ts,sigma,param_dict,Tadj=0, max_attempts = 1):
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
            # if(attempt>1):
            #     gp.setParam('Seed', attempt)
            try:
                solution = solver.solve(linear=model.get_objective(),minimize=False)
                if solution.status != reframed.solvers.solution.Status.OPTIMAL:
                    raise OptimizationError(f"Solver status is {solution.status.value}")
                r = solution.fobj
                #logging.info(f"Model solved successfully at temperature {T} at attempt {attempt}")
                success = True
                break
            except OptimizationError as err:
                if attempt == max_attempts:
                    logging.info(f'Attempt {attempt} failed to solve the problem, problem: {str(err)}. At Temperature {T}')
                pass
                #logging.info(f'Attempt {attempt} failed to solve the problem, problem: {str(err)}. At Temperature {T}')
        if success:
            rs.append(r)
        else:
            #logging.info(f"Failed to solve problem after {max_attempts} attempts at temperature {T}")
            rs.append(0) #Still returns 0 after failing to solve. Should fix later. For example: Return NaN, and stop checking if NaN is encountered in distance function
    return rs

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


oxygen_range = np.linspace(-20, 0, 20)  # uptake fluxes (negative)
glucose_range = np.linspace(-10, 0, 20)
temperature = 310  # Example temperature in K (can be a list too)

growth_rates = np.empty((len(oxygen_range), len(glucose_range)))

for i, o2 in enumerate(oxygen_range):
    for j, glc in enumerate(glucose_range):
        with model:
            # Set bounds for oxygen and glucose uptake
            model.reactions["r_1992"].lower_bound = o2
            model.reactions["r_1992"].upper_bound = 1000
            model.reactions["r_1714"].lower_bound = glc
            model.reactions["r_1714"].upper_bound = 1000

            # Run growth simulation at specified temperature(s)
            # simulate_growth expects a list of temperatures, so [temperature]
            r = simulate_growth(model, Ts=[temperature], sigma=sigma, param_dict=param_dict, Tadj=0, max_attempts=1)
            # r is a list with one element (growth at that temperature)
            growth_rates[i, j] = r[0]

# Plot
plt.figure(figsize=(8,6))
plt.imshow(growth_rates, origin='lower', aspect='auto',
           extent=[glucose_range.min(), glucose_range.max(), oxygen_range.min(), oxygen_range.max()],
           cmap='viridis')
plt.colorbar(label='Growth Rate')
plt.xlabel('Glucose uptake (mmol/gDW/hr)')
plt.ylabel('Oxygen uptake (mmol/gDW/hr)')
plt.title(f'Phenotype Phase Plane at T={temperature} K')
plt.show()
