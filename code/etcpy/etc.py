from typing import List
import numpy as np
import pandas as pd
import time
from reframed import CBModel
from reframed.solvers.solver import Solver
import reframed
import logging

import etcpy.reframed_mappers as reframed_mappers
from .thermal_parameters import calculate_thermal_params

from sympy import Float

T0 = 273.15
SLACK_FACTOR = 1.0      

class OptimizationError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)



def simulate_growth(model: CBModel, Ts,sigma,param_dict,Tadj=0):
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
    for T in Ts:
        # map temperature constraints
        mappers = reframed_mappers
        mappers.map_fNT(model,T,param_dict,solver_instance=solver)
        mappers.map_kcatT(model,T,param_dict,solver_instance=solver)
        mappers.set_NGAMT(solver,T)
        mappers.set_sigma(solver,sigma)
        solver.update()
        try:
            solution = solver.solve(linear=model.get_objective(),minimize=False)
            if solution.status != reframed.solvers.solution.Status.OPTIMAL:
                raise OptimizationError(f"Solver status is {solution.status.value}")
            r = solution.fobj
            logging.info("Model solved successfully")
        except OptimizationError as err:
            logging.info(f'Failed to solve the problem, problem: {str(err)}')
            r = 0
        print(T-273.15,r)
        rs.append(r)
    return rs


def simulate_chemostat(model: CBModel,dilu,param_dict,Ts,sigma,growth_id,glc_up_id,prot_pool_id):
    '''
    # Do simulation on a given dilution and a list of temperatures. 
    # model: reframed model
    # dilu, dilution rate
    # param_dict, a dictionary containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ts, a list of temperatures to simulate at. in K
    # sigma, saturation factor
    # growth_id, reaction of of growth
    # glc_up_id, reaction id of glucose uptake reaction
    # prot_pool_id, reaction id of prot_pool_exchange. --> prot_pool
    # working_model, If provided, warm-start FBA will be used for accelerating computations. The working model will be modified during this
    # process

    '''
    solutions = list() # corresponding to Ts. a list of solutions from model.optimize()
    mappers = reframed_mappers
    # Step 1: fix growth rate
    m0 = model.copy()
    m0.reactions[growth_id].lb = dilu
    solver: reframed.solvers.GurobiSolver = reframed.solver_instance(m0)
    for T in Ts:
        # Step 2: map temperature constraints. 
        mappers.map_fNT(m0,T,param_dict,solver_instance=solver)
        mappers.map_kcatT(m0,T,param_dict,solver_instance=solver)
        mappers.set_NGAMT(solver,T)
        mappers.set_sigma(solver,sigma)
        solver.update()
        try:
            # Step 3: set objective function as minimizing glucose uptake rate and protein useage
            glc_flux_solution = solver.solve(linear={glc_up_id: 1}, minimize=True)
            if glc_flux_solution.status != reframed.solvers.solution.Status.OPTIMAL:
                raise OptimizationError(f"Solver status is {glc_flux_solution.status.value}")
            glc_min_flux = glc_flux_solution.fobj
            solver.add_variable(var_id=glc_up_id,lb=0,ub=glc_min_flux*SLACK_FACTOR,update=True)
            solution = solver.solve(linear= {prot_pool_id: 1},minimize=True)
            if solution.status != reframed.solvers.solution.Status.OPTIMAL:
                raise OptimizationError(f"Solver status is {solution.status.value}")
            solutions.append(solution.values)
            logging.info('Model solved successfully')
        except OptimizationError as err:
            logging.info(f'Failed to solve the problem, problem: {str(err)}')
            #solutions.append(None)
            break # because model has been impaired. Further simulation won't give right output.
        finally:
            # Reset glucose flux bound
            glc_default_bound = model.reactions[glc_up_id].ub
            solver.add_variable(var_id=glc_up_id,lb=0,ub=glc_default_bound,update=True)
    return solutions

def sample_data_uncertainty(params,columns=None):
    '''
    # params is a dataframe with following columns:
    # Tm,Tm_std:  melting temperature. Given in K
    #
    # T90, temperature at which 90% of an enzyme is denatured. This is not mandentory. If missing, protein length will be used
    #      to calculate denaturation curve. Given in K
    # 
    # dCpt,dcpt_std: the heat capacity difference between transition state and ground state in enzyme catalyzed reaction. Given
    #                in J/mol/K
    # 
    # Topt,Topt_std: the optimal temprature at which the specific activity is maximized. Given in K
    #
    # Length, protein length. This will not be used in this function, but it will be used later in the calculation of thermal
    #         parameters.
    #
    # xx_std, corresponding uncertainty given by standard deviation.
    # 
    # columns: a list of columns to be sampled, could be any combination of ['Tm','dCpt','Topt]. 
    #          If it is None, then sample all three columns
    # 
    # The script will return an new dataframe with the same columns but with randomly sampled data
    '''
    sampled_params = params.copy()
    if columns is None: columns = ['Tm','dCpt','Topt']
    for col in columns:
        for ind in params.index: 
            sampled_params.loc[ind,col] = np.random.normal(params.loc[ind,col],params.loc[ind,col+'_std'])
            if col == 'Tm': 
                sampled_params.loc[ind,'T90'] = sampled_params.loc[ind,col] + params.loc[ind,'T90']-params.loc[ind,col]
    return sampled_params


def sample_data_uncertainty_with_constraint(inpt,columns=None):
    if type(inpt)==tuple:
        params,seed = inpt
        np.random.seed(seed+int(time.time()))
    else: params = inpt
    '''
    # params is a dataframe with following columns:
    # Tm,Tm_std:  melting temperature. Given in K
    #
    # T90, temperature at which 90% of an enzyme is denatured. This is not mandentory. If missing, protein length will be used
    #      to calculate denaturation curve. Given in K
    # 
    # dCpt,dcpt_std: the heat capacity difference between transition state and ground state in enzyme catalyzed reaction. Given
    #                in J/mol/K
    # 
    # Topt,Topt_std: the optimal temprature at which the specific activity is maximized. Given in K
    #
    # Length, protein length. This will not be used in this function, but it will be used later in the calculation of thermal
    #         parameters.
    #
    # xx_std, corresponding uncertainty given by standard deviation.
    # 
    # columns: a list of columns to be sampled, could be any combination of ['Tm','dCpt','Topt]. 
    #          If it is None, then sample all three columns
    # 
    # The script will return an new dataframe with the same columns but with randomly sampled data
    '''
    
    sampled_params = params.copy()
    if columns is None: columns = ['Tm','dCpt','Topt']
    for col in columns:
        lst = [np.random.normal(params.loc[ind,col],params.loc[ind,col+'_std']) for ind in sampled_params.index]
        sampled_params[col] = lst
          
    # resample those ones with Topt>=Tm
    for ind in sampled_params.index:
        tm,topt = sampled_params.loc[ind,'Tm'],sampled_params.loc[ind,'Topt']
        count = 0
        while topt>=tm:
            count += 1
            if 'Topt' in columns:
                topt = np.random.normal(params.loc[ind,'Topt'],params.loc[ind,'Topt_std'])
            if 'Tm' in columns: 
                tm = np.random.normal(params.loc[ind,'Tm'],params.loc[ind,'Tm_std'])
            if 'Topt' not in columns and 'Tm' not in columns:
                break
            if count>10: break
        sampled_params.loc[ind,'Tm'],sampled_params.loc[ind,'Topt'] = tm,topt
    
    # update Topt
    sampled_params['T90'] = sampled_params['Tm']+params['T90']-params['Tm']
    return sampled_params


def simulate_fva(model: CBModel,Ts: List[float],sigma: float,param_dict: dict,Tadj=0) -> pd.DataFrame:
    '''
    Simulate FVA on growth scenarios
    # model, reframed model
    # Ts, a list of temperatures in K
    # sigma, enzyme saturation factor
    # param_dict, a dictionary containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ensure that Topt is in K. Other parameters are in standard units.
    # Tadj, as descrbed in map_fNT
    # 
    '''
    rs: List[pd.DataFrame] = list()
    # This wierd DataFrame with zero columns is intended to prevent nasty errors when all
    # FVA solutions fail
    skeleton_df = pd.DataFrame({'reaction': [],'T': [], 'minimum': [], 'maximum': []})
    rs.append(skeleton_df)
    mappers = reframed_mappers
    for T in Ts:
        opt_model = model.copy()
        # map temperature constraints
        mappers.map_fNT(opt_model,T,param_dict)
        mappers.map_kcatT(opt_model,T,param_dict)
        mappers.set_NGAMT(opt_model,T)
        mappers.set_sigma(opt_model,sigma)
        try:
            # We have to take into consideration that the model might be infeasible as-is and
            # reframed has no proper way to handling this case, causing a TypeError (!) in such cases
            fva_solution = reframed.FVA(opt_model, obj_frac=0.999)
        except Exception as err:
            logging.info(f'Failed to solve the problem, problem: {str(err)}')
            res_frame = pd.DataFrame({'reaction': model.reactions.keys(), 'minimum': np.nan, 'maximum': np.nan})
        else:
            logging.info("Model solved successfully")
            reactions, values = zip(*fva_solution.items())
            minimum, maximum = zip(*values)
            res_frame = pd.DataFrame({'reaction': reactions, 'minimum': minimum, 'maximum': maximum})
        res_frame["T"] = T
        rs.append(res_frame)
    return pd.concat(rs)

        

def fva_chemostat(model: CBModel,dilu: float,param_dict: dict,Ts: List[Float],sigma: float,
growth_id: str,glc_up_id: str,prot_pool_id: str) -> pd.DataFrame:
    '''
    # Do FVA simulation on a given dilution and a list of temperatures. 
    # model, reframed model
    # dilu, dilution rate
    # param_dict, a dictionary containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ts, a list of temperatures to simulate at. in K
    # sigma, saturation factor
    # growth_id, reaction of of growth
    # glc_up_id, reaction id of glucose uptake reaction
    # prot_pool_id, reaction id of prot_pool_exchange. --> prot_pool
    '''
    solutions: List[pd.DataFrame] = list() # corresponding to Ts. a list of solutions from model.optimize()
    # This wierd DataFrame with zero columns is intended to prevent nasty errors when all
    # FVA solutions fail
    skeleton_df = pd.DataFrame({'reaction': [],'T': [], 'minimum': [], 'maximum': []})
    solutions.append(skeleton_df)
    mappers = reframed_mappers
    m0 = model.copy()
    # Step 1: fix growth rate, set objective function as minimizing glucose uptatke rate
    rxn_growth = m0.reactions[growth_id]
    rxn_growth.lb = dilu
    m0.set_objective({glc_up_id: 1})

    for T in Ts:
        m1 = m0.copy()  
        # Step 2: map temperature constraints. 
        mappers.map_fNT(m1,T,param_dict)
        mappers.map_kcatT(m1,T,param_dict)
        mappers.set_NGAMT(m1,T)
        mappers.set_sigma(m1,sigma)
        
        try: 
            # Step 3: minimize the glucose uptake rate. Fix glucose uptake rate, minimize enzyme usage
            solution1 = reframed.FBA(m1, minimize=True)
            if solution1.status != reframed.solvers.solution.Status.OPTIMAL:
                raise OptimizationError(f"Solver status is {solution1.status.value}")
            m1.reactions[glc_up_id].ub = solution1.fobj*1.001
            m1.set_objective({prot_pool_id: -1})
            
            fva_solution = reframed.FVA(m1, obj_frac=0.999)
            logging.info('FVA problems solved successfully')
            reactions, values = zip(*fva_solution.items())
            minimum, maximum = zip(*values)
            res_frame = pd.DataFrame({'reaction': reactions, 'minimum': minimum, 'maximum': maximum})
        except Exception as err:
            logging.info(f'Failed to solve the problem, problem: {str(err)}')
            res_frame = pd.DataFrame({'reaction': model.reactions.keys(), 'minimum': np.nan, 'maximum': np.nan})
        res_frame["T"] = T
        solutions.append(res_frame)
    return pd.concat(solutions)
