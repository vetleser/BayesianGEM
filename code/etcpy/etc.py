from typing import List
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
import time
import cobra
import reframed
from cobra import Model
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis
import logging

from sympy import Float

T0 = 273.15

def solve_unboundedness(model: Model, extreme_bound: float = 1000):
    """This function ensures no bound has an infinite magnitude, but are set to a threshold
    The input model is changed in-place

    Args:
        model (Model): A Cobra model to the modified
        extreme_bound (float): The extreme absolute value for bounds. Any bound exceeding this will be set
        to this value respecting the its sign
    """
    for reaction in model.reactions:
        if reaction.lower_bound < extreme_bound:
            reaction.lower_bound = - extreme_bound
        if reaction.upper_bound > extreme_bound:
            reaction.upper_bound = extreme_bound        






def simulate_fva(model: Model,Ts: List[float],sigma: float,df: pd.DataFrame,Tadj=0,processes: int=1) -> pd.DataFrame:
    '''
    Simulate FVA on growth scenarios
    # model, cobra model
    # Ts, a list of temperatures in K
    # sigma, enzyme saturation factor
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ensure that Topt is in K. Other parameters are in standard units.
    # Tadj, as descrbed in map_fNT
    # 
    '''
    rs: List[pd.DataFrame] = list()
    # This wierd DataFrame with zero columns is intended to prevent nasty errors when all
    # FVA solutions fail
    skeleton_df = pd.DataFrame({'reaction': [],'T': [], 'minimum': [], 'maximum': []})
    rs.append(skeleton_df)
    for T in Ts:
        with model:
            # map temperature constraints
            map_fNT(model,T,df)
            map_kcatT(model,T,df)
            set_NGAMT(model,T)
            set_sigma(model,sigma)
            r = flux_variability_analysis(model=model, fraction_of_optimum=0.99, processes=processes)
            logging.info("Model solved successfully")
            r["T"] = T
            r.reset_index(inplace=True)
            r.rename(columns={'index': 'reaction'},inplace=True)
            rs.append(r)
            """
            try:
                    r = flux_variability_analysis(model=model, fraction_of_optimum=0.99)
                    logging.info("Model solved successfully")
                    r["T"] = T
                    r.reset_index(inplace=True)
                    r.rename(columns={'index': 'reaction'},inplace=True)
                    rs.append(r)
            except OptimizationError as err:
                logging.info(f'Failed to solve the problem, problem: {str(err)}')
            """
    return pd.concat(rs)



def simulate_growth(model: Model,Ts,sigma,df,Tadj=0, working_model: Model=None):
    '''
    # model, cobra model
    # Ts, a list of temperatures in K
    # sigma, enzyme saturation factor
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # working_model, If provided, warm-start FBA will be used for accelerating computations. The working model will be modified during this
    # process
    # Ensure that Topt is in K. Other parameters are in standard units.
    # Tadj, as descrbed in map_fNT
    #
    '''
    warm_start = working_model is not None
    rs = list()
    for T in Ts:
        with model:
            # map temperature constraints
            if warm_start:
                map_fNT(working_model,T,df, reference_model=model)
                map_kcatT(working_model,T,df, reference_model=model)
                set_NGAMT(working_model,T)
                set_sigma(working_model,sigma)
            else:
                map_fNT(model,T,df)
                map_kcatT(model,T,df)
                set_NGAMT(model,T)
                set_sigma(model,sigma)

            try:
                if warm_start:
                    r = reframed.FBA(working_model)
                else:
                    r = model.optimize(raise_error=True).objective_value
                logging.info("Model solved successfully")
            except OptimizationError as err:
                logging.info(f'Failed to solve the problem, problem: {str(err)}')
                r = 0
            print(T-273.15,r)
            rs.append(r)
    return rs

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

def calculate_thermal_params(params):
    '''
    # params, a dataframe with at least following columns: Tm,T90,Length,dCpt,Topt,dCpt. All are in standard units.
    # 
    # The script will return a dataframe with following columns: dHTH,dSTS,dCpu,Topt,dCpt
    # 
    '''
    thermalparams = pd.DataFrame()
    
    # step 1: calculate dHTH,dSTS,dCpu from tm, t90/length
    for ind in params.index:
        tm,t90 = params.loc[ind,'Tm'],params.loc[ind,'T90']
        if t90 is not None and t90>tm and ~np.isnan(t90):
            dHTH,dSTS,dCpu = get_dH_dS_dCpu_from_TmT90(tm,t90)
            if dCpu <0: dHTH,dSTS,dCpu = get_dH_dS_dCpu_from_TmLength(tm,params.loc[ind,'Length'])
        else: dHTH,dSTS,dCpu = get_dH_dS_dCpu_from_TmLength(tm,params.loc[ind,'Length'])
         
        thermalparams.loc[ind,'dHTH'] = dHTH
        thermalparams.loc[ind,'dSTS'] = dSTS
        thermalparams.loc[ind,'dCpu'] = dCpu
        
        
    # step 2. copy columns Topt and dCpt
    thermalparams['Topt'] = params['Topt']
    thermalparams['dCpt'] = params['dCpt']
    
    return thermalparams
        

def fva_chemostat(model: Model,dilu: float,params: pd.DataFrame,Ts: List[Float],sigma: float,
growth_id: str,glc_up_id: str,prot_pool_id: str, processes: int=1) -> pd.DataFrame:
    '''
    # Do FVA simulation on a given dilution and a list of temperatures. 
    # model, cobra model
    # dilu, dilution rate
    # params: a dataframe containing Tm, T90, Length, dCpt, Topt. All temperatures are in K.
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
    df = calculate_thermal_params(params)
    with model as m0:
        # Step 1: fix growth rate, set objective function as minimizing glucose uptatke rate
        rxn_growth = m0.reactions.get_by_id(growth_id)
        rxn_growth.lower_bound = dilu

        m0.objective = glc_up_id
        m0.objective.direction = 'min'

        for T in Ts:
            with m0 as m1:  
                # Step 2: map temperature constraints. 
                map_fNT(m1,T,df)
                map_kcatT(m1,T,df)
                set_NGAMT(m1,T)
                set_sigma(m1,sigma)
                
                try: 
                    # Step 3: minimize the glucose uptake rate. Fix glucose uptake rate, minimize enzyme usage
                    solution1 = m1.optimize(raise_error=True)
                    m1.reactions.get_by_id(glc_up_id).upper_bound = solution1.objective_value*1.001
                    m1.objective = prot_pool_id
                    m1.objective.direction = 'min'
                    
                    fva_solution = flux_variability_analysis(m1,fraction_of_optimum=0.99, processes=processes)
                    logging.info('FVA problems solved successfully')
                    fva_solution["T"] = T
                    fva_solution.reset_index(inplace=True)
                    fva_solution.rename(columns={'index': 'reaction'},inplace=True)
                    solutions.append(fva_solution)
                except OptimizationError as err:
                    logging.info(f'Failed to solve the problem, problem: {str(err)}')
                    #solutions.append(None)
                    break # because model has been impaired. Further simulation won't give right output.
                
    return pd.concat(solutions)


def simulate_chomostat(model: Model,dilu,params,Ts,sigma,growth_id,glc_up_id,prot_pool_id, working_model: Model=None):
    '''
    # Do simulation on a given dilution and a list of temperatures. 
    # model, cobra model
    # dilu, dilution rate
    # params: a dataframe containing Tm, T90, Length, dCpt, Topt. All temperatures are in K.
    # Ts, a list of temperatures to simulate at. in K
    # sigma, saturation factor
    # growth_id, reaction of of growth
    # glc_up_id, reaction id of glucose uptake reaction
    # prot_pool_id, reaction id of prot_pool_exchange. --> prot_pool
    # working_model, If provided, warm-start FBA will be used for accelerating computations. The working model will be modified during this
    # process

    '''
    warm_start = working_model is not None
    reference_model = model if warm_start else None
    solutions = list() # corresponding to Ts. a list of solutions from model.optimize()
    df = calculate_thermal_params(params)
    model_to_optimize = working_model if warm_start else model
    with model_to_optimize as m0:
        # Step 1: fix growth rate
        rxn_growth = m0.reactions.get_by_id(growth_id)
        rxn_growth.lower_bound = dilu
        for T in Ts:
            with m0 as m1:
                opt_model = m0 if warm_start else m1
                # Step 2: map temperature constraints. 
                map_fNT(opt_model,T,df, reference_model=reference_model)
                map_kcatT(opt_model,T,df, reference_model=reference_model)
                set_NGAMT(opt_model,T)
                set_sigma(opt_model,sigma)
                try:
                    # Step 3: set objective function as minimizing glucose uptake rate and protein useage 
                    cobra.util.add_lexicographic_constraints(opt_model, [glc_up_id, prot_pool_id], ['min', 'min'])
                    solutions.append(opt_model.optimize())
                    logging.info('Model solved successfully')
                except OptimizationError as err:
                    logging.info(f'Failed to solve the problem, problem: {str(err)}')
                    #solutions.append(None)
                    break # because model has been impaired. Further simulation won't give right output.
    return solutions
