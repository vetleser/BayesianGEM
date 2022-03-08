from typing import List
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
import time
import cobra
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


def get_dH_dS_dCpu_from_TmT90(Tm,T90):
    '''
    # With knowing Tm and T90, get dHTH, dSTS and dCpu by solving
    # 
    # dHTH = slope*dSTS+intercept
    # deltaG(Tm) = 0
    # deltaG(T90) = -RTln9
    # 
    # to get dSTS, dCpu.
    # 
    # Tm, T90 are in K
    # dHTH, is in J/mol
    # dSTS is in J/mol/K
    # 
    '''
    TH = 373.5
    TS = 385
    R = 8.314
    slope = 299.58
    intercept = 20008
    
    a = np.array([[1,-slope,0],
                 [1,-Tm,Tm-TH-Tm*np.log(Tm/TS)],
                 [1,-T90,T90-TH-T90*np.log(T90/TS)]])
    b = np.array([intercept,0,-R*T90*np.log(9)])
    
    [dHTH,dSTS,dCpu] = list(np.linalg.solve(a,b))

    return dHTH,dSTS,dCpu

def get_dH_dS_dCpu_from_TmLength(Tm,N):
    '''
    # In case of negative obtained from get_dH_dS_dCpu_from_TmT90(Tm,T90), or this is no T90 data available, 
    # using the same euqations from  Sawle and Ghosh, Biophysical Journal, 2011 for delatH* and deltaS*.
    # Then caculate dCpu by solving deltaG(Tm) =0
    # Tm is in K
    # 
    # Tm, T90 are in K
    # dHTH, is in J/mol
    # dSTS is in J/mol/K
    # 
    '''
    dHTH = (4*N+143)*1000
    dSTS = 13.27*N+448
    
    TH = 373.5
    TS = 385
    
    def func(dCp):
        dGTm = dHTH + dCp*(Tm-TH)-Tm*dSTS -Tm*dCp*np.log(Tm/TS)
        return dGTm
    dCpu = fsolve(func,10000)[0]
    return dHTH,dSTS,dCpu

def change_rxn_coeff(rxn,met,new_coeff):
    '''
    # This is based on the rxn.add_metabolites function. If there the metabolite is already in the reaction,
    # new and old coefficients will be added. For example, if the old coeff of metA is 1, use
    # rxn.add_metabolites({metA:2}), After adding, the coeff of metA is 1+2 = 3
    #
    '''

    diff_coeff = new_coeff-rxn.metabolites[met]
    rxn.add_metabolites({met:diff_coeff})



def get_dGu(T: float,dHTH: float,dSTS: float,dCpu: float):
    '''
    # calculate the deltaG of unfolding process at temperature T
    # dHTH, dSTS are enthalpy and entropy at TH and TS, repspectively
    # dCpu is the heat capacity change unpoin unfolding
    '''
    TH = 373.5;
    TS = 385;

    dGu: float = dHTH +dCpu*(T-TH) -T*dSTS-T*dCpu*np.log(T/TS);
    return dGu


def get_fNT(T: float,dHTH: float,dSTS: float,dCpu: float):
    '''
    # Calculate the fraction of enzyme in native state
    # dHTH, dSTS are enthalpy and entropy at TH and TS, repspectively
    # dCpu is the heat capacity change unpoin unfolding
    # T, temeperature, a single value or a numpy array, in K
    '''
    R = 8.314;
    dGu = get_dGu(T,dHTH,dSTS,dCpu);
    f: float = 1/(1+np.exp(-dGu/R/T));
    return f


def map_fNT(model: Model,T: float,df: pd.DataFrame,Tadj: float=0, reference_model: Model=None):
    '''
    # apply the fraction of enzymes in native state to each protein.
    # model, cobra model
    # T, temperature, in K, float
    # Tadj, This is to adjust the orginal denaturation curve by moving to left by
    # Tadj degrees.
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu
    # reference_model In the case of warm-starting FBA, the parameters to modify are fetched from the
    # reference model cobra model.
    #
    #
    # Gang Li, 2019-05-03


    # in the enzyme constrained model, protein usage was describe as reaction
    # draw_prot_P0CW41 68.69778 prot_pool --> prot_P0CW41. 68.69778 is the molecular weight of protein P0CW41
    # the total enzyme pool was constrained by the upper bound of reaction prot_pool_exchange:
    #          --> prot_pool
    # To map fNT, change the coefficient of prot_pool in reactions like MW*prot_pool --> prot_P0CW41,
    #                                MW/fNT prot_pool --> prot_P0CW41

    '''

    met: cobra.Metabolite = model.metabolites.prot_pool
    if reference_model is not None:
        reference_met: cobra.Metabolite = reference_model.metabolites.prot_pool

    rxn: cobra.Reaction
    for rxn in met.reactions:
        # this is to ignore reaction 'prot_pool_exchange': --> prot_pool
        if len(rxn.metabolites)<2: continue

        uniprot_id: str = rxn.id.split('_')[-1]
        cols = ['dHTH', 'dSTS','dCpu','Topt']
        [dHTH, dSTS,dCpu,topt] = df.loc[uniprot_id,cols]
        fNT = get_fNT(T+Tadj,dHTH,dSTS,dCpu)
        if fNT < 1e-32: fNT = 1e-32
        if reference_model is not None:
            rxn_id = rxn.id
            reference_rxn: float = reference_model.reactions.get_by_id(rxn_id)
            new_coeff: float = reference_rxn.metabolites[reference_met]/fNT
        else: 
            new_coeff: float = rxn.metabolites[met]/fNT
        change_rxn_coeff(rxn,met,new_coeff)
        


def calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt):
    '''
    # Using Trainsition state theory to calculate kcat at temperature T.
    # dHTH, dSTS: entropy and enthalpy at comergence temperatures. Protein
    # unfolding process.
    # dCpu, heat capacity change unpon unfolding.
    # kcatTopt: kcat values at optimal temperature
    # Topt, optimal temperature of the enzyme, in K
    # T, temperature, in K
    #
    '''
    # Constants
    R = 8.314;
    TH = 373.5;
    TS = 385;
    T0 = 30+273.15;

    # Use the equation from solvedHT.m and re-organized
    dGuTopt = dHTH +dCpu*(Topt-TH) -Topt*dSTS-Topt*dCpu*np.log(Topt/TS);
    dHt = dHTH+dCpu*(Topt-TH)-dCpt*(Topt-T0)-R*Topt-(dHTH+dCpu*(Topt-TH))/(1+np.exp(-dGuTopt/(R*Topt)));

    # Calculate kcat at reference Temperautre
    kcat0 = kcatTopt/np.exp(np.log(Topt/T0)-(dHt+dCpt*(Topt-T0))/R/Topt+dHt/R/T0+dCpt*np.log(Topt/T0)/R);

    # Calculate kcat at given temperature
    kcatT = kcat0*np.exp(np.log(T/T0)-(dHt+dCpt*(T-T0))/R/T+dHt/R/T0+dCpt*np.log(T/T0)/R);

    return kcatT


def map_kcatT(model: Model,T: float,df: pd.DataFrame, reference_model: Model=None):
    '''
    # Apply temperature effect on enzyme kcat.
    # based on trainsition state theory
    # model, cobra model
    # T, temperature, in K
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # reference_model In the case of warm-starting FBA, the parameters to modify are fetched from the
    # reference model cobra model.
    # Ensure that Topt is in K. Other parameters are in standard units.
    #
    # Gang Li, 2019-05-03
    #
    '''
    met: cobra.Metabolite
    for met in model.metabolites:
        

        # look for those metabolites: prot_uniprotid
        if not met.id.startswith('prot_'): continue

        # ingore metabolite: prot_pool
        if met.id == 'prot_pool': continue
        uniprot_id = met.id.split('_')[1]

        if reference_model is not None:
            met_id = met.id
            reference_met = reference_model.metabolites.get_by_id(met_id)

        # Change kcat value.
        # pmet_r_0001 + 1.8518518518518518e-07 prot_P00044 + 1.8518518518518518e-07 prot_P32891 -->
        # 2.0 s_0710 + s_1399
        #
        # 1.8518518518518518e-07 is correponding to 1/kcat
        # change the kcat to kcat(T)
        # In some casese, this coefficient could be 2/kcat or some other values. This doesn't matter.
        #
        # a protein could be involved in several reactions
        cols = ['dHTH', 'dSTS','dCpu','Topt','dCpt']
        [dHTH, dSTS,dCpu,Topt,dCpt]=df.loc[uniprot_id,cols]

        rxn: cobra.Reaction
        for rxn in met.reactions:
            if rxn.id.startswith('draw_prot'): continue

            if reference_model is not None:
                # assume that Topt in the original model is measured at Topt
                rxn_id = rxn.id
                reference_rxn = reference_model.reactions.get_by_id(rxn_id)
                kcatTopt = -1/reference_rxn.metabolites[reference_met]
            else: 
                # assume that Topt in the original model is measured at Topt
                kcatTopt = -1/rxn.metabolites[met]
            kcatT = calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt)
            if kcatT < 1e-32: kcatT = 1e-32
            new_coeff = -1/kcatT

            change_rxn_coeff(rxn,met,new_coeff)

def getNGAMT(T):
    # T is in K, a single value
    def NGAM_function(T):
        return 0.740 + 5.893/(1+np.exp(31.920-(T-273.15))) + 6.12e-6*(T-273.15-16.72)**4

    lb = 5+273.15
    ub = 40+273.15
    if T < lb: NGAM_T = NGAM_function(lb)
    elif T > ub: NGAM_T =NGAM_function(ub)
    else: NGAM_T = NGAM_function(T)

    return NGAM_T


def set_NGAMT(model,T):
    # T is in K
    NGAM_T = getNGAMT(T)
    rxn = model.reactions.NGAM
    #ori_lb,ori_ub = rxn.lower_bound,rxn.upper_bound
    rxn.lower_bound = NGAM_T
    rxn.upper_bound = NGAM_T


def set_sigma(model,sigma):
    rxn = model.reactions.prot_pool_exchange
    #ori_ub_sigma = rxn.upper_bound
    rxn.upper_bound = 0.17866*sigma


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
                    r = working_model.optimize(raise_error=True).objective_value
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
    solutions = list() # corresponding to Ts. a list of solutions from model.optimize()
    df = calculate_thermal_params(params)
    with model as m0:
        # Step 1: fix growth rate, set objective function as minimizing glucose uptake rate
        rxn_growth = m0.reactions.get_by_id(growth_id)
        rxn_growth.lower_bound = dilu

        m0.objective = glc_up_id
        m0.objective.direction = 'min'
        if warm_start:
            working_model.objective = glc_up_id
            working_model.objective.direction = 'min'
        
        for T in Ts:
            with m0 as m1:  
                # Step 2: map temperature constraints. 
                if warm_start:
                    map_fNT(working_model,T,df, reference_model=m1)
                    map_kcatT(working_model,T,df, reference_model=m1)
                    set_NGAMT(working_model,T)
                    set_sigma(working_model,sigma)
                else:
                    map_fNT(m1,T,df)
                    map_kcatT(m1,T,df)
                    set_NGAMT(m1,T)
                    set_sigma(m1,sigma)
                
                try: 
                    # Step 3: minimize the glucose uptake rate. Fix glucose uptake rate, minimize enzyme usage
                    model_to_optimize = working_model if warm_start else m1
                    solution1 = model_to_optimize.optimize(raise_error=True)
                    model_to_optimize.reactions.get_by_id(glc_up_id).upper_bound = solution1.objective_value*1.001
                    model_to_optimize.objective = prot_pool_id
                    model_to_optimize.objective.direction = 'min'
                    
                    solution2 = model_to_optimize.optimize(raise_error=True)
                    solutions.append(solution2)
                    logging.info('Model solved successfully')
                except OptimizationError as err:
                    logging.info(f'Failed to solve the problem, problem: {str(err)}')
                    #solutions.append(None)
                    break # because model has been impaired. Further simulation won't give right output.
    return solutions
