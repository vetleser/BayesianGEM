from itertools import product
import logging
import os
import pickle
from typing import Dict, List
import pandas as pd
import reframed
from reframed import CBModel
from reframed.solvers.solver import Solver
import numpy as np
from sympy import solve
from scipy.optimize import fsolve
from GEMS import load_exp_batch_data

import etcpy.thermal_parameters as th

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")


# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def format_single_input(params, param_dict):
    """
    Format a single input parameter set for thermal parameters calculation.
    
    :param params: A dictionary containing the parameters.
    :param param_dict: A dictionary with keys as 'index_column_name_column_name' and values as new values.
    :return: A dictionary with formatted parameters.
    """
    new_params = params.copy()
    for key, val in param_dict.items():
        [ind, col] = key.split('_')
        new_params[ind][col] = val
    
    # Update T90
    new_params['T90'] = params['T90'] - params['Tm'] + new_params['Tm']
    
    return new_params

#From thermal_parameters:
def format_input(params, param_dict):
    new_params = params.copy()
    for key,val in param_dict.items():
        [ind,col] = key.split('_')
        new_params.loc[ind,col] = val
    
    # Update T90
    new_params['T90'] = params['T90']-params['Tm'] + new_params['Tm']
    df = calculate_thermal_params(new_params)

    return {id: {parameter: value for parameter, value in df.loc[id].iteritems()} for id in df.index}


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

    # step 3. Also copy Tm and T90, for convenience
    # thermalparams['Tm'] = params['Tm']
    # thermalparams['T90'] = params['T90']

    #logging.info(f"Thermal parameters calculated: \n {thermalparams}")
    
    return thermalparams

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
    R = 8.314
    TH = 373.5
    TS = 385
    T0 = 30+273.15

    # Use the equation from solvedHT.m and re-organized
    dGuTopt = dHTH +dCpu*(Topt-TH) -Topt*dSTS-Topt*dCpu*np.log(Topt/TS)
    dHt = dHTH+dCpu*(Topt-TH)-dCpt*(Topt-T0)-R*Topt-(dHTH+dCpu*(Topt-TH))/(1+np.exp(-dGuTopt/(R*Topt)))

    # Calculate kcat at reference Temperautre
    kcat0 = kcatTopt/np.exp(np.log(Topt/T0)-(dHt+dCpt*(Topt-T0))/R/Topt+dHt/R/T0+dCpt*np.log(Topt/T0)/R)

    # Calculate kcat at given temperature
    kcatT = kcat0*np.exp(np.log(T/T0)-(dHt+dCpt*(T-T0))/R/T+dHt/R/T0+dCpt*np.log(T/T0)/R)

    return kcatT


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

#Functions from thermal_parameters.py written above

def round_kcat(kcat, decimals=6):
    return round(kcat, decimals)

def save_results(kcat, tm, topt, dCpt, dict):
    """
    Save the results of kcat calculations to a file.
    
    :param kcat: The calculated kcat value.
    :param tm: The Tm value used in the calculation.
    :param topt: The Topt value used in the calculation.
    :param dCpt: The dCpt value used in the calculation.
    """

    kcat = float(kcat)  # Ensure kcat is a float for consistency
    kcat = round_kcat(kcat)  # Round kcat to 6 decimal places
    dCpt = float(dCpt)  # Ensure dCpt is a float for consistency
    if kcat in dict and kcat != 0:
        if abs(dict[kcat][0] - tm) > 0.1:
            logging.warning(f"Kcat {kcat} already exists in the dictionary. tm: {tm}, topt: {topt}, dCpt: {dCpt}")
            logging.warning(f"Existing entry: {dict[kcat]}")
        return
    dict[kcat] = [tm, topt, dCpt]


def func(tm, topt, dCpt, enzyme_id, T, reaction, T_dict):
    single_param = params[params.index.str.startswith(enzyme_id)].copy()
    new_single_param = single_param.copy()
    #logging.info(f"Shikimate parameters loaded: \n {shikimate_params}")

    new_single_param.loc[enzyme_id, 'Tm'] = tm
    new_single_param.loc[enzyme_id, 'Topt'] = topt
    new_single_param.loc[enzyme_id, 'dCpt'] = dCpt

    new_single_param['T90'] = single_param['T90'] - single_param['Tm'] + new_single_param['Tm']

    #logging.info(f"New shikimate parameters: \n {new_shikimate_params}")

    cols = ['dHTH', 'dSTS','dCpu','Topt','dCpt']

    df = calculate_thermal_params(new_single_param)
    #logging.info(f"Thermal parameters calculated: \n {df}")
    [dHTH, dSTS,dCpu,Topt,dCpt]= [df[parameter] for parameter in cols]
    #logging.info(f"Extracted parameters: dHTH={dHTH}, dSTS={dSTS}, dCpu={dCpu}, Topt={Topt}, dCpt={dCpt}")

    # logging.info(f"Thermal parameters calculated: \n {[dHTH, dSTS,dCpu,Topt,dCpt]}")
    kcatTopt = -1/model.reactions[reaction].stoichiometry[f'prot_{enzyme_id}']

    kcatT = calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt)
    #logging.info(f"Calculated kcatT for reaction {shikimate_kinase}. Tm {tm}, Topt {topt}, dCpt {(dCpt)}:\n {kcatT} for enzyme {enzyme_id} at T={T} K")
    #logging.info(f"Kcat: {float(kcatT)}")
    save_results(kcat=kcatT, tm=tm, topt=topt, dCpt=dCpt, dict=T_dict)

# 2025-06-10 18:12:54,811 P08566 Tm:    min = 311.74, max = 346.37
# 2025-06-10 18:12:54,811 P08566 Topt:  min = 271.20, max = 328.28
# 2025-06-10 18:12:54,811 P08566 dCpt:  min = -15273.10, max = 4263.70
tm_values = np.linspace(312, 346, 20)  # Example Tm values in K
topt_values = np.linspace(272, 328, 10)  # Example Tm values in K
dCpt_values = np.linspace(-15000, -1000, 10)  # Example dCpt values in J/mol/K


#logging.info(f"Model particle modified: \n {model_particle}")

path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
model = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
T = 310.0  # Temperature in K
dfae_batch,dfan_batch = load_exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))
sel_temp = dfae_batch.index + 273.15

kcat_dict = {}
rxn_dict : Dict[str, set]= {}
prot_dict : Dict[str, set] = {}
#Rxn r_0997No1, met prot_P08566: 1.8625408048767017e+30

rxns = model.reactions.values()


shikimate_kinase = 'r_0997No1'
enzyme_id = 'P08566'  # Example enzyme ID for testing
cols2 = ['Tm', 'Topt', 'dCpt']
kcat_dict = {}
counter = 0
# for tm in tm_values:
#     #logging.info(f"Processing Tm: {tm}")
#     counter += 1
#     logging.info(f"Counter: {counter}, Tm: {tm}")
#     for topt in topt_values:
#         for dCpt in dCpt_values:
#             if tm < topt: continue
#             func(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase)
#             shikimate_params = params[params.index.str.startswith(enzyme_id)].copy()
#             new_shikimate_params = shikimate_params.copy()
#             #logging.info(f"Shikimate parameters loaded: \n {shikimate_params}")

#             new_shikimate_params.loc[enzyme_id, 'Tm'] = tm
#             new_shikimate_params.loc[enzyme_id, 'Topt'] = topt
#             new_shikimate_params.loc[enzyme_id, 'dCpt'] = dCpt

#             new_shikimate_params['T90'] = shikimate_params['T90'] - shikimate_params['Tm'] + new_shikimate_params['Tm']

#             #logging.info(f"New shikimate parameters: \n {new_shikimate_params}")

#             cols = ['dHTH', 'dSTS','dCpu','Topt','dCpt']


#             df = calculate_thermal_params(new_shikimate_params)
#             #logging.info(f"Thermal parameters calculated: \n {df}")
#             [dHTH, dSTS,dCpu,Topt,dCpt]= [df[parameter] for parameter in cols]
#             #logging.info(f"Extracted parameters: dHTH={dHTH}, dSTS={dSTS}, dCpu={dCpu}, Topt={Topt}, dCpt={dCpt}")

#             # logging.info(f"Thermal parameters calculated: \n {[dHTH, dSTS,dCpu,Topt,dCpt]}")
#             kcatTopt = -1/model.reactions[shikimate_kinase].stoichiometry[f'prot_{enzyme_id}']

#             kcatT = calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt)
#             #logging.info(f"Calculated kcatT for reaction {shikimate_kinase}. Tm {tm}, Topt {topt}, dCpt {(dCpt)}:\n {kcatT} for enzyme {enzyme_id} at T={T} K")
#             #logging.info(f"Kcat: {float(kcatT)}")
#             save_results(kcat=kcatT, tm=tm, topt=topt, dCpt=dCpt, dict=kcat_dict)
                        

logging.info(f"Length of kcat_dict: {len(kcat_dict)}")
# Load the model and parameters
#param_dict = format_input(params, model_particle)
topt_values = [309.3333333333333]
dCpt_values = [-15000.0]

for T in sel_temp:
    T_dict = {}
    logging.info(f"Processing temperature: {T} K")
    for tm in tm_values:
        for topt, dCpt in product(topt_values, dCpt_values):
            if tm < topt: continue
            func(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase, T_dict=T_dict)





# rxn: reframed.CBReaction
# for rxn in model.reactions.values():
#     #if rxn.id.startswith('draw_prot'): continue
#     if rxn.id != shikimate_kinase: continue  # Only shikimate kinase for testing
#     for met in rxn.stoichiometry:
#         if not met.startswith('prot_'): continue
#         # ingore metabolite: prot_pool
#         if met == 'prot_pool': continue
#         uniprot_id = met.split('_')[1]
#         parameter_entries = param_dict[uniprot_id]
#         [dHTH, dSTS,dCpu,Topt,dCpt]= [parameter_entries[parameter] for parameter in cols]
#         [Tm, Topt, dCpt] = [parameter_entries[parameter] for parameter in cols2]
#         # Change kcat value.
#         # pmet_r_0001 + 1.8518518518518518e-07 prot_P00044 + 1.8518518518518518e-07 prot_P32891 -->
#         # 2.0 s_0710 + s_1399
#         #
#         # 1.8518518518518518e-07 is correponding to 1/kcat
#         # change the kcat to kcat(T)
#         # In some casese, this coefficient could be 2/kcat or some other values. This doesn't matter.
#         #
#         # a protein could be involved in several reactions
#         # assume that Topt in the original model is measured at Topt
#         kcatTopt = -1/rxn.stoichiometry[met]
#         kcatT = calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt)
#         #logging.info(f"Calculated kcatT for reaction {rxn.id}. Tm {Tm}, Topt {Topt}, dCpt {dCpt}:\n {kcatT} for {met} at T={T} K")
#         kcat_dict[f"Rxn {rxn.id}, met {met}"] = kcatT
#         if rxn.id not in rxn_dict:
#             rxn_dict[rxn.id] = set()
#         if met not in prot_dict:
#             prot_dict[met] = set()
#         rxn_dict[rxn.id].add(kcatT)
#         prot_dict[met].add(kcatT)


logging.info("DONE")