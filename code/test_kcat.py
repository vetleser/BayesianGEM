from collections import defaultdict
from itertools import product
import logging
import os
import pickle
import time
from typing import Dict, List, Tuple
import pandas as pd
import reframed
from reframed import CBModel
from reframed.solvers.solver import Solver
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics import euclidean_distances
from sklearn.preprocessing import StandardScaler
from sympy import solve
from scipy.optimize import fsolve
from GEMS import load_exp_batch_data
import matplotlib.pyplot as plt
import matplotlib.cm as cm



import etcpy.thermal_parameters as th

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
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

def merge_matches(matches):
    clusters = []
    seen = set()

    for match in matches:
        group = set(match)
        if not group & seen:
            clusters.append(group)
            seen.update(group)
        else:
            # Merge overlapping groups
            merged = False
            for i, c in enumerate(clusters):
                if c & group:
                    clusters[i] = c | group
                    seen.update(group)
                    merged = True
                    break
            if not merged:
                clusters.append(group)
                seen.update(group)
    return clusters

def analyze_cluster_diversity(cluster, parameter_list):
    tm_vals = [parameter_list[i][0] for i in cluster]
    topt_vals = [parameter_list[i][1] for i in cluster]
    dCpt_vals = [parameter_list[i][2] for i in cluster]

    return {
        "size": len(cluster),
        "Tm_range": max(tm_vals) - min(tm_vals),
        "Topt_range": max(topt_vals) - min(topt_vals),
        "dCpt_range": max(dCpt_vals) - min(dCpt_vals),
        "Tm_std": np.std(tm_vals),
        "Topt_std": np.std(topt_vals),
        "dCpt_std": np.std(dCpt_vals),
    }

def make_plot(T_dict, T, enzyme_id, reaction):
    plt.figure(figsize=(10, 6))

    for kcat, temps in T_dict.items():
        plt.plot(temps, [kcat] * len(temps), 'o-', label=f'kcat={kcat:.2f}' if len(temps) > 1 else "", alpha=0.7)

    plt.title(f"Calculated $k_{{cat}}$ for different $T_m$ at T={T} K")
    plt.xlabel("$T_m$ (K)")
    plt.ylabel("$k_{cat}$ (s$^{-1}$)")
    plt.grid(True, linestyle='--', alpha=0.5)

    # Only show legend for repeated kcat values (if any)
    handles, labels = plt.gca().get_legend_handles_labels()
    if labels:
        plt.legend()

    plt.tight_layout()
    plt.show()
    plt.savefig(f'../figures/analysis/kcat_vs_tm_{enzyme_id}_{reaction}_T{T}.png')

def make_diversity_plot(cluster_idx, sel_temp, max_legend_entries=10):
    """
    Plot the kcat (log10) curves for all parameter sets in a specific cluster.

    :param cluster_idx: Index of the cluster to plot.
    :param sel_temp: List of selected temperatures.
    """
    cluster = clusters[cluster_idx]
    
    cluster = list(cluster)

    plt.figure(figsize=(10, 6))
    
    for idx_in_cluster, i in enumerate(cluster):
        params = parameter_list[i]
        label = f"Tm={params[0]:.1f}, Topt={params[1]:.1f}, dCpt={params[2]:.2f}" if idx_in_cluster < max_legend_entries else None
        true_kcat = 10 ** kcat_matrix[i] -1

        plt.plot(sel_temp, true_kcat, alpha=0.5, label=label)
        if idx_in_cluster == max_legend_entries:
            break
         

    plt.xlabel("Temperature (°C)")
    plt.ylabel("kcat")
    plt.title(f"Kcat Profiles for Cluster {cluster_idx} ({len(cluster)} members)")
    plt.grid(True)
    
    # Only show legend if small enough
    if len(cluster) <= max_legend_entries:
        plt.legend(loc='best', fontsize='small', frameon=False)
    elif max_legend_entries > 0:
        plt.legend(loc='best', fontsize='small', frameon=False, title=f"First {max_legend_entries} parameter sets")
    
    plt.tight_layout()
    filename = f'../figures/analysis/cluster_{cluster_idx}_diversity_plot.png'
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved cluster diversity plot with legend to: {filename}")

def make_diversity_plot_from_simdata(result_dict, sel_temp, max_legend_entries=10):
    """
    Plot true kcat curves for multiple parameter sets from a result dictionary.
    
    :param result_dict: Dict with keys as (Tm, Topt, dCpt) and values as log10(kcat + 1) lists.
    :param sel_temp: List of temperatures corresponding to each kcat point.
    :param max_legend_entries: Maximum number of curves to include in the legend.
    """
    plt.figure(figsize=(10, 6))

    for idx, (params, kcat_list) in enumerate(result_dict.items()):
        # Convert log10(kcat + 1) back to true kcat
        
        label = f"Tm={params[0]:.1f}, Topt={params[1]:.1f}, dCpt={params[2]:.2f}" if idx < max_legend_entries else None
        plt.plot(sel_temp, kcat_list, alpha=0.6, label=label)
        if idx == max_legend_entries:
            break

    plt.xlabel("Temperature (K)")
    plt.ylabel("$k_{cat}$")
    plt.yscale('log')  # Log scale for kcat
    plt.title(rf"$k_{{cat}}$ Profiles for {(max_legend_entries)} Parameter Sets")
    plt.grid(True)

    # Add legend if there are few enough curves
    if len(result_dict) <= max_legend_entries:
        plt.legend(loc='best', fontsize='small', frameon=False)
    elif max_legend_entries > 0:
        plt.legend(loc='best', fontsize='small', frameon=False, title=f"First {max_legend_entries} parameter sets")

    plt.tight_layout()
    filename = '../figures/analysis/simdata_diversity_plot.png'
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved diversity plot for simdata to: {filename}")

        

def round_kcat(kcat, decimals=6):
    return round(kcat, decimals)

def plot_combined_kcats(all_T_dicts, sel_temp, enzyme_id, reaction):
    plt.figure(figsize=(10, 6))

    colors = cm.viridis(np.linspace(0, 1, len(sel_temp))) # type: ignore

    # For tm
    # for T, T_dict, color in zip(sel_temp, all_T_dicts, colors):
    #     for kcat, tm_list in T_dict.items():
    #         if len(tm_list) > 1:
    #             plt.plot(tm_list, [kcat] * len(tm_list), 'o-', color=color, alpha=0.7)
    #         else:
    #             plt.plot(tm_list, [kcat], 'o--', color=color, alpha=0.7)
    # #  Legend with one entry per T
    # for T, color in zip(sel_temp, colors):
    #     plt.plot([], [], 'o-', color=color, label=f"T = {T} K")


    for T, T_dict, color in zip(sel_temp, all_T_dicts, colors):
    # Extract all (Tm, kcat) pairs for this T
        points = []
        for kcat, param_list in T_dict.items():
            for param in param_list:
                points.append((param, kcat))
        # Sort points by Tm (x-axis)
        points.sort(key=lambda x: x[0])
        params, kcats = zip(*points)  # unzip into two lists
        
        plt.plot(params, kcats, 'o--', color=color, alpha=0.7, label=f"T = {T} K")


    
    plt.title(rf"Calculated $k_{{cat}}$ for different $\Delta C_p^\ddag$")
    plt.xlabel(r"$\Delta C_p^\ddag$ (J/mol/K)")
    plt.ylabel("$k_{cat}$ (s$^{{-1}}$)")
    #plt.yscale('log')  # Log scale for kcat
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(title="Fixed T", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

    # Optional: Save figure
    plt.savefig(f'../figures/analysis/kcat_vs_dcpt_combined_{enzyme_id}_{reaction}.png')



def save_results(kcat, tm, topt, dCpt, dict):
    """
    Save the results of kcat calculations to a file.
    
    :param kcat: The calculated kcat value.
    :param tm: The Tm value used in the calculation.
    :param topt: The Topt value used in the calculation.
    :param dCpt: The dCpt value used in the calculation.
    """

    kcat = float(kcat)  # Ensure kcat is a float for consistency
    kcat = round_kcat(kcat, 0)  # Round kcat to 6 decimal places
    dCpt = float(dCpt)  # Ensure dCpt is a float for consistency

    #param_to_evaluate = dCpt
    if kcat in dict:
        dict[kcat].append([tm, topt, dCpt])
        counter_dict[kcat] = counter_dict.get(kcat, 0) + 1
    else:
        dict[kcat] = [[tm, topt, dCpt]]
        counter_dict[kcat] = 1
    
def save_results2(kcat, tm, topt, dCpt, results_dict):
    """
    Save the results of kcat calculations to a file.
    
    :param kcat: The calculated kcat value.
    :param tm: The Tm value used in the calculation.
    :param topt: The Topt value used in the calculation.
    :param dCpt: The dCpt value used in the calculation.
    """

    kcat = float(kcat)  # Ensure kcat is a float for consistency
    kcat = round_kcat(kcat, 6)  # Round kcat to 6 decimal places
    dCpt = float(dCpt)  # Ensure dCpt is a float for consistency

    key = (tm, topt, dCpt)  # Use a tuple as the key for better readability

    if key in results_dict:
        logging.warning(f"Key {key} already exists in results_dict.")
    else:
        results_dict[key] = kcat


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

def func2(tm, topt, dCpt, enzyme_id, T, reaction, T_dict): #Identical to func, but with different save_results function
    
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
    save_results2(kcat=kcatT, tm=tm, topt=topt, dCpt=dCpt, results_dict=T_dict)

# 2025-06-10 18:12:54,811 P08566 Tm:    min = 311.74, max = 346.37
# 2025-06-10 18:12:54,811 P08566 Topt:  min = 271.20, max = 328.28
# 2025-06-10 18:12:54,811 P08566 dCpt:  min = -15273.10, max = 4263.70


#logging.info(f"Model particle modified: \n {model_particle}")

path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
model = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
T = 310.0  # Temperature in K
dfae_batch,dfan_batch = load_exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))
sel_temp = dfae_batch.index + 273.15


#Rxn r_0997No1, met prot_P08566: 1.8625408048767017e+30



shikimate_kinase = 'r_0997No1'
enzyme_id = 'P08566'  # Example enzyme ID for testing
cols2 = ['Tm', 'Topt', 'dCpt']
kcat_dict = {}
counter = 0

tm_values = np.linspace(312, 346, 50)  # Example Tm values in K
topt_values = np.linspace(272, 328, 50)  # Example Tm values in K
dCpt_values = np.linspace(-15000, 1000, 50) 
counter_dict = {}
result_dict = {}


# for tm, topt, dCpt in product(tm_values, topt_values, dCpt_values):
#     if tm < topt: continue  # Skip if Tm is less than Topt
#     counter += 1
#     if counter % 10000 == 0 or (counter %100 ==0 and counter < 1000):  # Log every 100 iterations
#         logging.info(f"Processed {counter} combinations so far.")
#         logging.info(f"Length of T_dict: {len(result_dict)}")
#         # Check for repeated kcat values
#         logging.info("Checking for repeated kcat values:")
#         # Log kcat values that have more than one entry
#         logging.info(f"Counter dictionary length: {len(counter_dict)}")
#         logging.info(f"Maximum counter value: {max(counter_dict.values()) if counter_dict else 0}")
#     #logging.info(f"Processing combination {counter}: Tm={tm}, Topt={topt}, dCpt={dCpt}")
        
#     # Call the function with the current parameters
#     func(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase, T_dict=result_dict)

# for kcat, values in counter_dict.items():
#     if values > 1:
#         logging.info(f"Kcat {kcat} has {values} entries in the dictionary.")
#         logging.info(f"Entries: {result_dict[kcat]}")

# dump_pickle(result_dict, f'../results/analysis/{enzyme_id}_kcat_dict.pkl')
# dump_pickle(counter_dict, f'../results/analysis/{enzyme_id}_counter_dict.pkl')


#------------------------------------------------------------------------------------
kcat_dict1 = load_pickle(f'../results/analysis/{enzyme_id}_kcat_dict.pkl')
counter_dict1 = load_pickle(f'../results/analysis/{enzyme_id}_counter_dict.pkl')

tm_tol = (346 - 312)/5 # = 6.8
topt_tol = (328 - 272)/5 # = 11.2

kcat_min = min(kcat_dict1.keys())
kcat_max = max(kcat_dict1.keys())
logging.info(f"Minimum kcat: {kcat_min}, Maximum kcat: {kcat_max}")


# for kcat, values in counter_dict1.items():
#     if kcat < 1000: continue  # Skip kcat values less than 1e-6
#     if values > 1:
#         entries = kcat_dict1[kcat]
#         tm_min = min(entry[0] for entry in entries)
#         tm_max = max(entry[0] for entry in entries)
#         topt_min = min(entry[1] for entry in entries)
#         topt_max = max(entry[1] for entry in entries)
#         dCpt_min = min(entry[2] for entry in entries)
#         dCpt_max = max(entry[2] for entry in entries)
#         tm_diff = tm_max - tm_min
#         topt_diff = topt_max - topt_min
#         dCpt_diff = dCpt_max - dCpt_min
#         if tm_tol < tm_diff and topt_tol < topt_diff :
#             logging.info(f"Tm and Topt difference for kcat {kcat} with {values} entries: {entries}")
#             logging.info(f"Tm_diff: {tm_max - tm_min}")
#             logging.info(f"Topt_diff: {topt_max - topt_min}")
#             logging.info(f"dCpt_diff: {dCpt_max - dCpt_min}")
#         # if topt_diff > 10:
#         #     logging.info(f"Topt difference for kcat {kcat} with {values} entries: {entries}")
#         #     logging.info(f"Topt_diff: {topt_max - topt_min}")
#         # if dCpt_diff > 1000:
#         #     logging.info(f"dCpt difference for kcat {kcat} with {values} entries: {entries}")
#         #     logging.info(f"dCpt_diff: {dCpt_max - dCpt_min}")




#--------------------------------------------------------------------------------------------------------------------
# all_T_dicts = []
# for T in sel_temp:
#     counter = 0
#     logging.info(f"Processing temperature: {T} K")
#     result_dict: Dict[Tuple, List[float]] = {}
#     for tm, topt, dCpt in product(tm_values, topt_values, dCpt_values):
#         if tm < topt: continue  # Skip if Tm is less than Topt
#         if counter % 10000 == 0 or (counter %100 ==0 and counter < 1000):  # Log every 100 iterations
#             logging.info(f"Processed {counter} combinations so far.")
#             logging.info(f"Maximum kcat: {max(result_dict.keys()) if result_dict else 'N/A'} = ")
#             # Check for repeated kcat values
#             #logging.info(f"Processing combination {counter}: Tm={tm}, Topt={topt}, dCpt={dCpt}")

#         # Call the function with the current parameters
#         func2(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase, T_dict=result_dict)
#         counter += 1
#     dump_pickle(result_dict, f'../results/analysis/{enzyme_id}_result_dict_T{T}.pkl')
#     all_T_dicts.append(result_dict)

# dump_pickle(all_T_dicts, f'../results/analysis/{enzyme_id}_all_T_dicts.pkl')

#-------------------------------------------------------------------------------------------------------------------

# all_T_dicts = load_pickle(f'../results/analysis/{enzyme_id}_all_T_dicts.pkl')
# combined_T = {key: [] for key in all_T_dicts[0].keys()}
# for d in all_T_dicts:
#     for params, kcat in d.items():
#         log_kcat = np.log10(kcat+1)  
#         combined_T[params].append(log_kcat)

# for values in combined_T.values():
#     if len(values) != len(sel_temp):
#         logging.warning(f"Length mismatch for key {params}: expected {len(sel_temp)}, got {len(values)}")

# logging.info(f"Combined T dictionary length: {len(combined_T)}")
# logging.info(f"Combined T dictionary keys: {list(combined_T.keys())[:10]}")  # Log first 10 keys for brevity
# logging.info(f"Combined T dictionary values: {list(combined_T.values())[:10]}")  # Log first 10 values for brevity

# logging.info("Convert combined_T to a scaled matrix for clustering")
# parameter_list = list(combined_T.keys())
# kcat_matrix = np.array(list(combined_T.values()))

# min_nonzero_fraction = 0.5  # Require at least 50% of points to be non-zero
# filter_threshold = 1e-3  # log10(kcat + 1) < 1e-3 ≈ kcat ≈ 0

# filtered_parameter_list = []
# filtered_kcat_matrix = []
# min_logkcat_temp3 = 4  # Because 10^4 = 10,000 ⇒ kcat > 10,000


# for i, row in enumerate(kcat_matrix):
#     nonzero_fraction = np.sum(row > filter_threshold) / len(row)
#     high_kcat_at_temp3 = row[3] > min_logkcat_temp3  # index 3 = 4th temperature

#     if nonzero_fraction >= min_nonzero_fraction and high_kcat_at_temp3:
#         filtered_parameter_list.append(parameter_list[i])
#         filtered_kcat_matrix.append(row)
    

# filtered_kcat_matrix = np.array(filtered_kcat_matrix)
# logging.info(f"Filtered down to {len(filtered_kcat_matrix)} parameter sets with sufficient non-zero kcat values.")

# kcat_matrix = filtered_kcat_matrix
# parameter_list = filtered_parameter_list

# kcat_matrix_scaled = StandardScaler().fit_transform(kcat_matrix)  # Transpose for scaling
# logging.info(f"Shape of kcat_matrix_scaled: {kcat_matrix_scaled.shape}")

# # Define a similarity threshold (tune as needed)
# threshold = 0.1  # very tight tolerance

# # Find all parameter sets with nearly identical kcat profiles
# logging.info(f"Finding parameter sets with kcat profiles within {threshold} of each other.")
# start = time.time()
# matches = []
# for i in range(len(kcat_matrix)):
#     if i % 10000 == 0:
#         logging.info(f"Processing parameter set {i} to find similar kcat profiles.")
#     dists = euclidean_distances(kcat_matrix[i:i+1], kcat_matrix)
#     similar = np.where(dists[0] < threshold)[0]
#     if len(similar) > 1:
#         matches.append(similar)
# end = time.time()
# logging.info(f"Found {len(matches)} sets of parameters with nearly identical kcat profiles in {end-start} seconds.")

# dump_pickle(matches, f'../results/analysis/{enzyme_id}_matches_treshold{threshold}.pkl')
# matches = load_pickle(f'../results/analysis/{enzyme_id}_matches_treshold{threshold}.pkl')

# logging.info(f"Merging matches into clusters")
# start = time.time()
# clusters = merge_matches(matches)
# end = time.time()
# logging.info(f"Found {len(clusters)} clusters of similar kcat profiles in {end-start} seconds.")

# diversity_stats = [analyze_cluster_diversity(cluster, parameter_list) for cluster in clusters]

# logging.info("Diversity statistics for each cluster:")
# diverse_clusters = sorted(
#     enumerate(diversity_stats),
#     key=lambda x: (x[1]["Tm_range"] + x[1]["Topt_range"] + x[1]["dCpt_range"]),
#     reverse=True
# )[:5]

# for idx, stats in diverse_clusters:
#     logging.info(f"\nCluster {idx} with {stats['size']} members:")
#     logging.info(stats)
#     logging.info("Input parameter sets:")
#     for i in clusters[idx]:
#         logging.info(parameter_list[i])


#make_diversity_plot(43, sel_temp, max_legend_entries=10)  # Example for cluster 2

# logging.info("Perform clustering using DBSCAN")
# db = DBSCAN(eps=0.8, min_samples=5, metric='euclidean')
# labels = db.fit_predict(kcat_matrix_scaled)

# clustered_params = defaultdict(list)
# for label, params in zip(labels, parameter_list):
#     clustered_params[label].append(params)

# logging.info(f"Number of clusters found: {len(set(labels)) - (1 if -1 in labels else 0)}")
# logging.info(f"Cluster 0 contains {len(clustered_params[0])} parameters.")
# logging.info(clustered_params[0][:5])  # Log first 5 parameters in cluster 0 for brevity



# -------------------------------------------------------------------------------------------------------------------
logging.info("Analyzing kcat profiles from evolutionary data")
df = load_pickle("../results/analysis/evo_combined_df_R098.pkl")
logging.info(f"DataFrame loaded with shape: {df.shape}")
df_tm_values = df[f'{enzyme_id}_Tm'].values
df_topt_values = df[f'{enzyme_id}_Topt'].values
df_dCpt_values = df[f'{enzyme_id}_dCpt'].values

candidates = [(tm, topt, dCpt) for tm, topt, dCpt in zip(df_tm_values, df_topt_values, df_dCpt_values) if tm >= topt]
logging.info(f"Total candidates found: {len(candidates)}")

unique_candidates = []  # Remove duplicates
for candidate in candidates:
    if candidate not in unique_candidates:
        unique_candidates.append(candidate)

candidates = unique_candidates  # Use the unique candidates list


logging.info(f"Total unique candidates found: {len(candidates)}")

logging.info(f"Tm values: min = {min(tm_values)}, max = {max(tm_values)}")
logging.info(f"Topt values: min = {min(topt_values)}, max = {max(topt_values)}")
logging.info(f"dCpt values: min = {min(dCpt_values)}, max = {max(dCpt_values)}")
logging.info(f"First 5 candidates: {candidates[:5]}")

logging.info(f"Calculating kcat for {len(candidates)} candidates at selected temperatures: {sel_temp}")
all_result_dicts = []
for T in sel_temp:
    evo_result_dict = {}  # Initialize result dictionary for storing results
    for tm, topt, dCpt in candidates:
        if tm < topt:
            logging.warning(f"Skipping candidate with Tm < Topt: Tm={tm}, Topt={topt}, dCpt={dCpt}")
            continue
        # Call the function with the current parameters
        func2(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase, T_dict=evo_result_dict)
    all_result_dicts.append(evo_result_dict)

combined_T = {key: [] for key in all_result_dicts[0].keys()}
for d in all_result_dicts:
    for params, kcat in d.items():
        log_kcat = np.log10(kcat+1)  
        combined_T[params].append(kcat)


make_diversity_plot_from_simdata(result_dict=combined_T, sel_temp=sel_temp, max_legend_entries=10)









#--------------------------------------------------------------------------------------------------------------------
# Load the model and parameters
#param_dict = format_input(params, model_particle)
# topt_values = [309.3333333333333]
# dCpt_values = [-15000.0]
tm_values = np.linspace(312, 346, 20)  # Example Tm values in K
topt_values = np.linspace(272, 328, 20)  # Example Tm values in K
dCpt_values = np.linspace(-15000, 1000, 20)  # Example dCpt values in J/mol/K

# test_tm = tm_values[10]
# test_opt = topt_values[10]
# logging.info(f"Test Tm: {test_tm}, topt: {test_opt}")

all_T_dicts = []
#For tm
# for T in sel_temp:
#     T_dict: Dict[float, List[float]] = {}
#     logging.info(f"Processing temperature: {T} K")
#     for tm in tm_values:
#         for dCpt, topt in product([-15000.0], [309.3333333333333]):
#             if tm < topt: continue
#             func(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase, T_dict=T_dict)
#     logging.info(f"Length of T_dict: {len(T_dict)}")
#     logging.info(f"T_dict[:10]: {list(T_dict.items())[:10]}")
#     all_T_dicts.append(T_dict)
#     #make_plot(T_dict=T_dict, T=T, enzyme_id=enzyme_id, reaction=shikimate_kinase)

# #For topt
# for T in sel_temp:
#     T_dict: Dict[float, List[float]] = {}
#     logging.info(f"Processing temperature: {T} K")
#     for topt in topt_values:
#         for dCpt, tm in product([-7222.222222222222], [329.89473684210526]):
#             if tm < topt: continue
#             func(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase, T_dict=T_dict)
#     logging.info(f"Length of T_dict: {len(T_dict)}")
#     logging.info(f"T_dict[:10]: {list(T_dict.items())[:10]}")
#     all_T_dicts.append(T_dict)
#     #make_plot(T_dict=T_dict, T=T, enzyme_id=enzyme_id, reaction=shikimate_kinase)

# #for dCpt
# for T in sel_temp:
#     T_dict: Dict[float, List[float]] = {}
#     logging.info(f"Processing temperature: {T} K")
#     for dCpt in dCpt_values:
#         for tm, topt in product([329.89473684210526], [301.4736842105263]):
#             if tm < topt: continue
#             func(tm=tm, topt=topt, dCpt=dCpt, enzyme_id=enzyme_id, T=T, reaction=shikimate_kinase, T_dict=T_dict)
#     logging.info(f"Length of T_dict: {len(T_dict)}")
#     logging.info(f"T_dict[:10]: {list(T_dict.items())[:10]}")
#     all_T_dicts.append(T_dict)

#plot_combined_kcats(all_T_dicts, sel_temp, enzyme_id, shikimate_kinase)



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