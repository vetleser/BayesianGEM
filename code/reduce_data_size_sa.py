#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import logging
import re

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# The purpose of this script is to trim the computational results from the simulations and subsequent PCA ordinations
# to the bare minimum required in the Jupyter notebook analyzing and visualizing the results.
# This is to make it fast and practical to run the notebook with laptop-grad hardware

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))


def extract_distances_from_simulation(filename):
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances

def extract_distances_and_population_and_acceptance_rates_from_simulation(filename):
    logging.info(f"Loading file: {filename}")
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances, sim_res.population, sim_res.acceptance_rates


sa_simulation_skeleton = load_pickle("../results/sa/simulation_skeleton.pkl")
columns = sa_simulation_skeleton.columns.tolist()

logging.info(f"Columns in the skeleton: {columns}")
#Columns in the skeleton: ['final_temp', 'move_type', 'step_size', 'simulation', 'outfile', 'random_seed']
#Newest code
outdir = '../results/sa'



file_1 = f'{outdir}/smcsa_gem_june9_ee0.25_normalize_0.pkl'
file_2 = f'{outdir}/smcsa_gem_june9_ee0.25_normalize_1.pkl'
file_3 = f'{outdir}/smcsa_gem_june9_ee0.5_normalize_2.pkl'
file_4 = f'{outdir}/smcsa_gem_june4_ee0.5_normalize_2.pkl'
file_5 = f'{outdir}/smcsa_gem_june9_ee0.75_normalize_3.pkl'
file_6 = f'{outdir}/smcsa_gem_june5_ee0.75_normalize_3.pkl'

desired_files = [
    file_1,
    file_2,
    file_3,
    file_4,
    file_5,
    file_6
]

data = [
    {'end_exploration': 0.25, 'initial_step_size': 1.0, 'final_step_size': 0.1, 'simulation' : 0, 'random_seed': 1234, 'outfile': file_1},
    {'end_exploration': 0.25, 'initial_step_size': 1.0, 'final_step_size': 0.1, 'simulation' : 1, 'random_seed': 10199, 'outfile': file_2},

    {'end_exploration': 0.5,  'initial_step_size': 1.0, 'final_step_size': 0.1, 'simulation' : 0, 'random_seed': 19164, 'outfile': file_3},
    {'end_exploration': 0.5,  'initial_step_size': 1.0, 'final_step_size': 0.1, 'simulation' : 1, 'random_seed': 23583, 'outfile': file_4},
    {'end_exploration': 0.75, 'initial_step_size': 1.0, 'final_step_size': 0.1, 'simulation' : 0, 'random_seed': 28129, 'outfile': file_5},
    {'end_exploration': 0.75, 'initial_step_size': 1.0, 'final_step_size': 0.1, 'simulation' : 1, 'random_seed': 66391, 'outfile': file_6}
]


df = pd.DataFrame(data)[['end_exploration', 'initial_step_size', 'final_step_size', 'simulation', 'random_seed', 'outfile']]
sa_simulation_skeleton = df.sort_values(by='end_exploration').reset_index(drop=True)



logging.info(f"DataFrame with desired files:\n{sa_simulation_skeleton}")
#newest code ended

# #New code:
# desired_files = [
#     "../results/sa/smcsa_gem_0.001_0.pkl",
#     "../results/sa/smcsa_gem_may19_0.0001_0.pkl",
#     "../results/sa/smcsa_gem_may19_0.0001_1.pkl",
#     "../results/sa/smcsa_gem_may22_0.1_0.pkl",
#     "../results/sa/smcsa_gem_may22_0.5_0.pkl",
#     "../results/sa/smcsa_gem_may22_0.5_1.pkl", 
# ]




# # Original filenames and data
# data = [
#     {'final_temp': 0.001,  'step_size': 0.1, 'move_type': 'gaussian', 'outfile': '../results/sa/smcsa_gem_0.001_0.pkl'},
#     {'final_temp': 0.0001, 'step_size': 1.0, 'move_type': 'normal',   'outfile': '../results/sa/smcsa_gem_may19_0.0001_1.pkl'},
#     {'final_temp': 0.0001, 'step_size': 1.0, 'move_type': 'gaussian', 'outfile': '../results/sa/smcsa_gem_may19_0.0001_0.pkl'},
#     {'final_temp': 0.0001, 'step_size': 0.1, 'move_type': 'normal',   'outfile': '../results/sa/smcsa_gem_may22_0.1_0.pkl'},
#     {'final_temp': 0.0001, 'step_size': 0.5, 'move_type': 'normal',   'outfile': '../results/sa/smcsa_gem_may22_0.5_0.pkl'},
#     {'final_temp': 0.0001, 'step_size': 0.5, 'move_type': 'gaussian', 'outfile': '../results/sa/smcsa_gem_may22_0.5_1.pkl'}
# ]

# # Extract simulation number from the filename
# for row in data:
#     match = re.search(r'_(\d+)\.pkl$', row['outfile'])
#     row['simulation'] = int(match.group(1)) if match else np.nan

# # Reorder columns
# df = pd.DataFrame(data)[['final_temp', 'move_type', 'step_size', 'simulation', 'outfile']]
# sa_simulation_skeleton = df.sort_values(by='step_size').reset_index(drop=True)


# logging.info(f"DataFrame with desired files:\n{sa_simulation_skeleton}")
# # evo_simulation_skeleton = evo_simulation_skeleton[evo_simulation_skeleton["outfile"].isin(desired_files)]
# # #New code ended




sa_simulation_skeleton["all_distances"], sa_simulation_skeleton["population"], sa_simulation_skeleton["acceptance_rates"] = zip(*list(map(extract_distances_and_population_and_acceptance_rates_from_simulation, sa_simulation_skeleton["outfile"])))

dump_pickle(sa_simulation_skeleton, "../results/sa/distance_frame_final.pkl")

logging.info(f"DataFrame with distances and population:\n{sa_simulation_skeleton}")
logging.info(f"Number of simulations: {len(sa_simulation_skeleton)}")
logging.info(f"Columns in the DataFrame: {sa_simulation_skeleton.columns.tolist()}")

# def aggregate_fva_results(result_df,simulation_attributes):
#     flattened_df_list = []
#     for _, row in result_df.drop(columns=["particle"]).iterrows():
#         raw_df = row["fva_res"]
#         df = raw_df
#         for attribute in simulation_attributes:
#             df[attribute] = row[attribute]
#         flattened_df_list.append(df)

#     combined_fva_frame = (
#         pd.concat(flattened_df_list).
#         assign(range= lambda df: df["maximum"] - df["minimum"],
#                                                             midpoint= lambda df: (df["maximum"] + df["minimum"]) / 2).
#         drop(columns=["minimum", "maximum"])
#             )
#     simulation_attributes.extend(["condition","reaction","T"])
#     aggregated_fva_res = (
#         combined_fva_frame.replace([np.inf, -np.inf],np.nan).
#         dropna(how="all").
#         groupby(simulation_attributes).
#         agg(["mean","min","max","std","count"])
#                         )
#     return aggregated_fva_res


# evo_fva_results = load_pickle("../results/crowdingDE/evo_fva.pkl")
# evo_aggregated_fva_results = aggregate_fva_results(evo_fva_results,["scaling_factor","crossover_prob","simulation"])
# dump_pickle(evo_aggregated_fva_results, "../results/crowdingDE/evo_aggregated_fva_res.pkl")

logging.info("DONE")
