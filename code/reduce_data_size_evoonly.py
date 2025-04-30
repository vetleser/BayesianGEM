#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import logging

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

def extract_distances_and_population_from_simulation(filename):
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances, sim_res.population


evo_simulation_skeleton = load_pickle("../results/crowdingDE/simulation_skeleton.pkl")

#New code:
desired_files = [
    "../results/crowdingDE/smcevo_gem_1.0_0.999_0.pkl",
    "../results/crowdingDE/smcevo_gem_1.0_0.999_1.pkl"
]

evo_simulation_skeleton = evo_simulation_skeleton[evo_simulation_skeleton["outfile"].isin(desired_files)]
#New code ended



evo_simulation_skeleton["all_distances"], evo_simulation_skeleton["population"] = zip(*list(map(extract_distances_and_population_from_simulation, evo_simulation_skeleton["outfile"])))

dump_pickle(evo_simulation_skeleton, "../results/crowdingDE/distance_frame.pkl")

def aggregate_fva_results(result_df,simulation_attributes):
    flattened_df_list = []
    for _, row in result_df.drop(columns=["particle"]).iterrows():
        raw_df = row["fva_res"]
        df = raw_df
        for attribute in simulation_attributes:
            df[attribute] = row[attribute]
        flattened_df_list.append(df)

    combined_fva_frame = (
        pd.concat(flattened_df_list).
        assign(range= lambda df: df["maximum"] - df["minimum"],
                                                            midpoint= lambda df: (df["maximum"] + df["minimum"]) / 2).
        drop(columns=["minimum", "maximum"])
            )
    simulation_attributes.extend(["condition","reaction","T"])
    aggregated_fva_res = (
        combined_fva_frame.replace([np.inf, -np.inf],np.nan).
        dropna(how="all").
        groupby(simulation_attributes).
        agg(["mean","min","max","std","count"])
                        )
    return aggregated_fva_res


evo_fva_results = load_pickle("../results/crowdingDE/evo_fva.pkl")
evo_aggregated_fva_results = aggregate_fva_results(evo_fva_results,["scaling_factor","crossover_prob","simulation"])
dump_pickle(evo_aggregated_fva_results, "../results/crowdingDE/evo_aggregated_fva_res.pkl")
