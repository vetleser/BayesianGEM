#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import logging
from evo_etc import CrowdingDE
from sa_etc import SimulatedAnnealing

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_posterior_particles(file, r2_threshold = 0.98):
    results: SimulatedAnnealing = load_pickle(file)
    columns = list(results.all_particles[0].keys())
    columns.sort()
    logging.info("Iterating over particles")
    data = list()
    for p in results.all_particles:
        data.append([p[k] for k in columns])
    logging.info("Creating Data Frame")
    df = pd.DataFrame(data=data,columns=columns)
    df['r2'] = results.all_distances
    # Running number assigned to particles to keep track of them when comparing with original data
    df['particle_ID'] = list(range(len(results.all_particles)))
    logging.info(df.shape)
    
    logging.info("Doing filtering and labelling of Data Frame")
    # We need to negate the results due to the fact that they
    # are orignally taken to mean distances which are to be minimized
    df['r2'] = -df['r2']
    df = df[df['r2'] > r2_threshold]
    logging.info(df.shape)
    return df

def combine_dataframes(df_dict):
    for label, df in df_dict.items():
        df["frame_ID"] = label
    return pd.concat(df_dict.values(),ignore_index=True)


outdir = '../results/sa'
model_frame: pd.DataFrame = load_pickle(f"{outdir}/simulation_skeleton.pkl")
model_frame["frame_ID"] = range(model_frame.shape[0])
# We are in this case only interested in the solutions with F=0.5 and CR=0.99
reduced_model_frame = (
    model_frame.
    set_index(["end"]).
    loc[(0.25, 0.5, 0.75), :].
)
logging.info("Loading data")
particle_dfs = list(map(build_a_dataframe_for_posterior_particles,reduced_model_frame.outfile))
logging.info("Augmenting data labeling")
df_dict = {}
for i, particle_df in zip(reduced_model_frame["frame_ID"],particle_dfs):
    df_dict[i] = particle_df

logging.info("Combining dataframes")
combined_df = combine_dataframes(df_dict)
dump_pickle(combined_df,f"{outdir}/evo_combined_df_R0985.pkl")