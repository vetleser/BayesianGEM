#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import logging
import evo_etc as CrowdingDE

#Files to look at:
#evo_pca, reduce_data_size

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

outdir_write_to = "../results/analysis"

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))



#Begin here
def build_a_dataframe_for_posterior_particles(file, r2_threshold = 0.98):
    results: CrowdingDE = load_pickle(file)
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
    logging.info(f"Shape of {file} before filtering: {df.shape}")
    
    logging.info("Doing filtering and labelling of Data Frame")
    # We need to negate the results due to the fact that they
    # are orignally taken to mean distances which are to be minimized
    df['r2'] = -df['r2']
    df = df[df['r2'] > r2_threshold]
    logging.info(f"Shape of {file} after filtering: {df.shape}")

    return df

def combine_dataframes(df_dict):
    for label, df in df_dict.items():
        df["frame_ID"] = label
    return pd.concat(df_dict.values(),ignore_index=True)



outdir_read_from = '../results/crowdingDE'
model_frame: pd.DataFrame = load_pickle(f"{outdir_read_from}/simulation_skeleton.pkl")
model_frame["frame_ID"] = range(model_frame.shape[0])
logging.info("Loading data")
particle_dfs = list(map(build_a_dataframe_for_posterior_particles,model_frame.outfile))
logging.info("Augmenting data labeling")
df_dict = {}
for i, particle_df in zip(model_frame["frame_ID"],particle_dfs):
    dump_pickle(particle_df, f"{outdir_write_to}/df_simulation_{i}_R098.pkl")
    df_dict[i] = particle_df


logging.info("Combining dataframes")
combined_df = combine_dataframes(df_dict)
dump_pickle(combined_df,f"{outdir_write_to}/evo_combined_df_R098.pkl")
logging.info("filtering_data.py DONE")
