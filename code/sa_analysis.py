#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import multiprocessing
import logging
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_all_particles(file, n_priors = 128, r2_threshold = 0.9):
    results = load_pickle(file)
    columns = list(results.all_particles[0].keys())
    columns.sort()
    logging.info("Iterating over particles")
    data = list()
    for p in results.all_particles:
        data.append([p[k] for k in columns])
    logging.info("Creating Data Frame")
    df = pd.DataFrame(data=data,columns=columns)
    df['r2'] = results.all_distances
    logging.info(df.shape)
    
    
    logging.info("Doing filtering and labelling of Data Frame")
    df['r2'] = -df['r2']
    df["period"] = "Intermediate"
    df.loc[:n_priors,"period"] = "Prior"
    df.loc[df["r2"] > r2_threshold,"period"] = 'Posterior'
    # Remove samples with a R2 score smaller than -3
    sel_index = df.index[df['r2']>-3]    
    df = df.loc[sel_index,:]
    logging.info(df.shape)

    return df

def inspect_acceptance_rate(filename, n_priors = 128):
    logging.info(f"Loading file: {filename}")
    results = load_pickle(filename)
    acceptance_rates = results.acceptance_rates

    return acceptance_rates 

outdir = '../results/sa'

file_0 = f'{outdir}/smcsa_gem_0.001_0.pkl'
file_1 = "../results/sa/smcsa_gem_may19_0.0001_1.pkl"
file_2 = f'{outdir}/smcsa_gem_may19_0.0001_0.pkl'

# df_0 = build_a_dataframe_for_all_particles(file_0)

# logging.info(f"Data Frame 0 shape: {df_0.shape}")
# logging.info(f"Data Frame 0 columns: {df_0.columns}")
# logging.info(f"Data Frame 0 head: {df_0.head()}")
# logging.info(f"Data Frame 0 tail: {df_0.tail()}")

acceptance_rates_0 = inspect_acceptance_rate(file_0)
acceptance_rates_1 = inspect_acceptance_rate(file_1)
acceptance_rates_2 = inspect_acceptance_rate(file_2)

plt.figure(figsize=(10, 5))
plt.plot(acceptance_rates_0, label='Acceptance Rate 0')
plt.plot(acceptance_rates_1, label='Acceptance Rate 1')
plt.plot(acceptance_rates_2, label='Acceptance Rate 2')
plt.xlabel('Generation')
plt.ylabel('Acceptance Rate')
plt.title('Acceptance Rate over all Generations')
plt.legend()
plt.show()
plt.savefig('../figures/acceptance_rate.png')
#     df_2 = build_a_dataframe_for_all_particles(file_2)

logging.info("DONE")

