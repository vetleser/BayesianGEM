#!/usr/bin/env python
# coding: utf-8

# This script is similar to evo_pca.py and has a direct dependence on it, but does t-SNE instead of PCA
# t-SNE is repeated on 12 different perplexities which results are intended to be compared afterwards. 
import multiprocessing
import pickle
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from functools import partial
import logging


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def perform_tsne_on_parameters(df, perplexity):
    # 1. normalize all columns to a standard normal distribution
    X = df.values[:,:-3]
    X_n = np.zeros_like(X)    
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/np.std(X[:,i])
    tsne = TSNE(2, perplexity=perplexity, learning_rate="auto", init="pca")
    return tsne.fit_transform(X_n)

N_CORES = multiprocessing.cpu_count()
cpu_pool = multiprocessing(processes=N_CORES)
perplexities = [5, 10, 30, 60, 100, 200, 300, 500, 1000, 2000, 5000, 10000]
combined_df = load_pickle("../results/evo_combined_particle_df.pkl")
tsne_frame = pd.DataFrame({"perplexity": perplexities})

tsne_frame["ordination"] = cpu_pool.map(partial(perform_tsne_on_parameters, combined_df),perplexities)
dump_pickle(tsne_frame,"../results/evo_tsne_full_ordination.pkl")
logging.info("DONE")