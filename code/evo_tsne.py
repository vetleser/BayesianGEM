#!/usr/bin/env python
# coding: utf-8

# This script is similar to evo_pca.py and has a direct dependence on it, but does t-SNE instead of PCA 
import pickle
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import logging

RANDOM_STATE = 5902


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def perform_tsne_on_parameters(df):
    # 1. normalize all columns to a standard normal distribution
    X = df.values[:,:-3]
    X_n = np.zeros_like(X)    
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/np.std(X[:,i])
    tsne = TSNE(2, learning_rate="auto", random_state=RANDOM_STATE)
    return tsne.fit_transform(X_n)


combined_df = load_pickle("../results/evo_combined_particle_df.pkl")
tsne_ordination = perform_tsne_on_parameters(combined_df)
dump_pickle(tsne_ordination,"../results/evo_tsne_full_ordination.pkl")
logging.info("DONE")
