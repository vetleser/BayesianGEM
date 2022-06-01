#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import logging

# This scripts depends on sample_pca.py, but only consideres the Unpermuted and Permuted 1 results

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))


def perform_pca_on_parameters(df):
    # 1. normalize all columns to a standard normal distribution
    X = df.values[:,:-3]
    X_n = np.zeros_like(X)    
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/np.std(X[:,i])
    pca = PCA(n_components=2)
    PCS = pca.fit_transform(X_n)
    logging.info(pca.explained_variance_ratio_)
    return PCS, pca.explained_variance_ratio_


combined_df = load_pickle("../results/permuted_smcabc_res/combined_particle_df.pkl")
reduced_df = combined_df[combined_df["origin"] in ("unpermuted","permuted_0")]
pca_ordination = perform_pca_on_parameters(reduced_df)
dump_pickle(pca_ordination,"../results/permuted_smcabc_res/pca_reduced_ordination.pkl")
