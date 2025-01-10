#!/usr/bin/env python
# coding: utf-8

# import pickle
# import pandas as pd
# import numpy as np
# from sklearn.decomposition import PCA
# import logging
# from evo_etc import CrowdingDE

import pickle
import pandas as pd
import numpy as np
import logging

from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from sklearn.decomposition import PCA


import seaborn as sns
import umap

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def perform_pca_on_parameters(df):
    epsilon = 1e-6

    logging.info(f"Shape of df is {df.shape}")

    # 1. normalize all columns to a standard normal distribution
    # df contains the trailing columns: r2, period, particle_ID, frame_ID, in total 4 columns to remove
    X = df.values[:,:-4]
    X_n = np.zeros_like(X)

    logging.info(f"Shape of X is {X.shape}")
    logging.info(f"Shape of X_n is {X_n.shape}")
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/(np.std(X[:,i]) + epsilon)
    pca = PCA(n_components=2)
    PCS = pca.fit_transform(X_n)
    logging.info(f"PCS shape is {PCS.shape}")
    logging.info(pca.explained_variance_ratio_)
    return PCS, pca.explained_variance_ratio_

#The above functions are taken from evo_pca.py

def perform_umap_on_parameters(df, n_components=2, n_neighbors=15, min_dist=0.1):
    epsilon = 1e-6
    logging.info(f"Shape of df is {df.shape}")
    X = df.values[:,:-4]
    #X = df.values[:10,:-4] #10 chooses only the first 10 rows for testing purposes
    X_n = np.zeros_like(X)
    logging.info(f"Shape of X is {X.shape}")
    logging.info(f"Shape of X_n is {X_n.shape}")

    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/(np.std(X[:,i]) + epsilon)

    reducer = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist)

    embedding = reducer.fit_transform(X_n)

    logging.info(f"Embedding shape is {embedding.shape}")

    return embedding




outdir = '../results/crowdingDE'
combined_df = load_pickle(f"{outdir}/evo_combined_df.pkl")

# logging.info("Performing PCA")
# pca_ordination = perform_pca_on_parameters(combined_df)

# logging.info("Performing UMAP")

n_components=2 #Dimensions of resulting data set, keep as 2
n_neighbors=200 #Can change
min_dist=0.1   #Can change

umap_ordination = perform_umap_on_parameters(df=combined_df, 
                                             n_components = n_components, 
                                             n_neighbors=n_neighbors, 
                                             min_dist=min_dist)
dump_pickle(umap_ordination,f"{outdir}/evo_umap.pkl")

logging.info("DONE")
