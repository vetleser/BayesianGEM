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

import seaborn as sns
import umap

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

#The above functions are taken from evo_pca.py

def perform_umap_on_parameters(df, reducer):
    epsilon = 1e-6
    X = df.values[:10,:-4] #10 chooses only the first 10 rows for testing purposes
    X_n = np.zeros_like(X)
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/(np.std(X[:,i]) + epsilon)

    embedding = reducer.fit_transform(X_n)

    logging.info(embedding.shape)

    return embedding, embedding




outdir = '../results/crowdingDE'
combined_df = load_pickle(f"{outdir}/evo_combined_df.pkl")

logging.info("Performing UMAP")

reducer = umap.UMAP()

umap_ordination = perform_umap_on_parameters(combined_df, reducer)
dump_pickle(umap_ordination,f"{outdir}/evo_umap.pkl")

logging.info("DONE")
