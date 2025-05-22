import os
import copy
import time
import math
import pickle
import random
import logging

import numpy as np
import pandas as pd
import numpy.typing as npt
from numpy.typing import NDArray
from typing import Dict, Iterable, List, Optional
from functools import partial
from itertools import starmap

from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import r2_score
from scipy.optimize import minimize

import evo_etc as evo
import evo_etc as CrowdingDE
import GEMS
from random_sampler import RV
from etcpy import etc
from etcpy.thermal_parameters import format_input
import multiprocessing
import reframed


# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))


outdir = "../results/analysis/flux_analysis"
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")

protein_IDs = pd.read_csv("../data/model_enzyme_params.csv").iloc[:,0].tolist()
params = pd.read_csv(('../data/model_enzyme_params.csv'),index_col=0)


df_flux = load_pickle(f"{outdir}/combined_flux_data_120677.0.pkl")
logging.info(f"Loaded df_r: \n {df_flux}")
logging.info(f"df_r shape: {df_flux.shape}")
logging.info(f"df_r columns: {df_flux.columns[:10]}")
#filtered_df = df_r[df_r.index.str.match(r"r_\d{4}No1")]

#param_dict = format_input(params, param_dict=df_flux.to_dict(orient='index'))


df_flux_importance = (df_flux > 0).sum(axis=0) / len(df_flux)
df_flux_importance = df_flux_importance.to_frame().T  # Make it a single-row DataFrame
#df_flux_importance = df_flux_importance.loc[:, df_flux_importance.ne(0.0).all()]


logging.info(f"df_flux_importance:\n{df_flux_importance}")

ae_model = load_pickle(f"../models/aerobic.pkl")
an_model = load_pickle(f"../models/anaerobic.pkl")

logging.info(f"Loaded ae_model: {type(ae_model)}")
# logging.info(f"Loaded an_model: {an_model}")

# logging.info(f"Reactions in ae_model: {len(reactions)}")
# logging.info(f"Loaded ReFramed ae_model: {ae_model}")
logging.info(f"Reactions in ae_model: {len(ae_model.reactions)}")

def convert_to_uniprot(gene: str) -> str:
    """
    Convert a gene ID to a UniProt ID.
    This is a placeholder function. The actual conversion logic should be implemented here.
    """
    # Example conversion logic (this should be replaced with actual logic)
    if gene.startswith("gene_"):
        return gene.replace("gene_", "uniprot_")
    return gene

reactions: set = set()
ids = set()
param_importance : Dict[str, float] = {ID: 0.0 for ID in protein_IDs}
logging.info(f"Parameter importance length: {len(param_importance)}")

for rxn_id, rxn in ae_model.reactions.items():
    genes: set = rxn.get_genes()
    if len(genes) > 0:
        uniprot_id = [convert_to_uniprot(gene) for gene in genes]
        #logging.info(f"Reaction {rxn_id} is associated with genes: {genes}")
    for met in rxn.stoichiometry:
        if not met.startswith('prot_'): continue
            # ingore metabolite: prot_pool
        if met == 'prot_pool': continue
        uniprot_id = met.split('_')[1]
        if uniprot_id in protein_IDs:
            parameter_entries = params.loc[uniprot_id]
            #logging.info(f"Reaction {rxn_id} is associated with uniprot_id: {uniprot_id}")
            reactions.add(rxn_id)
            ids.add(uniprot_id)
            param_importance[uniprot_id] += df_flux_importance[rxn_id].values[0]
   

logging.info(f"Enzyme calyzed reactions in ae_model: {len(reactions)}")
logging.info(f"Enzyme IDs in ae_model: {len(ids)}")

#logging.info(f"Parameter importance: {param_importance}")
non_zero_importance = {k: v for k, v in param_importance.items() if v > 0}
logging.info(f"Number of non-zero values in param_importance: {len(non_zero_importance)}")

sorted_param_importance = dict(sorted(param_importance.items(), key=lambda item: item[1], reverse=True))

logging.info(f"Sorted param importance: {list(sorted_param_importance)[:10]}")





# for col in df_r.columns:
#     this column is a reaction flux
#     need to check which enzyme catalyzes this reaction
#     convert the reaction importance to enzyme importance
