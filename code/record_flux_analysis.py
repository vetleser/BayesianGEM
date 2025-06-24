from collections import Counter
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
import matplotlib.pyplot as plt

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


logging.info(f"Loading relevant enzyme parameters models")
protein_IDs = pd.read_csv("../data/model_enzyme_params.csv").iloc[:,0].tolist()
params = pd.read_csv(('../data/model_enzyme_params.csv'),index_col=0)
ae_model = load_pickle(f"../models/aerobic.pkl")
logging.info(f"Loaded ae_model: {type(ae_model)}")
# logging.info(f"Loaded ae_model: {ae_model}")
an_model = load_pickle(f"../models/anaerobic.pkl")
# logging.info(f"Loaded an_model: {type(an_model)}")
# logging.info(f"Loaded an_model: {an_model}")
# for rxn_id, rxn in ae_model.reactions.items():
#     if rxn_id not in an_model.reactions:
#         logging.info(f"Reaction {rxn_id} in ae_model not in an_model: {rxn}")
#     elif str(rxn) != str(an_model.reactions[rxn_id]):
#         logging.info(f"Reaction {rxn_id} differs between models:")
#         logging.info(f"  ae_model: {rxn}")
#         logging.info(f"  an_model: {an_model.reactions[rxn_id]}")

logging.info(f"Finished loading models")

particle_IDs_to_evaluate = load_pickle(f"{outdir}/particle_IDs_to_evaluate.pkl")
logging.info(f"Loaded particle_IDs_to_evaluate: {particle_IDs_to_evaluate}")




def get_combined_parameter_importance(particle_ID):
    logging.info(f"Loading particle_ID: {particle_ID}")
    df_flux = load_pickle(f"{outdir}/combined_flux_data_{particle_ID}.pkl")
    df_flux_importance = (df_flux > 0).sum(axis=0) / len(df_flux)
    df_flux_importance = df_flux_importance.to_frame().T
    param_importance : Dict[str, float] = {ID: 0.0 for ID in protein_IDs}

    # for model in [ae_model, an_model]:
    #     if model == ae_model:
    #         logging.info(f"Loading ae_model")
    #     else:
    #         logging.info(f"Loading an_model")
    for rxn_id, rxn in ae_model.reactions.items():
        if not rxn_id.startswith('draw_'): continue
        for met in rxn.stoichiometry:
            if not met.startswith('prot_'): continue
                # ingore metabolite: prot_pool
            if met == 'prot_pool': continue
            uniprot_id = met.split('_')[1]
            #logging.info(f"Reaction {rxn_id}: {rxn} is associated with uniprot_id: {uniprot_id}")
            #reactions.add(rxn_id)
            #ids.add(uniprot_id)
            param_importance[uniprot_id] += df_flux_importance[rxn_id].values[0]
    return param_importance

def get_ae_parameter_importance(particle_ID):
    logging.info(f"Loading particle_ID: {particle_ID}")
    df_flux = load_pickle(f"{outdir}/combined_flux_data_{particle_ID}.pkl")
    df_flux = df_flux.iloc[0:8]
    df_flux_importance = (df_flux > 0).sum(axis=0) / len(df_flux)
    df_flux_importance = df_flux_importance.to_frame().T
    param_importance : Dict[str, float] = {ID: 0.0 for ID in protein_IDs}

    # for model in [ae_model, an_model]:
    #     if model == ae_model:
    #         logging.info(f"Loading ae_model")
    #     else:
    #         logging.info(f"Loading an_model")
    for rxn_id, rxn in ae_model.reactions.items():
        if not rxn_id.startswith('draw_'): continue
        for met in rxn.stoichiometry:
            if not met.startswith('prot_'): continue
                # ingore metabolite: prot_pool
            if met == 'prot_pool': continue
            uniprot_id = met.split('_')[1]
            #logging.info(f"Reaction {rxn_id}: {rxn} is associated with uniprot_id: {uniprot_id}")
            #reactions.add(rxn_id)
            #ids.add(uniprot_id)
            param_importance[uniprot_id] += df_flux_importance[rxn_id].values[0]
    return param_importance

def get_an_parameter_importance(particle_ID):
    logging.info(f"Loading particle_ID: {particle_ID}")
    df_flux = load_pickle(f"{outdir}/combined_flux_data_{particle_ID}.pkl")
    df_flux = df_flux.iloc[8:16]
    df_flux_importance = (df_flux > 0).sum(axis=0) / len(df_flux)
    df_flux_importance = df_flux_importance.to_frame().T
    param_importance : Dict[str, float] = {ID: 0.0 for ID in protein_IDs}

    # for model in [ae_model, an_model]:
    #     if model == ae_model:
    #         logging.info(f"Loading ae_model")
    #     else:
    #         logging.info(f"Loading an_model")
    for rxn_id, rxn in ae_model.reactions.items():
        if not rxn_id.startswith('draw_'): continue
        for met in rxn.stoichiometry:
            if not met.startswith('prot_'): continue
                # ingore metabolite: prot_pool
            if met == 'prot_pool': continue
            uniprot_id = met.split('_')[1]
            #logging.info(f"Reaction {rxn_id}: {rxn} is associated with uniprot_id: {uniprot_id}")
            #reactions.add(rxn_id)
            #ids.add(uniprot_id)
            param_importance[uniprot_id] += df_flux_importance[rxn_id].values[0]
    return param_importance


zero_importance_enzymes_list = []
combined_parameter_importance = {}
all_importance_values = []
for ID in particle_IDs_to_evaluate:
    parameter_importance = get_combined_parameter_importance(ID)
    ae_parameter_importance = get_ae_parameter_importance(ID)
    an_parameter_importance = get_an_parameter_importance(ID)
    # logging.info(f"Parameter importance for particle {ID}: \n {parameter_importance}")
    # logging.info(f"Aerobic parameter importance for particle {ID}: \n {ae_parameter_importance}")
    # logging.info(f"Anaerobic parameter importance for particle {ID}: \n {an_parameter_importance}")
    # dump_pickle(parameter_importance, f"../transfer/combined_parameter_importance_{ID}.pkl")
    # dump_pickle(ae_parameter_importance, f"../transfer/ae_parameter_importance_{ID}.pkl")
    # dump_pickle(an_parameter_importance, f"../transfer/an_parameter_importance_{ID}.pkl")


    all_importance_values.extend(parameter_importance.values())
    combined_parameter_importance = {k: combined_parameter_importance.get(k, 0) + v for k, v in parameter_importance.items()}
    #logging.info(f"Parameter importance for particle {ID}: \n {parameter_importance}")
    zero_importance_enzymes = {k for k, v in parameter_importance.items() if v == 0}
    zero_importance_enzymes_list.append(zero_importance_enzymes)
    logging.info(f"Number of zero values in param_importance for particle {ID}: {len(zero_importance_enzymes)}")

combined_zero_importance_enzymes = set.intersection(*zero_importance_enzymes_list)
#logging.info(f"Combined zero importance enzymes: {combined_zero_importance_enzymes}")
logging.info(f"Number of enzymes that have zero importance in all simulations: {len(combined_zero_importance_enzymes)}")

logging.info(f"Combined parameter importance: {combined_parameter_importance}")
combined_normalized_importance = {k: v / len(particle_IDs_to_evaluate) for k, v in combined_parameter_importance.items()}

dump_pickle(combined_normalized_importance, f"{outdir}/combined_normalized_importance.pkl")

# value_counts = Counter(combined_normalized_importance.values())

# # Optional: sort by value
# sorted_items = sorted(value_counts.items())  # list of (value, count)

# # Step 3: Plot
# values, counts = zip(*sorted_items)

# plt.figure(figsize=(10, 5))
# plt.bar(values, counts, width=0.01)  # you can adjust width based on how close the values are
# plt.yscale('log')  # Optional: log scale for better visibility
# plt.xlabel("Parameter importance value")
# plt.ylabel("Number of particles with this value")
# plt.title(f"Distribution of average parameter importance values across {len(particle_IDs_to_evaluate)} particles")
# plt.grid(True)
# plt.tight_layout()
# plt.show()
# plt.savefig(f"../figures/parameter_importance_distribution.png")



# df_flux = load_pickle(f"{outdir}/combined_flux_data_120677.0.pkl")
# logging.info(f"Loaded df_r: \n {df_flux}")
# logging.info(f"df_r shape: {df_flux.shape}")
# logging.info(f"df_r columns: {df_flux.columns[:10]}")
# #filtered_df = df_r[df_r.index.str.match(r"r_\d{4}No1")]

# #param_dict = format_input(params, param_dict=df_flux.to_dict(orient='index'))


# df_flux_importance = (df_flux > 0).sum(axis=0) / len(df_flux)
# df_flux_importance = df_flux_importance.to_frame().T  # Make it a single-row DataFrame
# #df_flux_importance = df_flux_importance.loc[:, df_flux_importance.ne(0.0).all()]


# logging.info(f"df_flux_importance:\n{df_flux_importance}")


# logging.info(f"Loaded ae_model: {type(ae_model)}")
# # logging.info(f"Loaded an_model: {an_model}")

# # logging.info(f"Reactions in ae_model: {len(reactions)}")
# # logging.info(f"Loaded ReFramed ae_model: {ae_model}")
# logging.info(f"Reactions in ae_model: {len(ae_model.reactions)}")


# reactions: set = set()
# ids = set()
# param_importance : Dict[str, float] = {ID: 0.0 for ID in protein_IDs}
# logging.info(f"Parameter importance length: {len(param_importance)}")

# for rxn_id, rxn in ae_model.reactions.items():
#     if not rxn_id.startswith('draw_'): continue
#     genes: set = rxn.get_genes()
#     for met in rxn.stoichiometry:
#         if not met.startswith('prot_'): continue
#             # ingore metabolite: prot_pool
#         if met == 'prot_pool': continue
#         uniprot_id = met.split('_')[1]
#         # if uniprot_id in protein_IDs:
#         logging.info(f"Reaction {rxn_id}: {rxn} is associated with uniprot_id: {uniprot_id}")
#         reactions.add(rxn_id)
#         ids.add(uniprot_id)
#         param_importance[uniprot_id] += df_flux_importance[rxn_id].values[0]
   

# logging.info(f"Enzyme calyzed reactions in ae_model: {len(reactions)}")
# logging.info(f"Enzyme IDs in ae_model: {len(ids)}")

# #logging.info(f"Parameter importance: {param_importance}")
# non_zero_importance = {k: v for k, v in param_importance.items() if v > 0}
# logging.info(f"Number of non-zero values in param_importance: {len(non_zero_importance)}")

# sorted_param_importance = dict(sorted(param_importance.items(), key=lambda item: item[1], reverse=True))

# logging.info(f"Sorted param importance: {list(sorted_param_importance)[:10]}")





# for col in df_r.columns:
#     this column is a reaction flux
#     need to check which enzyme catalyzes this reaction
#     convert the reaction importance to enzyme importance
