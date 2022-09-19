#!/usr/bin/env python
# coding: utf-8

# This script is supposed to be used with the results from gem_smcevo_at_three_conditions_permuted.py
# to perform FVA on the posterior distributions in order to compare them
from inspect import signature
import logging
from typing import Dict
import pandas as pd
import pickle
import evo_etc as evo
import numpy as np
import GEMS
from functools import partial
import multiprocessing

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

rng = np.random.default_rng(seed=16501)

N_SAMPLES = 20
N_CORES = multiprocessing.cpu_count()
cpu_pool = multiprocessing.Pool(processes=N_CORES)
signature_reactions = {'PDH': 'r_0961No1', 'FBA': 'r_0450No1', 'FCO': 'r_0438No1', 'PSP': 'r_0917No1', 'SHK': 'r_0997No1', 'GRW': 'r_2111'}

candidateType = Dict[str, float]

def read_posterior_particles(filename: str):
    model: evo.GA = pickle.load(open(file=filename,mode='rb'))
    return get_posterior_particles(model=model)

def get_posterior_particles(model: evo.GA):
    return [particle for particle, distance in
     zip(model.all_particles,model.all_distances) if distance < -0.90]

fva_functional = partial(GEMS.run_fva_at_three_conditions,reactions=list(signature_reactions.values()))


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("Reading data")
posterior_particles = read_posterior_particles('../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl')


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("Reading data")
full_simulation_skeleton: pd.DataFrame = pickle.load(open("../results/permuted_smcevo_res/simulation_skeleton.pkl",'rb'))
reduced_simulation_skeleton = full_simulation_skeleton.loc[:, ["prior_name","simulation","outfile"]]
reduced_simulation_skeleton["posterior_particles"] = list(map(read_posterior_particles,  
reduced_simulation_skeleton.outfile))
reduced_simulation_skeleton["sampled_particles"] = [rng.choice(a=particle_collection,size=N_SAMPLES,replace=False) for
 particle_collection in reduced_simulation_skeleton["posterior_particles"]]
logging.info("Running FVA")
fva_frame = (reduced_simulation_skeleton[["prior_name","simulation", "sampled_particles"]].
explode("sampled_particles",ignore_index=True).
rename(columns={"sampled_particles": "particle"}).
assign(fva_res = lambda df: cpu_pool.map(func=fva_functional,iterable=df.particle))
)
logging.info("Saving results")
pickle.dump(obj=fva_frame, file=open("../results/permuted_smcevo_res/evo_fva.pkl",'wb'))
logging.info("DONE")
