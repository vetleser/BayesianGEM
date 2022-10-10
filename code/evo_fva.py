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

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

rng = np.random.default_rng(seed=16501)

N_SAMPLES = 20
signature_reactions = {'PDH': 'r_0961No1', 'FBA': 'r_0450No1', 'FCO': 'r_0438No1', 'PSP': 'r_0917No1', 'SHK': 'r_0997No1', 'GRW': 'r_2111'}

candidateType = Dict[str, float]

def read_posterior_particles(filename: str):
    model: evo.GA = load_pickle(filename)
    return get_posterior_particles(model=model)

def get_posterior_particles(model: evo.GA):
    Yobs_batch = GEMS.aerobic_exp_data()
    dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
    sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
    Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}
    Yobs = {'rae':Yobs_batch['data'],
            'ran':Yobs_batch_an['data']}
    reduced_distances = [GEMS.distance_2(Yobs,res) for res in model.all_simulated_data]
    return [particle for particle, distance in
     zip(model.all_particles,reduced_distances) if distance < -0.90]

fva_functional = partial(GEMS.run_fva_at_three_conditions,reactions=list(signature_reactions.values()))


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("Reading data")
outdirs = {'tournament': '../results/evo_tournament', 'truncation': '../results/evo_truncation'}
model_frames = {method: load_pickle(f"{outdir}/simulation_skeleton.pkl") for method, outdir in outdirs.items()}
reduced_frames = {method: frame.loc[:, ["locality" if method == "tournament" else "num_elites", "simulation","outfile"]] for method, frame in model_frames.items()}

logging.info("Loading data")
for frame in reduced_frames.values():
    frame["posterior_particles"] = list(map(read_posterior_particles,frame.outfile))

for frame in reduced_frames.values():
    frame["sampled_particles"] = [rng.choice(a=particle_collection,size=min(N_SAMPLES,len(particle_collection)),replace=False) for
 particle_collection in frame["posterior_particles"]]
logging.info("Running FVA")
fva_frames = {}
for method, frame in reduced_frames.items():
    fva_frames[method] = (frame[["prior_name","simulation", "sampled_particles"]].
    explode("sampled_particles",ignore_index=True).
    rename(columns={"sampled_particles": "particle"}).
    assign(fva_res = lambda df: list(map(fva_functional,df.particle)))
    )
logging.info("Saving results")
pickle.dump(obj=fva_frames, file=open("../results/evo_fva.pkl",'wb'))
logging.info("DONE")
