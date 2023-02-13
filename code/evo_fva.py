#!/usr/bin/env python
# coding: utf-8

# This script is supposed to be used with the results from gem_smcevo_run.py
# to perform FVA on the best fitting particles in order to compare them
import logging
from typing import Dict
import pickle
import evo_etc as evo
import numpy as np
import GEMS
import multiprocessing
from functools import partial

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
    model: evo.CrowdingDE = load_pickle(filename)
    return get_posterior_particles(model=model)

def get_posterior_particles(model: evo.CrowdingDE, r2_threshold):
    return [particle for particle, distance in
     zip(model.all_particles,model.all_distances) if distance < -r2_threshold]

fva_functional = partial(GEMS.run_fva_at_two_conditions,reactions=list(signature_reactions.values()))


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("Reading data")
outdir = '../results/crowdingDE'
model_frame = load_pickle(f"{outdir}/simulation_skeleton.pkl")
reduced_frame = model_frame.loc[:, ["scaling_factor","crossover_prob", "simulation","outfile"]].pipe(lambda df: df[(df["scaling_factor"] == 0.5) & (df["crossover_prob"] == 0.99)])

logging.info("Loading data")
posterior_particles = map(read_posterior_particles,reduced_frame.outfile)
reduced_frame["sampled_particles"] = [rng.choice(a=particle_collection,size=min(N_SAMPLES,len(particle_collection)),replace=False) for
    particle_collection in posterior_particles]


dump_pickle(reduced_frame,f"../results/evo_fva_frame.pkl")

logging.info("Running FVA")
with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as p:
    fva_frame = (reduced_frame[["scaling_factor","crossover_prob", "simulation", "sampled_particles"]].
    explode("sampled_particles",ignore_index=True).
    dropna().
    rename(columns={"sampled_particles": "particle"}).
    assign(fva_res = lambda df: list(p.map(fva_functional,df.particle)))
    )
logging.info("Saving results")
pickle.dump(obj=fva_frame, file=open(f"{outdir}/evo_fva.pkl",'wb'))
logging.info("DONE")
