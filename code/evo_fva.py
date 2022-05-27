#!/usr/bin/env python
# coding: utf-8

# This script is supposed to be used with the results from gem_smcevo_at_three_conditions_permuted.py
# to perform FVA on the posterior distributions in order to compare them
import logging
from typing import Dict
import pandas as pd
import pickle
import abc_etc as etc
import numpy as np
import GEMS
import multiprocessing

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

rng = np.random.default_rng(seed=16501)

N_SAMPLES = 20
N_CORES = multiprocessing.cpu_count()
cpu_pool = multiprocessing.Pool(processes=N_CORES)

candidateType = Dict[str, float]

def read_posterior_particles(filename: str):
    model: etc.SMCABC = pickle.load(open(file=filename,mode='rb'))
    return get_posterior_particles(model=model)

def get_posterior_particles(model: etc.SMCABC):
    return [particle for particle, distance in
     zip(model.all_particles,model.all_distances) if distance < -0.90]


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("Reading data")
posterior_particles = read_posterior_particles('../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl')
sampled_particles = rng.choice(a=posterior_particles,size=N_SAMPLES,replace=False)
fva_frame = pd.DataFrame({"particle": sampled_particles})
logging.info("Running FVA")
fva_res = cpu_pool.map(func=GEMS.run_fva_at_three_conditions,iterable=sampled_particles)
fva_frame["fva_res"] = list(fva_res)
logging.info("Saving results")
pickle.dump(obj=fva_frame, file=open("../results/evo_fva.pkl",'wb'))
logging.info("DONE")