#!/usr/bin/env python
# coding: utf-8

# This script is supposed to be used with the results from gem_smcabc_at_three_conditions_permuted.py
# to perform FVA on the posterior distributions in order to compare them
import logging
from typing import Dict
import pandas as pd
import pickle
import abc_etc as etc
import numpy as np
import GEMS
import multiprocessing

rng = np.random.default_rng(seed=46465)

N_SAMPLES = 20
N_CORES = multiprocessing.cpu_count()
cpu_pool = multiprocessing.Pool(processes=N_CORES)

candidateType = Dict[str, float]

def read_posterior_particles(filename: str):
    model: etc.SMCABC = pickle.load(open(file=filename,mode='wb'))
    return get_posterior_particles(model=model)

def get_posterior_particles(model: etc.SMCABC):
    return [particle for particle, distance in
     zip(model.all_particles,model.all_distances) if distance < -0.90]

def run_fva_on_particles(particle_list):
    pass


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("Reading data")
full_simulation_skeleton: pd.DataFrame = pickle.load(open("../results/permuted_smcabc_res/simulation_skeleton.pkl",'rb'))
reduced_simulation_skeleton = full_simulation_skeleton.loc[full_simulation_skeleton.status == "original", ["origin","outfile"]]
reduced_simulation_skeleton["posterior_particles"] = map(read_posterior_particles,  
reduced_simulation_skeleton.outfiles)
reduced_simulation_skeleton["sampled_particles"] = [rng.choice(a=particle_collection,size=N_SAMPLES,replace=True) for
 particle_collection in reduced_simulation_skeleton["posterior_particles"]]
logging.info("Running FVA")
fva_frame = (reduced_simulation_skeleton[["origin","sampled_particles"]].
explode("sampled_particles",ignore_index=True).
rename({"sampled_particle": "particle"}).
assign(fva_res = lambda df: cpu_pool.map(func=GEMS.run_fva_at_three_conditions,iterable=df.particle))
)
logging.info("Saving results")
pickle.dump(obj=fva_frame)
logging.info("DONE")
