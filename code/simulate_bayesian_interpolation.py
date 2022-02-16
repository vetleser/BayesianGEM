#!/usr/bin/env python
# coding: utf-8


# This script is similar to similate_interpolation.py, but uses Bayesian comparsion of posterior
# parameter sets instead of constructing a concensus model for each distribution

# In[1]
from typing import List
from sympy import Float
import abc_etc as abc
import numpy as np
import GEMS
import os
import pandas as pd
import pickle
import permute_parameters
import logging
import numpy.typing as npt
from itertools import combinations
from multiprocessing import Process,cpu_count,Manager

# In[2]

# Preliminary definitions and settings
n_comparisons = 10
simulator = GEMS.simulate_at_three_conditions_2
distance_function = GEMS.distance_2
Yobs_batch = GEMS.aerobic_exp_data()
Yobs_chemo = GEMS.chemostat_exp_data()
#Yobs_batch_an = GEMS.anaerobic_exp_data()
dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}
Yobs = {'rae':Yobs_batch['data'],
        'chemostat':Yobs_chemo['data'],
        'ran':Yobs_batch_an['data']}


def calculate_distances_parallel(particles):
        def simulate_one(particle,index,Q):
            '''
            particle:  parameters 
            Q:      a multiprocessing.Queue object
            index:  the index in particles list
            '''
            res = simulator(particle)
            # ysim = {simulated}

            Q.put((index,res))

        Q = Manager().Queue()
        jobs = [Process(target=simulate_one,args=(particle,index,Q)) 
                               for index,particle in enumerate(particles)]
        
        for p in jobs: p.start()
        for p in jobs: p.join()
        
        distances = [None for _ in range(len(particles))]
        simulated_data = [None for _ in range(len(particles))]

        while not Q.empty():
            index,res = Q.get(timeout=1)
            distances[index] = distance_function(Yobs,res)
            simulated_data[index] = res
        
        # Q may not always contain the result of all jobs we passed to it,
        # this must be handled carefully
        missing_indicies = [i for i, val in enumerate(distances) if val is None]
        # Great care must be taken, deleting indicies must take place in
        # decreasing order in order to not shift the indicies
        missing_indicies.sort(reverse=True)
        # We solve the problem by removing elements corresponding to missing
        # values
        for i in missing_indicies:
            del simulated_data[i]
            del distances[i]
        return distances,simulated_data

# In[3]:

def extract_posterior_particles(model: abc.SMCABC, r2_threshold = 0.9):
    posterior_idxs = np.nonzero(np.array(model.all_distances) < -r2_threshold)[0]
    return [model.all_particles[idx] for idx in posterior_idxs]


n_permutations = 3
infiles = ['../results/smcabc_gem_three_conditions_save_all_particles.pkl'] + [f'../results/smcabc_gem_three_conditions_permuted_{i}_save_all_particles.pkl' for i in [0, 1]]
model_frame = pd.DataFrame({'origin' : ['unpermuted', 'permuted_0', 'permuted_1'], 'file_path': infiles})
model_frame.set_index('origin', inplace=True)
modeling_results = [pickle.load(open(infile,'rb')) for infile in model_frame['file_path']]
model_frame['modeling_results'] = modeling_results
model_frame['posterior_particles'] = [extract_posterior_particles(model=model) for model in modeling_results] 

def create_intermediate_model(from_model: abc.candidateType, to_model: abc.candidateType, ratio: Float):
    # A ratio of 0 will yield the from_model, whereas a ratio of 1 will yield the to_model
    intermediate_model: abc.candidateType = dict()
    parameter_names = from_model.keys()
    for parameter in parameter_names:
        intermediate_value: float = (1 - ratio)*from_model[parameter] + ratio * to_model[parameter]
        intermediate_model[parameter] = intermediate_value
    return intermediate_model

ratios = np.linspace(0,1,11)
result_frame = pd.DataFrame()
model_combinations = list(combinations(model_frame.index,2))
result_frame['from'] = [combination[0] for combination in model_combinations]
result_frame['to'] = [combination[1] for combination in model_combinations]

# In[]
results: List[List[pd.DataFrame]] = []

for from_model_name, to_model_name in zip(result_frame['from'], result_frame['to']):
    case_result = []
    from_posterior_particles: pd.Series[abc.candidateType] = model_frame['posterior_particles'][from_model_name]
    to_posterior_particles: pd.Series[abc.candidateType] = model_frame['posterior_particles'][to_model_name]
    from_particles = from_posterior_particles.sample(n=n_comparisons, replace=True)
    to_particles = to_posterior_particles.sample(n=n_comparisons, replace=True)
    for from_particle, to_particle in zip(from_particle, to_particle):
        itermediate_models = list(map(lambda ratio: create_intermediate_model(from_model=from_particle,to_model=to_particle,ratio=ratio), ratios))
        distances, simulated_data = calculate_distances_parallel(itermediate_models)
        case_result.append(pd.DataFrame({'ratio' : ratios,'distances' : distances, 'simulated_data' : simulated_data}))
    results.append(case_result)

result_frame['results'] = results

pickle.dump(obj=result_frame, file=open('../results/interpolation_results.pkl','wb'))
