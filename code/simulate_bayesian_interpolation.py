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
import pandas as pd
import pickle
import abc_etc as etc
import numpy.typing as npt
from itertools import combinations
from multiprocessing import Process,cpu_count,Manager

# In[2]

# Preliminary definitions and settings
n_comparisons = 10
random_seed = 249
rng = np.random.default_rng(random_seed)
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
        # The timeout is an emergency hatch designed to catch 
        # processes which for some reason are caught in a deadlock
        for p in jobs: p.join(timeout=1000)
        
        distances = [np.inf for _ in range(len(particles))]
        simulated_data = [None for _ in range(len(particles))]

        while not Q.empty():
            index,res = Q.get(timeout=1)
            distances[index] = distance_function(Yobs,res)
            simulated_data[index] = res
        
        # Q may not always contain the result of all jobs we passed to it,
        # this must be handled carefully
        return distances,simulated_data

# In[3]:

def read_posterior_particles(filename: str):
    model: etc.SMCABC = pickle.load(open(file=filename,mode='rb'))
    return get_posterior_particles(model=model)

def get_posterior_particles(model: etc.SMCABC):
    return [particle for particle, distance in
     zip(model.all_particles,model.all_distances) if distance < -0.90]

def extract_posterior_particles(model: abc.SMCABC, r2_threshold = 0.9):
    posterior_idxs = np.nonzero(np.array(model.all_distances) < -r2_threshold)[0]
    return [model.all_particles[idx] for idx in posterior_idxs] 


model_frame: pd.DataFrame = pickle.load(open("../results/permuted_smcabc_res/simulation_skeleton.pkl",'rb'))

model_frame["posterior_particles"] = list(map(read_posterior_particles,
model_frame.outfile))
model_frame.set_index(["origin","status"], inplace=True)

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
    from_posterior_particles: List[abc.candidateType] = model_frame['posterior_particles'][from_model_name]
    to_posterior_particles: List[abc.candidateType] = model_frame['posterior_particles'][to_model_name]
    
    from_particles = rng.choice(from_posterior_particles, size=n_comparisons, replace=True)
    to_particles = rng.choice(to_posterior_particles, size=n_comparisons, replace=True)
    for from_particle, to_particle in zip(from_particles, to_particles):
        itermediate_models = list(map(lambda ratio: create_intermediate_model(from_model=from_particle,to_model=to_particle,ratio=ratio), ratios))
        distances, simulated_data = calculate_distances_parallel(itermediate_models)
        case_result.append(pd.DataFrame({'ratio' : ratios,'distances' : distances, 'simulated_data' : simulated_data}))
    results.append(case_result)

result_frame['results'] = results

pickle.dump(obj=result_frame, file=open('../results/permuted_smcabc_res/bayesian_interpolation_results.pkl','wb'))
