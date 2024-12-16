#!/usr/bin/env python
# coding: utf-8

# This script is intended to compare the Bayesian calculation method and evolutionary algorithm on a bimodal objective function

import random_sampler
import numpy as np
import abc_etc as abc
import evo_etc as evo
import os
import math
import logging
import copy

# logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
random_seed = 5354 #Changing the seed gives new evolutionary population. Bayesian population is new even for same seed
maxiter = 200
Yobs = None
min_epsilon = -1
population_size = 128
outdir = "./results/toy_example" 
if not os.path.exists(outdir):
    os.makedirs(outdir)


def simulator(candidate):
    # There is really nothing to similate in this toy example, so we just return the answer
    logging.debug(candidate)
    return candidate



# Due to the API for which the function is used, it must have two arguments even though we only care about the first one
def fitness_function(dummy,candidate):
    denomenator = 1+((candidate["x"]-1)**2+(candidate["y"]-1)**2)*((candidate["x"]+1)**2+(candidate["y"]+1)**2)
    R2 = 1/denomenator
    logging.debug(f"R2={R2}")
    return -R2


def shifted_fitness_function(dummy, candidate):
    denomenator = 1+math.sqrt((candidate["x"]-1)**2+(candidate["y"]-1)**2)*(math.sqrt((candidate["x"]+1)**2+(candidate["y"]+1)**2)+.5)
    R2 = 1/denomenator
    logging.debug(f"R2={R2}")
    return -R2


rng = np.random.default_rng(random_seed)

n_iterations = 4

priors = {var: random_sampler.RV(dist_name='normal',loc=0,scale=0.2,rng=rng) for var in ("x","y")}

def distribution_is_bimodal(model: evo.CrowdingDE, tol = 10e-3):
    final_population = model.population[-1]
    final_particles = [model.all_particles[i] for i in final_population]
    upper_optimum_reached = False
    lower_optimum_reached = False
    for particle in final_particles:
        if abs(particle["x"] - 1) < tol and abs(particle["y"] - 1) < tol:
            upper_optimum_reached = True
        if abs(particle["x"] + 1) < tol and abs(particle["y"] + 1) < tol:
            lower_optimum_reached = True
    print("Upper optimum reached" if upper_optimum_reached else "Upper optimum NOT reached")
    print("Lower optimum reached" if lower_optimum_reached else "Lower optimum NOT reached")
    return upper_optimum_reached and lower_optimum_reached


# for i in range(n_iterations):
#     bayesian_model = abc.SMCABC(simulator=simulator,
#                                     priors=copy.deepcopy(priors),
#                                     min_epsilon=-1,
#                                     population_size=32,
#                                     distance_function=fitness_function,
#                                     Yobs=Yobs,
#                                     outfile=f"{outdir}/bayesian_{i}.pkl",
#                                     generation_size=32,
#                                     cores=1,
#                                     maxiter=maxiter)
#     bayesian_model.run_simulation()

for j in range(4):
    random_seed += 1
    rng = np.random.default_rng(random_seed)
    for i in range(n_iterations):
        crowdingDE_model = evo.CrowdingDE(simulator=simulator,
                                    priors=copy.deepcopy(priors),
                                    min_epsilon=-1,
                                    generation_size=32,
                                    distance_function=fitness_function,
                                    Yobs=Yobs,
                                    outfile=f"{outdir}/crowdingDE_{j}_{i}.pkl",
                                    maxiter=maxiter,
                                    rng=rng,
                                    cores=1,
                                    crossover_prob=.5,
                                    n_children=16,
                                    scaling_factor=0.5,
                                    save_intermediate=False
                                    )
        crowdingDE_model.run_simulation()



# for i in range(n_iterations):
#     crowdingDE_model = evo.CrowdingDE(simulator=simulator,
#                                 priors=copy.deepcopy(priors),
#                                 min_epsilon=-1,
#                                 generation_size=32,
#                                 distance_function=fitness_function,
#                                 Yobs=Yobs,
#                                 outfile=f"{outdir}/crowdingDE_{i}.pkl",
#                                 maxiter=maxiter,
#                                 rng=rng,
#                                 cores=1,
#                                 crossover_prob=.5,
#                                 n_children=16,
#                                 scaling_factor=0.5,
#                                 save_intermediate=False
#                                 )
#     crowdingDE_model.run_simulation()
