#!/usr/bin/env python
# coding: utf-8

# This script is intended to compare the Bayesian calculation method and evolutionary algorithm on a bimodal objective function

import random_sampler
import numpy as np
import abc_etc as abc
import evo_etc as evo
import sa_etc as sa
import os
import math
import logging
import copy
import time


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
random_seed = 5354 #Changing the seed gives new evolutionary (and simulated annealing) population . Bayesian population is new even for same seed
maxiter = 200
Yobs = None
min_epsilon = -1
population_size = 128
outdir = "../results/sa" 
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

# start = time.time()
# for j in range(4):
#     random_seed += 1
#     rng = np.random.default_rng(random_seed)
#     for i in range(n_iterations):
#         crowdingDE_model = evo.CrowdingDE(simulator=simulator,
#                                     priors=copy.deepcopy(priors),
#                                     min_epsilon=-1,
#                                     generation_size=32,
#                                     distance_function=fitness_function,
#                                     Yobs=Yobs,
#                                     outfile=f"{outdir}/crowdingDE_{j}_{i}.pkl",
#                                     maxiter=maxiter,
#                                     rng=rng,
#                                     cores=1,
#                                     crossover_prob=.5,
#                                     n_children=16,
#                                     scaling_factor=0.5,
#                                     save_intermediate=False
#                                     )
#         crowdingDE_model.run_simulation()

# end = time.time()
# logging.info(f"Time for evo: {end-start} seconds")
# start = time.time()
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
# end = time.time()
# logging.info(f"Time for evo: {end-start} seconds")

min_layer_list = [1, 1, 3, 3]
max_layer_list = [1, 3, 5, 10]

logging.info("Attempting Simulated Annealing")
start1 = time.time()
for j in range(4):
    rng = np.random.default_rng(random_seed)  # fresh RNG per sim
    for i in range(n_iterations):
        start = time.time()
        logging.info(f"Simulated Annealing plot {j}, simulation {i} started")
        simanneal_model = sa.SimulatedAnnealing(
                            simulator=simulator,
                            priors=copy.deepcopy(priors),
                            min_epsilon=-1,
                            distance_function=fitness_function,
                            Yobs=Yobs,
                            maxiter=maxiter,
                            generation_size = 32,
                            outfile=f"{outdir}/simanneal_minmax_{j}_{i}.pkl",
                            initial_temp=100,
                            cooling_rate=0.95,
                            final_temp=1,
                            rng=rng,
                            cores=1,
                            version=1,
                            min_layers=min_layer_list[j],
                            max_layers=max_layer_list[j]
                            )
        simanneal_model.run_simulation()
        end = time.time()
        logging.info(f"Simulated Annealing {j} finished in {end-start} seconds")


end1 = time.time()

logging.info(f"Simulated annealing: {end1-start1} seconds")

# logging.info("Attempting Simulated Annealing")
# start = time.time()
# for i in range(n_iterations):
#     rng = np.random.default_rng(random_seed)  # fresh RNG per sim
#     simanneal_model = sa.SimulatedAnnealing(
#                         simulator=simulator,
#                         priors=copy.deepcopy(priors),
#                         min_epsilon=-1,
#                         distance_function=fitness_function,
#                         Yobs=Yobs,
#                         maxiter=maxiter,
#                         generation_size = 32,
#                         outfile=f"{outdir}/simanneal_{i}.pkl",
#                         initial_temp=100,
#                         cooling_rate=0.95,
#                         final_temp=1,
#                         rng=rng,
#                         cores=1
#                         )
#     simanneal_model.run_simulation_2()

logging.info("DONE")

