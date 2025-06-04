#!/usr/bin/env python
# coding: utf-8

# This script is intended to compare the Bayesian calculation method and evolutionary algorithm on a bimodal objective function

from functools import partial
import random_sampler
import numpy as np
import abc_etc as abc
import evo_etc as evo
import sa_etc as sa
import sa_etc_new as sa_new
import os
import math
import logging
import copy
import time
import pandas as pd
import dill


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return dill.load(open(file=filename,mode='rb'))
        
def dump_pickle(obj,filename):
    return dill.dump(obj=obj,file=open(file=filename, mode='wb'))

random_seed = 5354 #Changing the seed gives new evolutionary (and simulated annealing) population . Bayesian population is new even for same seed
maxiter = 200
Yobs = None
min_epsilon = -1
population_size = 32
outdir = "../results/toy_example_new" 
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

def fitness_function2(dummy,candidate):
    denomenator = 1+((candidate["x"]-1)**2+(candidate["y"]-1)**2)*((candidate["x"]+1)**2+(candidate["y"]+1)**2)*((candidate["x"]+1)**2+(candidate["y"]-1)**2)*((candidate["x"]-1)**2+(candidate["y"]+1)**2)
    R2 = 1/denomenator
    logging.debug(f"R2={R2}")
    return -R2


def fitness_function3(dummy, candidate):
    minima_list = [(-1, -1), (1, 1), (-1, 1), (1, -1)]
    sigma = 0.6  # Standard deviation of each Gaussian
    R2 = 0

    for min_x, min_y in minima_list:
        dx = candidate["x"] - min_x
        dy = candidate["y"] - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma**2)
        R2 += math.exp(exponent)

    return -R2/1.0077467856174704

def fitness_function4(dummy, candidate):
    minima_list = [(-1, 1), (1, 1), (-1, -1), (1, -1)]
    sigma = 0.6  # Standard deviation of each Gaussian
    R2 = 0
    scaling_list = [1, 0.9, 0.8, 0.7]

    for i, (min_x, min_y) in enumerate(minima_list):
        dx = candidate["x"] - min_x
        dy = candidate["y"] - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma**2)
        R2 += math.exp(exponent) * scaling_list[i]

    return -R2/1.006582525974071

def general_fitness_function(dummy, candidate, minima_list, scaling_list, sigma=0.6):
    def compute_R2(x, y):
        R2 = 0
        for i, (min_x, min_y) in enumerate(minima_list):
            dx = x - min_x
            dy = y - min_y
            exponent = -(dx**2 + dy**2) / (2 * sigma**2)
            R2 += math.exp(exponent) * scaling_list[i]
        return R2

    max_R2 = max(compute_R2(x, y) for (x, y) in minima_list)
    value = compute_R2(candidate["x"], candidate["y"])
    return -value / max_R2

def rastrigin_function(dummy, candidate):
    A = 10
    n = 2
    x, y = candidate["x"], candidate["y"]
    # Calculate the Rastrigin function for x and y
    sum1 = (x**2 - A * np.cos(2 * np.pi * x))
    sum2 = (y**2 - A * np.cos(2 * np.pi * y))
    
    # Sum them together with the constant A * n
    return -(A * n + sum1 + sum2)/80.70658039


test1 = fitness_function4(None, {"x": 0, "y": 0})
test2 = fitness_function4(None, {"x": 1, "y": 1})
test3 = fitness_function4(None, {"x": -1, "y": -1})
test4 = fitness_function4(None, {"x": 1, "y": -1})
test5 = fitness_function4(None, {"x": -1, "y": 1})
logging.info(f"Test fitness function 4, (0, 0): {test1}")
logging.info(f"Test fitness function 4, (1, 1): {test2}")
logging.info(f"Test fitness function 4, (-1, -1): {test3}")
logging.info(f"Test fitness function 4, (1, -1): {test4}")
logging.info(f"Test fitness function 4, (-1, 1): {test5}")

def fitness_function5(dummy, candidate):
    minima_list = [(-1, 1), (1, 1), (-1, -1), (-3, 3)]
    sigma = 0.6  # Standard deviation of each Gaussian
    R2 = 0
    scaling_list = [1, 0, 0, 1]

    for i, (min_x, min_y) in enumerate(minima_list):
        dx = candidate["x"] - min_x
        dy = candidate["y"] - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma**2)
        R2 += math.exp(exponent) * scaling_list[i]

    return -R2/1.0054182663306717



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

def testing_fitness_function(fitness_function, minima_list):
    for candidate in minima_list:
        fitness = fitness_function(None, {"x": candidate[0], "y": candidate[1]})
        logging.warning(f"Candidate {candidate}: fitness = {fitness:.4f}")


def perform_sa(generation_size, 
               final_temp_list, 
               filename, 
               step_size, 
               minima_list, 
               scaling_list, 
               n_plots=4, 
               n_iterations=4,
               fitness_function = None,
               min_epsilon: float = -1,):
    
    if fitness_function is None:
        fitness_function = partial(general_fitness_function, minima_list=minima_list, scaling_list=scaling_list)
    #testing_fitness_function(fitness_function, minima_list)
    cooling_rate_list = [(final_temp / initial_temp) ** (1 / maxiter) for final_temp in final_temp_list]
    for j in range(n_plots):
        rng = np.random.default_rng(random_seed)  # fresh RNG per sim
        for i in range(n_iterations):
            start = time.time()
            logging.warning(f"Simulated Annealing plot {j} of {n_plots}, simulation {i} of {n_iterations} started")
            simanneal_model = sa.SimulatedAnnealing(
                                simulator=simulator,
                                priors=copy.deepcopy(priors),
                                min_epsilon=min_epsilon,
                                distance_function=fitness_function,
                                Yobs=Yobs,
                                maxiter=maxiter,
                                generation_size = generation_size,
                                outfile=f"{outdir}/{filename}_{j}_{i}.pkl",
                                initial_temp=initial_temp,
                                cooling_rate=cooling_rate_list[j],
                                final_temp=final_temp_list[j],
                                rng=rng,
                                cores=1,
                                version=2,
                                step_size=step_size,
                                )
            simanneal_model.run_simulation()
            end = time.time()
            #logging.info(f"Time for SA {j}: {end-start} seconds")

def perform_sa_new(generation_size, 
               final_temp_list, 
               filename, 
               step_size, 
               minima_list, 
               scaling_list, 
               n_plots=4, 
               n_iterations=4,
               fitness_function = None,
               min_epsilon: float = -1,):
    
    if fitness_function is None:
        fitness_function = partial(general_fitness_function, minima_list=minima_list, scaling_list=scaling_list)
    testing_fitness_function(fitness_function, minima_list)
    cooling_rate_list = [(final_temp / initial_temp) ** (1 / maxiter) for final_temp in final_temp_list]
    for j in range(n_plots):
        rng = np.random.default_rng(random_seed)  # fresh RNG per sim
        for i in range(n_iterations):
            start = time.time()
            logging.warning(f"Simulated Annealing plot {j} of {n_plots}, simulation {i} of {n_iterations} started")
            simanneal_model = sa_new.SimulatedAnnealing(
                                simulator=simulator,
                                priors=copy.deepcopy(priors),
                                min_epsilon=min_epsilon,
                                distance_function=fitness_function,
                                Yobs=Yobs,
                                maxiter=maxiter,
                                generation_size = generation_size,
                                outfile=f"{outdir}/{filename}_new_{j}_{i}.pkl",
                                initial_temp=initial_temp,
                                cooling_rate=cooling_rate_list[j],
                                final_temp=final_temp_list[j],
                                rng=rng,
                                cores=1,
                                version=2,
                                step_size=step_size,
                                )
            simanneal_model.run_simulation()
            end = time.time()
            #logging.info(f"Time for SA {j}: {end-start} seconds")
        


figures = [0, 1, 2, 3, 4, 5]
final_temp_list1 = [1.0, 0.1, 0.01, 0.001]
final_temp_list2 = [0.01, 0.001, 0.0001, 0.00001]
minima_list1 = [(-1, 1), (1, 1), (-1, -1), (1, -1)]
minima_list2 = [(1, 1), (2, 2)]
scaling_list1 = [1, 1, 1, 1]
scaling_list2 = [1, 0.8, 0.6, 0.4]
scaling_list3 = [1, 0.8, 0.6, 0]
scaling_list4 = [1, 1]
step_size1 = 0.1
step_size2 = 1
x_lim1 = (-2.5, 2.5)
y_lim1 = (-2.5, 2.5)
x_lim2 = (-1, 5)
y_lim2 = (-1, 5)
population_size = 32


toy_example_df = pd.DataFrame({"Figure": figures})

toy_example_df["Final_temp_list"] = [final_temp_list1] * 5 + [final_temp_list2] * 1
toy_example_df["Minima_list"] = [minima_list1] * 3 + [minima_list2] * 3
toy_example_df["Scaling_list"] = [scaling_list1] + [scaling_list2] + [scaling_list3] + [scaling_list4] * 3
toy_example_df["Step_size"] = [step_size1] * 4 + [step_size2] * 2
toy_example_df["X_lim"] = [x_lim1] * 3 + [x_lim2] * 3
toy_example_df["Y_lim"] = [y_lim1] * 3 + [y_lim2] * 3
toy_example_df["Generation_size"] = [population_size] * 6

with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
    logging.warning(f"Toy example dataframe:\n{toy_example_df}")

#dump_pickle(toy_example_df, f"{outdir}/toy_example_df.pkl")

toy_example_df2 = pd.DataFrame({"Figure": [6, 7, 8]})
toy_example_df2["Final_temp_list"] = [final_temp_list1] * 2 + [final_temp_list2] * 1
toy_example_df2["Minima_list"] = [[(1, 1), (3, 3)]] * 3
toy_example_df2["Scaling_list"] = [scaling_list4] * 3
toy_example_df2["Step_size"] = [step_size1] * 1 + [step_size2] * 2
toy_example_df2["X_lim"] = [x_lim2] * 3
toy_example_df2["Y_lim"] = [y_lim2] * 3
toy_example_df2["Generation_size"] = [population_size] * 3
with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
    logging.warning(f"Toy example dataframe 2:\n{toy_example_df2}")

dump_pickle(toy_example_df2, f"{outdir}/toy_example_df2.pkl")

initial_temp = 100
#final_temp_list = [1.0, 0.1, 0.01, 0.001]
# final_temp_list = [0.0001, 0.00001, 0.000001, 0.0000001]
#cooling_rate_list = [(final_temp / initial_temp) ** (1 / maxiter) for final_temp in final_temp_list]

n_plots = 4
n_simulations = 4

perform_sa(
    generation_size=32,
    final_temp_list=[1.0, 0.1, 0.01, 0.001],
    filename="simanneal_rastr",
    step_size=1,
    minima_list=[4.5229936666666666666666666, 4.5229936666666666666666666],
    scaling_list=[],
    n_plots=n_plots,
    n_iterations=n_simulations,
    fitness_function=rastrigin_function,
    min_epsilon=-1

)

# for figure in toy_example_df["Figure"]:
#     row = toy_example_df.loc[toy_example_df["Figure"] == figure].iloc[0]  # safely get the matching row
#     logging.warning(f"Performing simulated annealing for figure {figure}")
#     perform_sa_new(
#         #fitness_function= rastrigin_function,
#         generation_size=row["Generation_size"],
#         final_temp_list=row["Final_temp_list"],
#         filename=f"simanneal_fig{figure}",
#         step_size=row["Step_size"],
#         minima_list=row["Minima_list"],
#         scaling_list=row["Scaling_list"],
#         n_plots=n_plots,
#         n_iterations=n_simulations
#     )
    
    

# for figure in toy_example_df2["Figure"]:
#     row = toy_example_df2.loc[toy_example_df2["Figure"] == figure].iloc[0]  # safely get the matching row
#     logging.warning(f"Performing simulated annealing for figure {figure}")
#     perform_sa(
#         generation_size=row["Generation_size"],
#         final_temp_list=row["Final_temp_list"],
#         filename=f"simanneal_fig{figure}",
#         step_size=row["Step_size"],
#         minima_list=row["Minima_list"],
#         scaling_list=row["Scaling_list"],
#         n_plots=n_plots,
#         n_iterations=n_simulations
#     )


logging.warning("DONE")



































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
# start1 = time.time()
# for i in range(n_iterations):
#     start = time.time()
#     crowdingDE_model = evo.CrowdingDE(simulator=simulator,
#                                 priors=copy.deepcopy(priors),
#                                 min_epsilon=-1,
#                                 generation_size=32,
#                                 distance_function=fitness_function3,
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
#     end = time.time()
#     logging.info(f"Time for evo {i}: {end-start} seconds")

# end1 = time.time()
# logging.info(f"Time for all evo: {end1-start1} seconds")


# logging.info("Attempting Simulated Annealing")
# start1 = time.time()
# for j in range(n_plots):
#     rng = np.random.default_rng(random_seed)  # fresh RNG per sim
#     for i in range(n_iterations):
#         start = time.time()
#         logging.info(f"Simulated Annealing plot {j}, simulation {i} started")
#         simanneal_model = sa.SimulatedAnnealing(
#                             simulator=simulator,
#                             priors=copy.deepcopy(priors),
#                             min_epsilon=-1,
#                             distance_function=fitness_function3,
#                             Yobs=Yobs,
#                             maxiter=maxiter,
#                             generation_size = 100,
#                             outfile=f"{outdir}/simanneal_test_{j}_{i}.pkl",
#                             initial_temp=initial_temp,
#                             cooling_rate=cooling_rate_list[j],
#                             final_temp=final_temp_list[j],
#                             rng=rng,
#                             cores=1,
#                             version=2
#                             )
#         simanneal_model.run_simulation()
#         end = time.time()
#         logging.info(f"Time for SA {j}: {end-start} seconds")


# end1 = time.time()


# logging.info(f"Simulated annealing: {end1-start1} seconds")

# logging.info("Attempting Simulated Annealing, fitness function 4")
# start1 = time.time()
# for j in range(4):
#     rng = np.random.default_rng(random_seed)  # fresh RNG per sim
#     for i in range(n_iterations):
#         start = time.time()
#         logging.info(f"Simulated Annealing plot {j}, simulation {i} started")
#         simanneal_model = sa.SimulatedAnnealing(
#                             simulator=simulator,
#                             priors=copy.deepcopy(priors),
#                             min_epsilon=-1,
#                             distance_function=fitness_function5,
#                             Yobs=Yobs,
#                             maxiter=maxiter,
#                             generation_size = 32,
#                             outfile=f"{outdir}/simanneal_step1_{j}_{i}.pkl",
#                             initial_temp=initial_temp,
#                             cooling_rate=cooling_rate_list[j],
#                             final_temp=final_temp_list[j],
#                             rng=rng,
#                             cores=1,
#                             version=2
#                             )
#         simanneal_model.run_simulation()
#         end = time.time()
#         logging.info(f"Time for SA {j}: {end-start} seconds")


# end1 = time.time()

# logging.info(f"Simulated annealing: {end1-start1} seconds")

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

