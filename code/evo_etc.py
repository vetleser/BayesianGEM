#!/usr/bin/env python
# coding: utf-8

# #### Add constraint: Topt<Tm

# In[1]:


import logging
from math import sqrt
from typing import Callable, Dict
from inspyred.ec import replacers, variators
import numpy as np
import numpy.typing as npt
from numpy.random.mtrand import seed
import scipy
import scipy.stats as ss
from scipy.stats import poisson
from multiprocessing import Process,cpu_count,Manager
from decimal import Decimal
from random_sampler import RV
import inspyred
import pickle
import inspyred.ec
import inspyred.ec.variators as variators
import inspyred.ec.evaluators as evaluators
import inspyred.ec.selectors as selectors
import os
import GEMS
import random as rand
simResultType = Dict[str, Dict[str, npt.NDArray[np.float64]]]
individualType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]


# In[2]:

random_seed = 12168


class GA:
    def __init__(self,simulator: Callable[[candidateType], simResultType], priors : individualType, min_epsilon: float,
    distance_function: Callable[[distanceArgType, distanceArgType], float],
                 Yobs: distanceArgType, outfile: str,cores: int = cpu_count(),generation_size: int =128, mutation_frequency: float = 1.,
                  selection_proportion: float = 0.5, maxiter: int = 100000):
        '''
        simulator:       a function that takes a dictionary of parameters as input. Ouput {'data':Ysim}
        priors:          a dictionary which use id of parameters as keys and RV class object as values
        min_epsilon:     minimal epsilon
        generation_size: the size of each population
        distance_function: a function that calculate the distance between observed data and simulated data
        Yobs:            observed data
        outfile:         unique id for the experiment. This will be also used to continue a simulation that 
                         is partly done
        cores:           number of treads
        mutation_frequency: The poisson lambda parameter to specify the number of mutations for each differentiation
        selection_proportion: The proportion of the population in each generation used to produce offspring
        maxiter: The maximum number of generations to simulate

        
        !!!Important: distance is to be minimized!!!
        '''
        self.simulator = simulator
        self.priors = priors
        self.distance_function = distance_function
        self.min_epsilon = min_epsilon
        self.Yobs = Yobs
        self.outfile = outfile
        self.population = []  # a list of populations [p1,p2...]
        self.cores = cores    
        self.epsilons = [np.inf]          # min distance in each generation
        self.generation_size = generation_size   # number of particles to be simulated at each generation
        self.mutation_frequency = mutation_frequency
        self.selection_proportion = selection_proportion
        self.maxiter = maxiter

    
        def archiver(random, population, archive, args):
            self.population.append(population)
            max_generation_epsilon = max(p.fitness for p in population)
            self.epsilons.append(max_generation_epsilon)
            logging.info(f"Model epsilon {max_generation_epsilon}")
            pickle.dump(self,self.outfile, 'wb')


        def terminator(population, num_generations, num_evaluations, args) -> bool:
            return self.epsilons[-1] <= self.min_epsilon or num_generations > self.maxiter

        def evaluator(candidate: individualType, args: dict) -> float:
            res = self.simulator(args = {key: value.loc for key, value in candidate.items()})
            distance = self.distance_function(Yobs, res)
            return distance


        def generator(random, args):
            def correct_validity(seed: individualType):
                seed = self.priors
                # check Topt < Tm, if false, resample from prior
                protein_ids = set(entry.split('_')[0] for entry in seed)
                for id in protein_ids:
                    Tm_key = id + "_Tm"
                    Topt_key = id + "_Topt"
                    for _ in range(10):
                        # We try to get things right 10 times before we give up
                        Tm = seed[Tm_key]
                        Topt = seed[Topt_key]
                        if Tm < Topt:
                            seed[Tm_key] = Tm.mutate()
                            seed[Topt_key] = Topt.mutate()
                        else:
                            break
            return correct_validity(self.priors)
        

        @variators.mutator
        def mutate(random: int, candidate: individualType, args: dict):
            def correct_validity(entry):
                # As we only change one parameter at a time, we only need to check
                # the validity of the parameters of one enzyme
                protein_id = entry.split('_')[0]
                Tm_key = protein_id + "_Tm"
                Topt_key = protein_id + "_Topt"
                for _ in range(10):
                    # We try to get things right 10 times before we give up
                    Tm = candidate[Tm_key]
                    Topt = candidate[Topt_key]
                    if Tm < Topt:
                        candidate[Tm_key] = Tm.mutate()
                        candidate[Topt_key] = Topt.mutate()
                    else:
                        break

            n_mutations = poisson.rvs(mu=self.mutation_frequency)
            entries: list[str] = list(candidate.keys())
            for _ in range(n_mutations):
                entry = rand.choice(entries)
                candidate[entry] = candidate[entry].mutate()
                correct_validity(entry)


        @variators.crossover
        def cross(random: int, mom: individualType, dad: individualType, args: dict):
            def merge_distributions(mom: RV, dad: RV):
                if mom.dist_name != dad.dist_name:
                    raise ValueError("Attemt to merge distributions of different types")
                dist_name = mom.dist_name
                # For normally distributed variables, the take the average of the variances
                scale = sqrt((mom.scale**2 + dad.scale**2) / 2) if dist_name == 'normal' else (mom.scale + dad.scale) / 2
                loc = (mom.loc + dad.loc) / 2
                return RV(dist_name, loc, scale)
                

            common_keys = set(mom.keys())
            common_keys.update(dad.keys())
            return {key: merge_distributions(mom[key], dad[key]) for key in common_keys}


        self.archiver = archiver
        self.terminator = terminator
        self.evaluator = evaluator
        self.generator = generator
        self.mutator = mutate
        self.crossover = cross


    
    def run_simulation(self):
        algorithm = inspyred.ec.EvolutionaryComputation(random_seed)
        algorithm.terminator = self.terminator
        algorithm.archiver = self.archiver
        algorithm.variator = [self.crossover, self.mutator]
        algorithm.selector = selectors.tournament_selection
        algorithm.replacer = replacers.truncation_replacement


        kwargsdict = {'mp_evaluator' : self.evaluator, 'mp_nprocs': self.cores, 'num_selected': int(self.generation_size * self.selection_proportion)}
        algorithm.evolve(generator=self.generator, evaluator=evaluators.parallel_evaluation_mp, maximize=False, pop_size=self.generation_size, args = kwargsdict)
