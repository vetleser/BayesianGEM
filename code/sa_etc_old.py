#!/usr/bin/env python
# coding: utf-8

# #### Add constraint: Topt<Tm

# In[1]:

from itertools import repeat
import logging
from typing import Callable, Dict, Iterable, List, Set
import dill
import numpy as np
import numpy.typing as npt
import time
from multiprocessing import cpu_count
from pebble.concurrent.process import TimeoutError
import pebble
from random_sampler import RV
import copy

simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]


class SimulatedAnnealing():
    def __init__(self, 
                 simulator: Callable[[candidateType], simResultType], 
                 priors : priorType, 
                 min_epsilon: float,
                 distance_function: Callable[[distanceArgType, distanceArgType], float],
                 Yobs: distanceArgType, 
                 outfile: str,
                 cores: int = cpu_count(),
                 maxiter: int = 100000, 
                 rng: np.random.Generator = None, 
                 save_intermediate: bool = False,
                 
                 generation_size: int = 128,
                 cooling_rate: float = 0.95,
                 initial_temp: float = 100,
                 final_temp: float = 0.1 
                 ):
        """Implements the Simulated Annealing algorithm designed to detect multiple optima in the fitness landscape

        Args:
            simulator (Callable[[candidateType], simResultType]): a function that takes a dictionary of parameters as input. Ouput {'data':Ysim}
            priors (priorType): a dictionary which use id of parameters as keys and RV class object as values
            min_epsilon (float): minimal epsilon
            distance_function (Callable[[distanceArgType, distanceArgType], float]): a function that calculate the distance between observed data and simulated data
            Yobs (distanceArgType): observed data
            outfile (str): unique id for the experiment. This will be also used to continue a simulation that 
                         is partly done
            cores (int, optional): number of treads. Defaults to cpu_count().
            generation_size (int, optional): the size of each population. Defaults to 128.
            maxiter (int, optional): The maximum number of generations to simulate. Defaults to 100000.
            rng (np.random.Generator, optional): Random number generator in order to ensure reproducibility. Defaults to None.
            n_children (int, optional): The number of children of generate for each generation. Defaults to 2.
            scaling_factor (float, optional): The scaling factor for the Differential Evolution, meaning how much weight is applied to the two secondary individuals. Defaults to 0.5.
            crossover_prob (float, optional): The crossover probability for the two secondary parents. Defaults to 0.9.
            save_intermediate (bool, optional): Should intermediate results be saved for each iteration? If toggled on, computations can be resumed if interrupted prematurly, but this will come at a performance penalty which parallelism cannot alleviate.. Defaults to False.
            
            !!! Distance is to be minimized!!!
        """
        self.simulator = simulator
        self.priors = priors
        self.parameter_names = list(priors.keys())
        self.n_parameters = len(self.parameter_names)
        self.param_std: Dict[str, float] = {param: prior.scale for param, prior in priors.items()}
        if rng is None:
            default_seed = 1952
            self.rng = np.random.Generator(np.random.PCG64(default_seed))
        else:
            self.rng = rng
        self.distance_function = distance_function
        self.min_epsilon = min_epsilon
        self.Yobs = Yobs
        self.outfile = outfile
        # Compared to SMC-ABC this seems a bit odd. The reationale is indirection.
        # The indicies in the list specify which of the particles in self.all_particles is part of the current population
        self.population: List[candidateType] = []
        self.cores = cores    
        self.epsilons: List[float] = []          # min distance in each generation
        self.generation_size = generation_size   # number of particles to be simulated at each generation
        self.maxiter = maxiter
        # NOTE: all_simulated_data, all_distances, all_particles and birth_generation MUST be aligned
        self.all_simulated_data: List[simResultType] = []
        self.all_distances: List[float] = []
        self.all_particles: List[candidateType] = []
        # The number of times the particle has been challenged by an offspring
        self.times_challenged: List[int] = []
        # Specifies the generation each of the particles are born in
        self.birth_generation: List[int] = []
        self.generation = 0
        self.save_intermediate = save_intermediate

        self.cooling_rate = cooling_rate
        self.initial_temp = initial_temp
        self.current_temp = initial_temp
        self.final_temp = final_temp

    def choose_particle_2(self, particle1: candidateType, particle2: candidateType)-> candidateType:
        index_p1 = self.get_index(particle1)
        index_p2 = self.get_index(particle2)

        d1 = self.all_distances[index_p1]
        d2 = self.all_distances[index_p2]
        if self.energy_function(d1, d2):
            return particle2
        else:
            return particle1
        
    def get_index(self, particle):
        return self.all_particles.index(particle)
        
    def change_all_parameters_2(self, current_population: List[candidateType]):
        logging.info("Changing all parameters of a particle")
        new_population = copy.deepcopy(current_population)
        for particle in new_population:
            for key in particle:
                old_parameter_value = particle[key]
                particle[key] += 0.1 * self.rng.normal(0, 1)
                if not self.check_validity(particle, key):
                    #self.all_particles.append(particle) # Skal kanskje ikke være her, men i evaluate_candidates
                    particle[key] = old_parameter_value

                # else:
                #     # This new parameter value violates our constraints, so we must revert the change
                #     particle[key] = old_parameter_value
        return new_population
    
    def change_all_parameters(self, old_particle):
        logging.info("Changing all parameters of a particle")
        new_particle = copy.deepcopy(old_particle)
        for key in new_particle:
            old_parameter_value = new_particle[key]
            new_particle[key] += 0.1 * self.rng.normal(0, 1)
            if not self.check_validity(new_particle,key):
                # This new parameter value violates our constraints, so we must revert the change
                new_particle[key] = old_parameter_value
        return new_particle
    

    def generate_candidates_2(self, indices):
        candidates: List[candidateType] = []
        for index in indices:
            candidate = {parameter: value for parameter, value in self.all_particles[index].items()}
            for key in candidate:
                old_parameter_value = candidate[key]
                candidate[key] += 0.1 * self.rng.normal(0, 1)
                if not self.check_validity(candidate, key):
                    #self.all_particles.append(particle) # Skal kanskje ikke være her, men i evaluate_candidates
                    candidate[key] = old_parameter_value
            candidates.append(candidate)
        logging.info("Evaluating fitness of candidates")
        self.evaluate_candidates(candidates)

    def generate_candidates_old2(self, indices): #Unused
        candidates: List[candidateType] = []
        for index in indices:
            # if index is None:
            #     continue
            candidate = {parameter: value for parameter, value in self.all_particles[index].items()}
            for key in candidate:
                old_parameter_value = candidate[key]
                candidate[key] += 0.1 * self.rng.normal(0, 1)
                if not self.check_validity(candidate, key):
                    #self.all_particles.append(particle) # Skal kanskje ikke være her, men i evaluate_candidates
                    candidate[key] = old_parameter_value
            candidates.append(candidate)
        logging.info("Evaluating fitness of candidates")
        self.evaluate_candidates(candidates)


    def replace_population_2(self,original_population: npt.NDArray[np.int64], candidates: npt.NDArray[np.int64]):
       logging.info(f"Replacing population with children")
       current_population : Set[int] = set()
       for i, (old_particle, candidate_particle) in enumerate(zip(original_population, candidates)):
            chosen_particle = self.choose_particle(old_particle, candidate_particle)
            if chosen_particle == old_particle:
                current_population.add(old_particle)
                #Keep index in current population as it is
                self.inner_iterations_list[i] += 1 #Add 1 to inner_iterations_list
                logging.info(f"Keeping particle {i}") #For checking, remove later
                
            else:
                current_population.add(candidate_particle)
                #Add index of this particle to self.population
                self.inner_iterations_list[i] = self.min_layers #Set inner_iterations_list to 1
                #self.is_improved_list[i] = True
                logging.info(f"Replacing particle {i}") #For checking, remove later
       logging.info(f"Updating population")
       self.population.append(list(current_population))

    def replace_population_4(self,original_population: npt.NDArray[np.int64], candidates, candidates_idxs): #Unused
       logging.info(f"Replacing population with children")
       current_population : Set[int] = set()
       for i,  candidate_particle in zip(candidates_idxs, candidates):
            old_particle = original_population[i]
            chosen_particle = self.choose_particle(old_particle, candidate_particle)
            if chosen_particle == old_particle:
                current_population.add(old_particle)
                #Keep index in current population as it is
                self.inner_iterations_list[i] += 1 #Add 1 to inner_iterations_list
                logging.info(f"Keeping particle {i}") 
                
            else:
                current_population.add(candidate_particle)
                #Add index of this particle to self.population
                self.inner_iterations_list[i] = 1 #Set inner_iterations_list to 1
                self.is_improved_list[i] = True
                logging.info(f"Replacing particle {i}") 
       logging.info(f"Updating population")
       self.population.append(list(current_population))
    
    def simulate_generation_2(self):
        current_population = np.array(list(self.population[-1]))
        self.generate_candidates_2(current_population)
        candidates_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
        self.replace_population_2(current_population, candidates_idxs)

        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")
        

    def simulate_generation_3(self):
        current_index_population = np.array(list(self.population[-1]))
        #To be replaced by generate_candidates()
        current_particles = [self.all_particles[x] for x in current_index_population]
        candidates = self.change_all_parameters_2(current_particles)
        self.evaluate_candidates(candidates)
        #----


        new_index_population = []
        for i, (old_particle, candidate_particle) in enumerate(zip(current_particles, candidates)):
            chosen_particle = self.choose_particle_2(old_particle, candidate_particle)
            if chosen_particle == old_particle:
                new_index_population.append(self.get_index(old_particle))
                #Keep index in current population as it is
                self.inner_iterations_list[i] += 1 #Add 1 to inner_iterations_list
                logging.info(f"Keeping particle {i}") #For checking, remove later
                
            else:
                new_index_population.append(self.get_index(candidate_particle))
                #Add index of this particle to self.population
                self.inner_iterations_list[i] = 1 #Set inner_iterations_list to 1
                logging.info(f"Replacing particle {i}")

        self.population.append(new_index_population)
        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")


                



    def simulate_generation_old(self):
        if self.cores == 1:
            start = time.time()
            for particle in self.population_old:
                for _ in range(self.inner_iterations): #rewrite with generator to get rid of nested loop
                    new_particle = self.change_all_parameters(particle)
                    chosen_particle = self.choose_particle_2(particle, new_particle)
                    if chosen_particle != particle:
                        # We have a new candidate
                        index = self.population_old.index(particle)
                        self.population_old[index] = chosen_particle
                        #self.birth_generation[index] = (self.generation)
                        #self.times_challenged[index] = 0 #Trenger ikke denne?
                        break
            
        else:
            with pebble.ProcessPool(self.cores) as p:
                res_map = p.map(self.simulator, self.population_old)
            simulated_data = list(res_map) # type: ignore

    def run_simulation_old(self) -> None:
        # Ensures random state is respected
        for item in self.priors.values():
            item.set_rng(rng=self.rng)

        logging.info(f"")
       
        logging.info(f"Generating initial population with {self.generation_size} particles")
        self.population_old = [self.generator() for _ in range(self.generation_size)]
        logging.info(f"Evaluating initial population")
        for particle in self.population_old:
            distance, simulated_data = self.evaluate_single_candidate(particle)
            self.all_distances.append(distance)
            self.all_simulated_data.append(simulated_data)
    
        max_generation_epsilon = max(self.all_distances)
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")
        self.generation += 1
        self.current_temp *= self.cooling_rate

        while self.generation <= self.maxiter:
            if max_generation_epsilon < self.min_epsilon:
                logging.info(f"Fitness objective reached at generation {self.generation}")
                logging.info(f"Exiting simulated annealing")
                break
            if self.current_temp < self.final_temp:
                logging.info(f"Temperature has reached a minimum at generation {self.generation}. Fitness objective not reached.")
                logging.info(f"Exiting simulated annealing")
                #break
            
            logging.info(f"Running generation {self.generation} of {self.maxiter}. Current temperature: {self.current_temp}")
            self.simulate_generation_old()

            self.generation += 1
            self.current_temp *= self.cooling_rate
            if self.save_intermediate:
                    dill.dump(self,open(self.outfile,'wb'))
        else:
            # This else-clause belongs to the main simulated annealing loop
            logging.info("Fitness objective not reached after maximum number of generations")
            logging.info("Exiting simulated annealing")
        
        logging.info(f"Saving results to {self.outfile}")
        dill.dump(self, file=open(self.outfile,mode='wb'))


    def evaluate_single_candidate(self, candidate: candidateType): #Unused/don't want to use this function
        # Specifying timeout of 30 minutes
        timeout = 30*60
        # This function both evaluates newly born individuals and store them into the archive
        start = time.time()
        # Candidates for which evaluating fitness was successful
        successfull_candidates = set()
        candidate_counter = 0
        # No need for creating a parallel cluster in this case
        try:
            simulated_data: simResultType = self.simulator(candidate)
            logging.info("Evaluation of candidate ran successfully")
            success = True
        except Exception as e:
            logging.error(f"Candidate evaluation failed: {e}")

        distance = self.distance_function(self.Yobs, simulated_data)

        return distance, simulated_data
    

    def get_distance(self, particle):
        return self.all_distances[self.all_particles.index(particle)]
