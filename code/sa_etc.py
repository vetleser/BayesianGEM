#!/usr/bin/env python
# coding: utf-8

# #### Add constraint: Topt<Tm

# In[1]:

from itertools import repeat
import logging
from typing import Callable, Dict, Iterable, List, Set, Tuple, Optional
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
                 rng: np.random.Generator = None,  # type: ignore
                 save_intermediate: bool = False,
                 
                 generation_size: int = 128,
                 cooling_rate: float = 0.95,
                 initial_temp: float = 100,
                 final_temp: float = 0.1,
                 inner_iterations: int = 1,
                 min_layers: int = 1,
                 max_layers: int = 10,
                 version: int = 1,
                 step_size: float = 0.1
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
        self.population_old: List[candidateType] = []
        self.population: List[List[int]] = []
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
        #self.cooling_rate = (final_temp / initial_temp) ** (1 / maxiter) #cooling rate lines up with iterations
        self.initial_temp = initial_temp
        self.current_temp = initial_temp
        self.final_temp = final_temp
        self.inner_iterations = inner_iterations
        self.inner_iterations_list: List[int] = []
        self.is_improved_list : List[bool] = []
        self.current_layer : int = 1
        self.birth_generation_layer : List[Tuple[int,int]] = []
        self.min_layers = min_layers
        self.max_layers = max_layers

        self.acceptance_rates: List[float] = []
        self.version = version

        self.log_ef = True
        self.log_ef_list: List[float] = []
        self.step_size = step_size


    def generator(self) -> candidateType:
            candidate = {param: self.priors[param].rvfv() for param in self.parameter_names}
            [self.correct_validity(candidate=candidate,entry=entry) for entry in candidate]
            return candidate


    def correct_validity(self, candidate, entry: str) -> None:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt > Tm > 0 (in Kelvin of course) in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')
        # We assume that entries are of the form PROTID_{Tm,Topt,dCpt}
        # If this is not the case, we assume that the algorithm is used for another kind of inference problem,
        # so we skip this domain-specific check. This also applies to the dCPt as mutatating them does not violate the constraint
        if len(split_entry) != 2 or split_entry[1] not in ("Tm","Topt"):
            return
        protein_id = split_entry[0]
        Tm_key = protein_id + "_Tm"
        Topt_key = protein_id + "_Topt"
        for _ in range(10):
            # We try to get things right 10 times before we give up
            Tm = candidate[Tm_key]
            Topt = candidate[Topt_key]
            if not Tm > Topt > 0:
                self.mutate_param(candidate,Tm_key)
                self.mutate_param(candidate,Topt_key)
            else:
                break
    
    
    def update_std(self):
        """
        This method updates the enzyme parameter preceived standard deviations. Corresponds to update_posterior() in abc_etc 
        """
        logging.info('Updating standard deviations to parameters')
        parameters = dict()   # {'Protein_Tm':[]}
        for particle_idx in self.population[-1]:
            particle = self.all_particles[particle_idx]
            for p,val in particle.items(): 
                lst = parameters.get(p,[])
                lst.append(val)
                parameters[p] = lst        
        for p, lst in parameters.items():
            self.param_std[p] = np.std(lst)

    def evaluate_candidates(self, candidates: List[candidateType]) -> None: #Taken straight from evo_etc
        # Specifying timeout of 30 minutes
        timeout = 30*60
        # This function both evaluates newly born individuals and store them into the archive
        start = time.time()
        simulated_data = []
        # Candidates for which evaluating fitness was successful
        successfull_candidates = set()
        candidate_counter = 0
        if self.cores == 1:
            # No need for creating a parallel cluster in this case
            res_iter: Iterable[simResultType] = map(self.simulator,candidates)
            while True:
                try:
                    raw_res = res_iter.__next__()
                except StopIteration:
                    # We have now iterated over all particles
                    break
                else:
                    logging.info("Evaluation of particle ran successfully")
                    simulated_data.append(raw_res)
                    successfull_candidates.add(candidate_counter)
                finally:
                    candidate_counter += 1
        else:
                try:
                    with pebble.ProcessPool(self.cores) as p:
                        #Wrapping simulator
                        serialized_simulator = dill.dumps(self.simulator)
                        simulator_func = dill.loads(serialized_simulator)
                        #Old attempt
                        res_iter: Iterable[simResultType] = p.map(self.simulator, candidates,timeout=timeout).result()
                        while True:
                            try:
                                raw_res = res_iter.next() # type: ignore
                            except StopIteration:
                                # We have now iterated over all particles
                                break
                            except TimeoutError:
                                logging.info("Evaluation of particle timed out")
                            else:
                                logging.info("Evaluation of particle ran successfully")
                                simulated_data.append(raw_res)
                                successfull_candidates.add(candidate_counter)
                            finally:
                                candidate_counter += 1
                except (OSError, RuntimeError) as e:
                    logging.error('failed parallel_evaluation_mp: {0}'.format(str(e)))
                    raise
        
        distances = [self.distance_function(self.Yobs, res) for res in simulated_data]

        # save all simulated results
        self.all_simulated_data.extend(simulated_data)
        self.all_distances.extend(distances)
        # This deals with the problem of candidates failing evaluation
        self.all_particles.extend([candidate for counter, candidate in enumerate(candidates) if counter in successfull_candidates])
        self.birth_generation.extend(repeat(self.generation,len(simulated_data))) #Can probably be removed, replaced by birth_generation_layer

        self.birth_generation_layer.extend(repeat((self.generation,self.current_layer),len(simulated_data))) 

        #The four lines above are in use. The one lines below are not. Must look further into it, maybe remove, maybe implement, maybe add more lines 
        self.times_challenged.extend(repeat(0,len(simulated_data))) #This is not in use in evo_etc either, just recorded as information. Or not updated either it seems
        end = time.time()
        logging.info('Completed parallel evaluation of candidates in {0} seconds'.format(end - start))
        logging.debug(f"Length of all_simulated_data is {len(self.all_simulated_data)}")
        logging.debug(f"Length of all_distances is {len(self.all_distances)}")
        logging.debug(f"Length of all_particles is {len(self.all_particles)}")




    def mutate_param(self,candidate, entry: str):
            candidate[entry] = self.priors[entry].rvfv()

    def check_validity(self, candidate, entry: str) -> bool:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt > Tm in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')
        # We assume that entries are of the form PROTID_{Tm,Topt,dCpt}
        # If this is not the case, we assume that the algorithm is used for another kind of inference problem,
        # so we skip this domain-specific check. This also applies to the dCPt as mutatating them does not violate the constraint
        if len(split_entry) != 2 or split_entry[1] not in ("Tm","Topt"):
            return True
        protein_id = split_entry[0]
        Tm_key = protein_id + "_Tm"
        Topt_key = protein_id + "_Topt"
        Tm = candidate[Tm_key]
        Topt = candidate[Topt_key]
        return Tm > Topt > 0



    def energy_function(self, d1: float, d2: float) -> bool:
        if self.log_ef:
            logging.debug(f"Energy function is : {np.exp(-(d2-d1)/self.current_temp)}")
            self.log_ef_list.append(np.exp(-(d2-d1)/self.current_temp))
            self.log_ef = False
        if d2 < d1:
            return True
        else:
            return self.rng.random() < np.exp(-(d2-d1)/self.current_temp)
        

        
    def choose_particle(self, idx1: int, idx2: int)-> int: 
        d1 = self.all_distances[idx1]
        d2 = self.all_distances[idx2]
        if self.energy_function(d1, d2):
            return idx2
        else:
            return idx1
        


    def generate_candidates(self, particle_idxs: npt.NDArray[np.int64])-> None: #Kan bruke change_all_parameters her istedenfor å skrive det eksplisitt
        logging.info("Generating candidates")
        candidates: List[candidateType] = []
        for idx in particle_idxs.tolist():
            candidate = {parameter: value for parameter, value in self.all_particles[idx].items()}
            for key in candidate:
                old_parameter_value = candidate[key]
                candidate[key] += self.step_size * self.rng.normal(0, 1) #Endre til å bruke change_all_parameters, og/eller måte på å endre verdiene
                if not self.check_validity(candidate, key):
                    candidate[key] = old_parameter_value
            candidates.append(candidate)
        logging.info("Evaluating fitness of candidates")
        self.evaluate_candidates(candidates)


    def replace_population(self,original_population: npt.NDArray[np.int64], candidates: npt.NDArray[np.int64], filtered_indices: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
        logging.info(f"Replacing population with candidates")
        current_population = original_population.copy()
        for candidate, i in zip(candidates, filtered_indices):
            original_particle = original_population[i]
            chosen_particle = self.choose_particle(original_particle, candidate)
            if chosen_particle == original_particle:
                current_population[i] = original_particle
                if self.current_layer == 1:
                    self.inner_iterations_list[i] += 1
                self.is_improved_list[i] = False
                logging.debug(f"Keeping particle {i}")
            else:
                current_population[i] = candidate
                self.inner_iterations_list[i] = self.min_layers
                self.is_improved_list[i] = True
                logging.debug(f"Replacing particle {i}")

        return current_population

        


    def filter_candidates(self, current_population: npt.NDArray[np.int64])-> Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]:
        logging.info("Finding candidates to improve") 
        improvable_candidates: List[int] = []
        improvable_indices: List[int] = []
        for i, ID in enumerate(current_population):
            if self.current_layer <= self.inner_iterations_list[i] and not self.is_improved_list[i]:
                improvable_candidates.append(ID)
                improvable_indices.append(i)
                if self.current_layer > 1:
                    logging.info(f"Found particle to improve at index {i} in layer {self.current_layer}")

        return (
                np.array(improvable_candidates, dtype=np.int64),
                np.array(improvable_indices, dtype=np.int64)
                )

    def simulate_generation(self) -> None: #Kan lage en if layer ==1, og resten i en annen løkke. Dropper det, funker nå
        current_population = np.array(list(self.population[-1]))
        current_max_layer = min(max(self.inner_iterations_list), self.max_layers) #Akkurat nå er max_layer for en partikkel begrenset av inner_iterations_list. Må kanskje endre det
        current_max_layer = max(current_max_layer, self.min_layers)
        self.is_improved_list = [False for _ in range(self.generation_size)]
        self.current_layer = 1
        logging.info(f"Max layer is {current_max_layer}")
        for layer in range(1, current_max_layer+1):
            self.current_layer = layer
            logging.info(f"Checking layer {layer}")
            particles_to_improve, positions = self.filter_candidates(current_population) #Check if particles should be improved, and keep track of which particles and where they are in current_population. Replacement for commented lines below
            self.generate_candidates(particles_to_improve)
            candidates_idxs = np.array([
                i for i, (gen, layer) in enumerate(self.birth_generation_layer)
                if gen == self.generation and layer == self.current_layer
            ], dtype=np.int64)
            logging.debug(f"Candidates idxs: {candidates_idxs}")
            current_population = self.replace_population(current_population, candidates_idxs, positions)

        logging.info(f"Current population is {current_population}")
        self.population.append(list(current_population))

        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")

    

    def replace_population_2(self,original_population: npt.NDArray[np.int64], candidates: npt.NDArray[np.int64]):
        logging.info(f"Replacing population with candidates")
        current_population : Set[int] = set()
        acceptance_counter: int = 0
        for i, (old_particle, candidate_particle) in enumerate(zip(original_population, candidates)):
            chosen_particle = self.choose_particle(old_particle, candidate_particle)
            if chosen_particle == old_particle:
                current_population.add(old_particle)
                #Keep index in current population as it is
                self.inner_iterations_list[i] += 1 #Add 1 to inner_iterations_list
                logging.debug(f"Keeping particle {i}") #For checking, remove later
                
            else:
                current_population.add(candidate_particle)
                #Add index of this particle to self.population
                acceptance_counter += 1
                self.inner_iterations_list[i] = self.min_layers #Set inner_iterations_list to 1
                #self.is_improved_list[i] = True
                logging.debug(f"Replacing particle {i}") #For checking, remove later
        acceptance_rate: float = acceptance_counter / len(candidates)
        logging.info(f"Acceptance rate in generation {self.generation}: {acceptance_rate}")
        self.acceptance_rates.append(acceptance_rate)
        self.population.append(list(current_population))

    def generate_candidates_2(self, indices):
        candidates: List[candidateType] = []
        for index in indices:
            candidate = {parameter: value for parameter, value in self.all_particles[index].items()}
            for key in candidate:
                old_parameter_value = candidate[key]
                candidate[key] += self.step_size * self.rng.normal(0, 1)
                if not self.check_validity(candidate, key):
                    #self.all_particles.append(particle) # Skal kanskje ikke være her, men i evaluate_candidates
                    candidate[key] = old_parameter_value
            candidates.append(candidate)
        logging.info("Evaluating fitness of candidates")
        self.evaluate_candidates(candidates)

    def simulate_generation_2(self):
        current_population = np.array(list(self.population[-1]))
        self.generate_candidates_2(current_population)
        candidates_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
        self.replace_population_2(current_population, candidates_idxs)

        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")
        


    def run_simulation(self) -> None:
        # Ensures random state is respected
        for item in self.priors.values():
            item.set_rng(rng=self.rng)

        logging.info(f"Using version {self.version} of Simulated Annealing")
        if self.generation == 0:
            logging.info(f"Generating initial population with {self.generation_size} particles")
            initial_population = [self.generator() for _ in range(self.generation_size)]
            logging.info(f"Evaluating initial population")
            self.evaluate_candidates(initial_population)


            self.population.append(list(range(len(self.all_particles))))


            max_generation_epsilon = max(self.all_distances)
            self.epsilons.append(max_generation_epsilon)
            self.update_std()
            self.inner_iterations_list = [self.min_layers for _ in range(self.generation_size)]
            logging.info(f"Model epsilon {max_generation_epsilon}")
            self.generation += 1

        while self.generation <= self.maxiter:
            start_generation = time.time()
            if max_generation_epsilon < self.min_epsilon:
                logging.info(f"Fitness objective reached at generation {self.generation}")
                logging.info(f"Exiting simulated annealing")
                break
            if self.current_temp < self.final_temp:
                logging.info(f"Temperature has reached a minimum at generation {self.generation}. Fitness objective not reached.")
                logging.info(f"Exiting simulated annealing")
                #break
            
            logging.info(f"Running generation {self.generation} of {self.maxiter}. Current temperature: {self.current_temp}")
            self.log_ef = True
            if self.version == 1:
                self.simulate_generation()
            elif self.version == 2:
                self.simulate_generation_2()
            # else:
            #     self.simulate_generation_3()
            
            
            end_generation = time.time()
            logging.info(f"Generation {self.generation} completed in {end_generation-start_generation} seconds")

            self.generation += 1
            if self.current_temp > self.final_temp:
                self.current_temp *= self.cooling_rate
            logging.info(f" Inner_iterations_list: {self.inner_iterations_list}")
            if self.save_intermediate:
                    dill.dump(self,open(self.outfile,'wb'))
        else:
            # This else-clause belongs to the main simulated annealing loop
            logging.info("Fitness objective not reached after maximum number of generations")
            logging.info("Exiting simulated annealing")
            logging.info(f"Self.log_ef_list: {self.log_ef_list}")
            logging.info(f"Self.acceptance_rates: {self.acceptance_rates}")
            logging.info(f"Final temperature: {self.current_temp}")

        logging.info(f"Saving results to {self.outfile}")
        dill.dump(self, file=open(self.outfile,mode='wb'))
        


