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
                 final_temp: float = 0.1,
                 inner_iterations: int = 1
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
        self.population2: List[List[int]] = []
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
        for particle in self.population:
            #particle = self.all_particles[particle_idx]
            for p,val in particle.items(): 
                lst = parameters.get(p,[])
                lst.append(val)
                parameters[p] = lst        
        for p, lst in parameters.items():
            self.param_std[p] = np.std(lst)

    def evaluate_candidates(self, candidates: List[candidateType]):
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
                                raw_res = res_iter.next()
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
        self.birth_generation.extend(repeat(self.generation,len(simulated_data)))
        #The four lines above are in use. The one lines below are not. Must look further into it, maybe remove, maybe implement, maybe add more lines 
        self.times_challenged.extend(repeat(0,len(simulated_data))) #This is not in use in evo_etc either, just recorded as information. Or not updated either it seems
        end = time.time()
        logging.debug('Completed parallel evaluation of candidates in {0} seconds'.format(end - start))
        logging.debug(f"Length of all_simulated_data is {len(self.all_simulated_data)}")
        logging.debug(f"Length of all_distances is {len(self.all_distances)}")
        logging.debug(f"Length of all_particles is {len(self.all_particles)}")
        #return distances, simulated_data


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


    def particle_distance(self, idx_1: int, idx_2: int)-> float:
        epsilon = 1e-8
        particle_1 = self.all_particles[idx_1]
        particle_2 = self.all_particles[idx_2]
        return sum([(particle_1[key] - particle_2[key])**2 / (self.param_std[key]**2 + epsilon) for key in self.param_std.keys()])




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
    
    def get_distance(self, particle):
        return self.all_distances[self.all_particles.index(particle)]

    def energy_function(self, d1, d2):
        if d2 < d1:
            return True
        else:
            return self.rng.random() < np.exp(-(d2-d1)/self.current_temp)
        
    def get_index(self, particle):
        return self.all_particles.index(particle)
        
    
    

    def simple_evaluation(self, particle1, particle2): #Change name, not simple anymore. Compare_particles() f ex
        index_p1 = self.get_index(particle1)
        index_p2 = self.get_index(particle2)

        d1 = self.all_distances[index_p1]
        d2 = self.all_distances[index_p2]
        if self.energy_function(d1, d2):
            return particle2
        else:
            return particle1
        
    def simple_evaluation_2(self, particle1, particle2): #Change name, not simple anymore. Compare_particles() f ex
        d1 = self.all_distances[particle1]
        d2 = self.all_distances[particle2]
        if self.energy_function(d1, d2):
            return particle2
        else:
            return particle1
        


    def update_population(self, original_population: npt.NDArray[np.int64]): #Brukes ikke
        logging.info(f"Applying Simulated Annealing to current population")
        new_population: Set[int] = set(original_population)
        for particle in original_population:
            new_particle = self.change_all_parameters(particle)
            chosen_particle = self.simple_evaluation(particle, new_particle)
            new_population.remove(particle)
            new_population.add(chosen_particle)
        return

    def generate_candidates(self, indices):
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

    def generate_candidates_2(self, indices):
        candidates: List[candidateType] = []
        for index in indices:
            if index is None:
                continue
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

    def generate_candidates_4(self, particle_indices):
        candidates: List[candidateType] = []
        candidate_indices = []
        for i, index in enumerate(particle_indices):
            if index is None:
                continue
            candidate = {parameter: value for parameter, value in self.all_particles[index].items()}
            for key in candidate:
                old_parameter_value = candidate[key]
                candidate[key] += 0.1 * self.rng.normal(0, 1)
                if not self.check_validity(candidate, key):
                    #self.all_particles.append(particle) # Skal kanskje ikke være her, men i evaluate_candidates
                    candidate[key] = old_parameter_value
            candidates.append(candidate)
            candidate_indices.append(i)
        logging.info("Evaluating fitness of candidates")
        self.evaluate_candidates(candidates)
        return candidates, candidate_indices

    def replace_population(self,original_population: npt.NDArray[np.int64], candidates: npt.NDArray[np.int64]):
       logging.info(f"Replacing population with children")
       current_population : Set[int] = set()
       for i, (old_particle, candidate_particle) in enumerate(zip(original_population, candidates)):
            chosen_particle = self.simple_evaluation_2(old_particle, candidate_particle)
            if chosen_particle == old_particle:
                current_population.add(old_particle)
                #Keep index in current population as it is
                self.inner_iterations_list[i] += 1 #Add 1 to inner_iterations_list
                logging.info(f"Keeping particle {i}") #For checking, remove later
                
            else:
                current_population.add(candidate_particle)
                #Add index of this particle to self.population
                self.inner_iterations_list[i] = 1 #Set inner_iterations_list to 1
                #self.is_improved_list[i] = True
                logging.info(f"Replacing particle {i}") #For checking, remove later
       logging.info(f"Updating population")
       self.population2.append(list(current_population))

    def replace_population_4(self,original_population: npt.NDArray[np.int64], candidates, candidates_idxs):
       logging.info(f"Replacing population with children")
       current_population : Set[int] = set()
       for i,  candidate_particle in zip(candidates_idxs, candidates):
            old_particle = original_population[i]
            chosen_particle = self.simple_evaluation_2(old_particle, candidate_particle)
            if chosen_particle == old_particle:
                current_population.add(old_particle)
                #Keep index in current population as it is
                self.inner_iterations_list[i] += 1 #Add 1 to inner_iterations_list
                logging.info(f"Keeping particle {i}") #For checking, remove later
                
            else:
                current_population.add(candidate_particle)
                #Add index of this particle to self.population
                self.inner_iterations_list[i] = 1 #Set inner_iterations_list to 1
                self.is_improved_list[i] = True
                logging.info(f"Replacing particle {i}") #For checking, remove later
       logging.info(f"Updating population")
       self.population2.append(list(current_population))
    
        
    def simulate_generation_3(self):
        current_population = np.array(list(self.population2[-1]))
        self.generate_candidates(current_population)
        candidates_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
        self.replace_population(current_population, candidates_idxs)

        max_generation_epsilon = max(self.all_distances[p] for p in self.population2[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")

    def particles_to_idxs(self, particles):
        idxs = [self.get_index(particle) for particle in particles]
        return idxs


    def simulate_generation_4(self):
        current_population = np.array(list(self.population2[-1]))
        max_layer = min(max(self.inner_iterations_list), 10)
        self.is_improved_list = [False for _ in range(self.generation_size)]
        logging.info(f"Max layer is {max_layer}")
        for layer in range(1, max_layer+1):
            logging.info(f"Checking layer {layer}")
            filtered_candidates = []
            filtered_indices = []
            for i, ID in enumerate(current_population):
                if layer <= self.inner_iterations_list[i] and not self.is_improved_list[i]:
                    filtered_candidates.append(ID)
                    filtered_indices.append(i)
                    logging.info(f"Another candidate for index {i}")
                else:
                    filtered_candidates.append(None)
                    filtered_indices.append(None)
            candidates, candidates_idxs = self.generate_candidates_4(filtered_candidates)
            candidates = self.particles_to_idxs(candidates)
            #candidates_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
            #logging.info(f"Length of candidates_idxs is {len(candidates_idxs)}")
            self.replace_population_4(current_population, candidates, candidates_idxs)

        

        max_generation_epsilon = max(self.all_distances[p] for p in self.population2[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")


        

    def simulate_generation_2(self):
        current_index_population = np.array(list(self.population2[-1]))
        #To be replaced by generate_candidates()
        current_particles = [self.all_particles[x] for x in current_index_population]
        candidates = self.change_all_parameters_2(current_particles)
        self.evaluate_candidates(candidates)
        #----


        new_index_population = []
        for i, (old_particle, candidate_particle) in enumerate(zip(current_particles, candidates)):
            chosen_particle = self.simple_evaluation(old_particle, candidate_particle)
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

        self.population2.append(new_index_population)
        max_generation_epsilon = max(self.all_distances[p] for p in self.population2[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        logging.info(f"Model epsilon {max_generation_epsilon}")


                



    def simulate_generation(self):
        if self.cores == 1:
            start = time.time()
            for particle in self.population:
                for _ in range(self.inner_iterations): #rewrite with generator to get rid of nested loop
                    new_particle = self.change_all_parameters(particle)
                    chosen_particle = self.simple_evaluation(particle, new_particle)
                    if chosen_particle != particle:
                        # We have a new candidate
                        index = self.population.index(particle)
                        self.population[index] = chosen_particle
                        #self.birth_generation[index] = (self.generation)
                        #self.times_challenged[index] = 0 #Trenger ikke denne?
                        break
            
        else:
            with pebble.ProcessPool(self.cores) as p:
                res_map = p.map(self.simulator, self.population)
            simulated_data = list(res_map)


    def run_simulation_2(self) -> None:
        # Ensures random state is respected
        for item in self.priors.values():
            item.set_rng(rng=self.rng)

        logging.info(f"")
        if self.generation == 0:
            logging.info(f"Generating initial population with {self.generation_size} particles")
            initial_population = [self.generator() for _ in range(self.generation_size)]
            logging.info(f"Evaluating initial population")
            #distances, simulated_data = self.evaluate_candidates(initial_population)
            self.evaluate_candidates(initial_population)

            # self.all_particles.extend(initial_population) #These lines are included inside self.evaluate_candidates()
            # self.all_distances.extend(distances)
            # self.all_simulated_data.extend(simulated_data)
            self.population2.append(list(range(len(self.all_particles))))


            max_generation_epsilon = max(self.all_distances)
            self.epsilons.append(max_generation_epsilon)
            self.update_std()
            self.inner_iterations_list = [1 for _ in range(self.generation_size)]
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

            #self.simulate_generation_2()
            #self.simulate_generation_3()
            self.simulate_generation_4()
            logging.info("Using version 4")


            self.generation += 1
            self.current_temp *= self.cooling_rate
            logging.debug(f" Inner iterations_list: {self.inner_iterations_list}")
            if self.save_intermediate:
                    dill.dump(self,open(self.outfile,'wb'))
        else:
            # This else-clause belongs to the main simulated annealing loop
            logging.info("Fitness objective not reached after maximum number of generations")
            logging.info("Exiting simulated annealing")

        logging.info(f"Saving results to {self.outfile}")
        dill.dump(self, file=open(self.outfile,mode='wb'))
        

    def run_simulation(self) -> None:
        # Ensures random state is respected
        for item in self.priors.values():
            item.set_rng(rng=self.rng)

        logging.info(f"")
       
        logging.info(f"Generating initial population with {self.generation_size} particles")
        self.population = [self.generator() for _ in range(self.generation_size)]
        logging.info(f"Evaluating initial population")
        for particle in self.population:
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
            self.simulate_generation()

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
