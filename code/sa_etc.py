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
                 cooling_rate: float = None,
                 initial_temp: float = 100,
                 final_temp: float = 0.1,
                 inner_iterations: int = 1,
                 min_layers: int = 1,
                 max_layers: int = 10,
                 version: int = 1,
                 temp_step_size: float = 0.1,
                 normalize: bool = False,
                 move_type : str= 'normal',
                 dCpt_step_size: float = 1.0,
                 step_size: float = 0.1,  # This is the step size for the normal move type
                 end_exploration: float = 1.0,  # This is the fraction of the maximum number of generations that is used for exploration
                 initial_step_size: float = 10,  # This is the initial step size for the normal move type
                 final_step_size: float = 0.1,  # This is the final step size for the normal move type, after cooling down
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

        # self.cooling_rate = cooling_rate
        if cooling_rate is None:
            self.cooling_rate = (final_temp / initial_temp) ** (1 / (maxiter*end_exploration)) #cooling rate lines up with iterations
        else:
            self.cooling_rate = cooling_rate
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
        self.step_size = step_size  # This is the step size for the normal move type
        self.temp_step_size = temp_step_size
        self.dCpt_step_size = dCpt_step_size
        self.param_min : dict[str, float] = {}
        self.param_max : dict[str, float] = {}
        self.normalize = normalize
        self.move_type = move_type
        self.end_exploration = end_exploration

        self.initial_step_size = initial_step_size  # This is the initial step size for the normal move type
        self.final_step_size = final_step_size  # This is the final step size for the normal move type, after cooling down

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

    def update_minmax(self):
        for p in self.parameter_names:
            # Get the minimum and maximum values of the parameter
            min_val = min(particle[p] for particle in self.all_particles)
            max_val = max(particle[p] for particle in self.all_particles)
            # Update the prior with the new minimum and maximum values
            self.param_min[p] = min_val
            self.param_max[p] = max_val

    def indexed_simulator(self, index_candidate: Tuple[int, candidateType]) -> Tuple[int, simResultType]:
        """
        This method is used to evaluate the fitness of a candidate. It is used in the parallel evaluation of candidates.
        """
        idx, candidate = index_candidate
        max_attempts = 1
        for attempt in range(max_attempts):
            try:
                result = self.simulator(candidate)
                return idx, result
            except Exception as e:
                logging.warning(f"Simulator failed for candidate {idx} on attempt {attempt+1}: {e}")
                time.sleep(0.5)  # Optional delay
        logging.warning(f"Simulation failed after {max_attempts} attempts for candidate {idx}")
        standard_simdata = {
        'rae': np.zeros(8, dtype=np.float64),
        'ran': np.zeros(8, dtype=np.float64)
            }
        return idx, standard_simdata  # Return a default value or raise an error


    def evaluate_candidates(self, candidates: List[candidateType]) -> None: #Taken straight from evo_etc
        #TO DO: make sure new candidates are lined up correctly with the old ones. Indexed_simulater takes too long
        indexed_candidates = list(enumerate(candidates))
        results_dict = {}
        # Specifying timeout of 30 minutes
        timeout = 30*60
        # This function both evaluates newly born individuals and store them into the archive
        start = time.time()
        simulated_data = []
        # Candidates for which evaluating fitness was successful
        successfull_candidates = set()
        candidate_counter = 0
        timed_out = False
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
                        # serialized_simulator = dill.dumps(self.simulator)
                        # simulator_func = dill.loads(serialized_simulator)
                        #Old attempt
                        res_iter: Iterable[simResultType] = p.map(self.simulator, candidates,timeout=timeout).result()
                        #New attempt
                        # res_iter: Iterable[simResultType] = p.map(self.indexed_simulator, indexed_candidates,timeout=timeout).result()
                        while True:
                            try:
                                raw_res = res_iter.next() # type: ignore
                                # idx, raw_res = res_iter.next()
                            except StopIteration:
                                # We have now iterated over all particles
                                break
                            except TimeoutError:
                                timed_out = True
                                logging.info("Evaluation of particle timed out")
                            else:
                                logging.info("Evaluation of particle ran successfully")
                                # results_dict[idx] = raw_res
                                simulated_data.append(raw_res)
                                successfull_candidates.add(candidate_counter)
                            finally:
                                candidate_counter += 1
                except (OSError, RuntimeError) as e:
                    logging.error('failed parallel_evaluation_mp: {0}'.format(str(e)))
                    raise

        if timed_out:
            final_simulated_data = []
            counter = 0
            for i in range(len(candidates)):
                if i in successfull_candidates:
                    final_simulated_data.append(simulated_data[counter])
                    counter += 1
                else:
                    standard_simdata = {
            'rae': np.zeros(8, dtype=np.float64),
            'ran': np.zeros(8, dtype=np.float64)
                }
                    final_simulated_data.append(standard_simdata)

            simulated_data = final_simulated_data
            




        distances = [self.distance_function(self.Yobs, res) for res in simulated_data]
        logging.info(f"Candidate distances are {distances}")

        # save all simulated results
        self.all_simulated_data.extend(simulated_data)
        self.all_distances.extend(distances)
        # This deals with the problem of candidates failing evaluation
        # self.all_particles.extend([candidate for counter, candidate in enumerate(candidates) if counter in successfull_candidates])
        # self.all_particles.extend(candidates)
        self.all_particles.extend(candidates)

        self.birth_generation.extend(repeat(self.generation,len(simulated_data))) #Can probably be removed, replaced by birth_generation_layer

        self.birth_generation_layer.extend(repeat((self.generation,self.current_layer),len(simulated_data))) #Unused, but might be useful for later

        #The four lines above are in use. The one lines below are not. Must look further into it, maybe remove, maybe implement, maybe add more lines
        self.times_challenged.extend(repeat(0,len(simulated_data))) #This is not in use in evo_etc either, just recorded as information. Or not updated either it seems
        end = time.time()
        logging.info('Completed parallel evaluation of candidates in {0} seconds'.format(end - start))
        # logging.warning(f"Length of all_simulated_data is {len(self.all_simulated_data)}")
        # logging.warning(f"Length of all_distances is {len(self.all_distances)}")
        # logging.warning(f"Length of all_particles is {len(self.all_particles)}")




    def mutate_param(self,candidate, entry: str):
            candidate[entry] = self.priors[entry].rvfv()

    def check_validity(self, candidate, entry: str) -> bool:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt < Tm in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')
        # We assume that entries are of the form PROTID_{Tm,Topt,dCpt}
        # If this is not the case, we assume that the algorithm is used for another kind of inference problem,
        # so we skip this domain-specific check. This also applies to the dCPt as mutatating them does not violate the constraint
        # if split_entry[1] == 'dCpt':
        #     if candidate[entry] < -29675.04 or candidate[entry] > 19962.98:
        #         logging.info(f"Invalid dCpt value {candidate[entry]} for {entry}. Can work, but inspect carefully.") #Values from evo simulations
        if len(split_entry) != 2 or split_entry[1] not in ("Tm","Topt"):
            return True
        protein_id = split_entry[0]
        Tm_key = protein_id + "_Tm"
        Topt_key = protein_id + "_Topt"
        Tm = candidate[Tm_key]
        Topt = candidate[Topt_key]
        return Tm > Topt > 0

    def check_scaled_validity(self, scaled_candidate, entry: str) -> bool:
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
        Tm = scaled_candidate[Tm_key] * (self.param_max[Tm_key] - self.param_min[Tm_key]) + self.param_min[Tm_key]
        Topt = scaled_candidate[Topt_key] * (self.param_max[Topt_key] - self.param_min[Topt_key]) + self.param_min[Topt_key]

        return Tm > Topt > 0

    def adaptive_step_size(self) -> float:
        if self.generation > self.maxiter * self.end_exploration:
            return self.final_step_size
        decay_rate = np.log(self.final_step_size/ self.initial_step_size) / (self.maxiter*self.end_exploration)
        return self.initial_step_size * np.exp(decay_rate * self.generation)

    def energy_function(self, current_dist: float, candidate_dist: float) -> bool:
        delta = candidate_dist - current_dist
        exponent = -delta / self.current_temp
        exponent = min(exponent, 700)  # Prevent overflow in exp
        acceptance_probability = np.exp(exponent)
        if self.log_ef:
            logging.debug(f"Energy function value is : {np.exp(-(delta)/self.current_temp)}")
            self.log_ef_list.append(np.exp(-(delta)/self.current_temp))
            self.log_ef = False
        if candidate_dist < current_dist:
            return True
        else:
            return self.rng.random() < acceptance_probability



    def choose_particle(self, idx1: int, idx2: int)-> int:
        d1 = self.all_distances[idx1]
        d2 = self.all_distances[idx2]
        if self.energy_function(d1, d2):
            return idx2
        else:
            return idx1



    def replace_population(self,original_population: npt.NDArray[np.int64], candidates: npt.NDArray[np.int64]):
        logging.info(f"Replacing population with candidates")
        current_population : List[int] = []
        acceptance_counter: int = 0
        for i, (old_particle, candidate_particle) in enumerate(zip(original_population, candidates)):
            chosen_particle = self.choose_particle(old_particle, candidate_particle)
            if chosen_particle == old_particle:
                current_population.append(old_particle)
                #Keep index in current population as it is
                self.inner_iterations_list[i] += 1 #Add 1 to inner_iterations_list
                logging.debug(f"Keeping particle {i}") #For checking, remove later

            else:
                current_population.append(candidate_particle)
                #Add index of this particle to self.population
                acceptance_counter += 1
                self.inner_iterations_list[i] = self.min_layers #Set inner_iterations_list to 1
                #self.is_improved_list[i] = True
                logging.debug(f"Replacing particle {i}") #For checking, remove later
        acceptance_rate: float = acceptance_counter / len(candidates)
        logging.info(f"Acceptance rate in generation {self.generation}: {acceptance_rate}")
        self.acceptance_rates.append(acceptance_rate)
        self.population.append(list(current_population))

    def get_bounds(self, param: str) -> Tuple[float, float]:
        """
        This method returns the bounds of the parameters. It is used to normalize the parameters
        """
        if param.endswith('_Tm'):
            min_val = 0
            max_val = 400
        elif param.endswith('_Topt'):
            min_val = 0
            max_val = 400
        elif param.endswith('_dCpt'):
            min_val = -12000
            max_val = -4000
        else:
            min_val = -10
            max_val = 10
        return min_val, max_val

    def normalize_particle(self, particle_idx: np.int64) -> candidateType:
        min_val = -10
        max_val = 10

        #scaled_candidate = {parameter: (value - self.param_min[parameter])/(self.param_max[parameter]- self.param_min[parameter]) for parameter, value in self.all_particles[particle_idx].items()} #Max/min from current parameters
        scaled_candidate = {
        parameter: (value - self.get_bounds(parameter)[0]) / (self.get_bounds(parameter)[1] - self.get_bounds(parameter)[0])
        for parameter, value in self.all_particles[particle_idx].items()
        }
        return scaled_candidate

    def denormalize_particle(self, scaled_candidate: candidateType) -> candidateType:
        min_val = -10
        max_val = 10
        #denormalized_candidate = {parameter: value * (self.param_max[parameter]- self.param_min[parameter]) + self.param_min[parameter] for parameter, value in scaled_candidate.items()} #Max/min from current parameters
        denormalized_candidate = {
        parameter: value * (self.get_bounds(parameter)[1] - self.get_bounds(parameter)[0]) + self.get_bounds(parameter)[0]
        for parameter, value in scaled_candidate.items()
        }
        return denormalized_candidate


    def generate_candidates(self, particle_idxs: npt.NDArray[np.int64]) -> None: #Kan bruke change_all_parameters her istedenfor å skrive det eksplisitt
        logging.info("Generating candidates")
        candidates: List[candidateType] = []
        if self.normalize:
            for idx in particle_idxs:
                candidate = {parameter: value for parameter, value in self.all_particles[idx].items()}
                for key in candidate:
                    old_parameter_value = candidate[key]
                    step_size = self.adaptive_step_size()
                    if key.endswith('_dCpt'):
                        rel_step_size = step_size/30
                    else:
                        rel_step_size = step_size/old_parameter_value
                    candidate[key] *= 1 + rel_step_size * (self.rng.random()-0.5) #Endre til å bruke change_all_parameters, og/eller måte på å endre verdiene
                    if not self.check_validity(candidate, key):
                        candidate[key] = old_parameter_value
                candidates.append(candidate)
        else:
            for idx in particle_idxs:
                candidate = {parameter: value for parameter, value in self.all_particles[idx].items()}
                for key in candidate:
                    old_parameter_value = candidate[key]
                    if self.end_exploration < 1.0:
                        #logging.info("Using adaptive step size")
                        step_size = self.adaptive_step_size()
                        if key.endswith('_dCpt'):
                            candidate[key] += step_size * (self.rng.random() - 0.5)
                        elif key.endswith('_Tm') or key.endswith('_Topt'):
                            candidate[key] += step_size * (self.rng.random() - 0.5)
                        else:
                            candidate[key] += step_size * (self.rng.random() - 0.5)
                    else:
                        if key.endswith('_dCpt'):
                            candidate[key] += self.dCpt_step_size * (self.rng.random() - 0.5)
                        elif key.endswith('_Tm') or key.endswith('_Topt'):
                            candidate[key] += self.temp_step_size * (self.rng.random() - 0.5)
                        else:
                            candidate[key] += self.step_size * (self.rng.random() - 0.5)
                            candidate[key] = np.clip(candidate[key], -5.12, 5.12)  # Clip to avoid extreme values
                        if not self.check_validity(candidate, key):
                            candidate[key] = old_parameter_value
                candidates.append(candidate)
        logging.info("Evaluating fitness of candidates")
        self.evaluate_candidates(candidates)

    def simulate_generation(self):
        current_population = np.array(list(self.population[-1]))
        self.generate_candidates(current_population)
        candidates_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
        self.replace_population(current_population, candidates_idxs)

        logging.info(f"Current population is {self.population[-1]}")
        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        min_generation_epsilon = min(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        #self.update_std() Takes a long time, and not used
        #self.update_minmax()
        logging.info(f"Model epsilon {max_generation_epsilon}")
        logging.info(f"Model min epsilon {min_generation_epsilon}")



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
            #self.update_std() takes a long time, and not used
            self.update_minmax()
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
            self.simulate_generation()

            self.generation += 1
            if self.current_temp > self.final_temp:
                self.current_temp *= self.cooling_rate

            # #Adaptive cooling
            # recent_rate = np.mean(self.acceptance_rates[-5:])
            # if recent_rate > upper_threshold:
            #     # Too many moves accepted → cool down faster (exploit)
            #     self.current_temp = max(self.current_temp * self.cooling_rate, self.final_temp)
            # elif recent_rate < lower_threshold:
            #     # Too few moves accepted → reheat slightly (explore)
            #     self.current_temp = min(self.current_temp / self.cooling_rate, self.initial_temp)

            logging.info(f" Inner_iterations_list: {self.inner_iterations_list}")
            if self.save_intermediate and self.generation % 250 == 0:
                logging.info(f"Saving intermediate results to {self.outfile}")
                dill.dump(self,open(self.outfile,'wb'))
            
            end_generation = time.time()
            logging.info(f"Generation {self.generation} completed in {end_generation-start_generation} seconds")

        else:
            # This else-clause belongs to the main simulated annealing loop
            logging.info("Fitness objective not reached after maximum number of generations")
            logging.info("Exiting simulated annealing")
            logging.info(f"Self.log_ef_list: {self.log_ef_list}")
            logging.info(f"Self.acceptance_rates: {self.acceptance_rates}")
            logging.info(f"Final temperature: {self.current_temp}")

        logging.info(f"Saving results to {self.outfile}")
        dill.dump(self, file=open(self.outfile,mode='wb'))



