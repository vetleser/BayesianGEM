#!/usr/bin/env python
# coding: utf-8

# In[1]:


import evo_etc as evo
import numpy as np
import GEMS
import os
import pandas as pd
import logging
import dill

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

# In[2]:


def main():

    task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    outdir = '../results/permuted_smcevo_res'
    candidate_frame: pd.DataFrame = pickle.load(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='rb'))
    entry = candidate_frame.iloc[task_idx]
    priors, simulation,prior_name, outfile, random_seed = entry[["priors", "simulation", "prior_name", "outfile","random_seed"]]
    maxiter = 2000

    def model_init():
        Yobs_batch = GEMS.aerobic_exp_data()
        Yobs_chemo = GEMS.chemostat_exp_data()
        #Yobs_batch_an = GEMS.anaerobic_exp_data()
        dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
        sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
        Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}

        Yobs = {'rae':Yobs_batch['data'],
                'chemostat':Yobs_chemo['data'],
                'ran':Yobs_batch_an['data']}

                    
        Yobs = {'rae':Yobs_batch['data'],
                'chemostat':Yobs_chemo['data'],
                'ran':Yobs_batch_an['data']}


        # #### Define model settings

        # In[ ]:
        rng = np.random.default_rng(random_seed)
        min_epsilon = -1.0 # equivalent to r2 score of 1
        population_size = 128
        

        logging.info('Initialize model')
        # Parameters are set to generate 128 children per generation in order to provide 
        # comparable results with the Bayesian fitting algorithm
        model = evo.GA(simulator= GEMS.simulate_at_three_conditions_2,
                                priors=priors,
                                min_epsilon=min_epsilon,
                                generation_size=population_size,
                                distance_function=GEMS.distance_2,
                                Yobs=Yobs,
                                outfile=outfile,
                                maxiter=maxiter,
                                rng=rng,
                                mutation_frequency=100,
                                mutation_prob=.5,
                                selection_proportion=0.5,
                                n_children=1,
                                save_intermediate=False)
        
        return model

    if os.path.exists(outfile):
        model: evo.GA = dill.load(open(outfile,"rb"))
        model.maxiter = maxiter
        logging.info(f'Resuming evolutionary simulations with prior parameter set {prior_name}, simulation {simulation}')
    else:
        model = model_init()
        logging.info(f'Start evolutionary simulations with prior parameter set {prior_name}, simulation {simulation}')

    model.run_simulation()
    logging.info("DONE")


if __name__ == '__main__':
    main()
