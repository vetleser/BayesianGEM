#!/usr/bin/env python
# coding: utf-8

# In[1]:


import evo_etc as evo
import numpy as np
import GEMS
import os
import pandas as pd
import logging
import pickle

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

# In[2]:


def main():
    task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

    Yobs_batch = GEMS.aerobic_exp_data()
    Yobs_chemo = GEMS.chemostat_exp_data()
    #Yobs_batch_an = GEMS.anaerobic_exp_data()
    dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
    sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
    Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}


    # In[5]:


    Yobs = {'rae':Yobs_batch['data'],
            'chemostat':Yobs_chemo['data'],
            'ran':Yobs_batch_an['data']}


    # In[ ]:


    path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
    params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)


    # #### Define priors

    # In[ ]:




                    
    Yobs = {'rae':Yobs_batch['data'],
            'chemostat':Yobs_chemo['data'],
            'ran':Yobs_batch_an['data']}

    outdir = '../results/permuted_smcevo_res'
    candidate_frame: pd.DataFrame = pickle.load(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='rb'))
    entry = candidate_frame.iloc[task_idx]



    # In[]

    min_epsilon = -1.0 # equivalent to r2 score of 1
    population_size = 100
    logging.info('Initialize model')


    priors, simulation,prior_name, outfile, random_seed = entry[["priors", "simulation", "prior_name", "outfile","random_seed"]]


    # #### Define model settings

    # In[ ]:
    rng = np.random.default_rng(random_seed)
    min_epsilon = -1.0 # equivalent to r2 score of 1
    population_size = 128
    outfile = '../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl'


    # In[ ]:


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
                        maxiter=500,
                        rng=rng,
                        mutation_frequency=100,
                        mutation_prob=.5,
                        selection_proportion=0.5,
                        n_children=1,
                        save_intermediate=True)

    # #### Run simulation

    # In[ ]:


    logging.info(f'Start evolutionary simulations with prior parameter set {prior_name}, simulation {simulation}')
    model.run_simulation()
    logging.info("DONE")


if __name__ == '__main__':
    main()
