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
from random_sampler import RV
from multiprocessing import cpu_count

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

# In[2]:


def main():
    task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    outdir = '../results/crowdingDE'
    candidate_frame: pd.DataFrame = pickle.load(file=open(file=f'{outdir}/simulation_skeleton_test.pkl',mode='rb'))
    entry = candidate_frame.iloc[task_idx]
    simulation, outfile, random_seed, scaling_factor, crossover_prob = entry[["simulation", "outfile","random_seed",
    "scaling_factor","crossover_prob"]]
    maxiter = 10
    Yobs_batch = GEMS.aerobic_exp_data()
    #Yobs_batch_an = GEMS.anaerobic_exp_data()
    dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
    sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
    Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}

    Yobs = {'rae':Yobs_batch['data'],
            'ran':Yobs_batch_an['data']}

    path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
    params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
    priors = dict()
    for ind in params.index: 
        for col in ['Tm','Topt','dCpt']: 
            priors['{0}_{1}'.format(ind,col)] = RV('normal',
                                                        loc=params.loc[ind,col],
                                                        scale=params.loc[ind,col+'_std'])


    # #### Define model settings

    # In[ ]:
    rng = np.random.default_rng(random_seed)
    min_epsilon = -1.0 # equivalent to r2 score of 1
    population_size = 256
    cores = int(os.getenv("SLURM_CPUS_PER_TASK"))
    

    logging.info('Initialize model')
    model = evo.CrowdingDE(simulator= GEMS.simulate_at_two_conditions_2,
                            priors=priors,
                            min_epsilon=min_epsilon,
                            generation_size=population_size,
                            distance_function=GEMS.distance_2,
                            Yobs=Yobs,
                            outfile=outfile,
                            maxiter=maxiter,
                            rng=rng,
                            scaling_factor=scaling_factor,
                            crossover_prob=crossover_prob,
                            n_children=128,
                            save_intermediate=False,
                            cores=cores
                            )
    
    
    logging.info(f"""Start evolutionary simulations with CrowdingDE,
     scaling factor {scaling_factor}, crossover probability {crossover_prob}, simulation {simulation}, population size {population_size}, cores: {cores}""")

    model.run_simulation()
    logging.info("DONE")


if __name__ == '__main__':
    main()
