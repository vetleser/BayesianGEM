#!/usr/bin/env python
# coding: utf-8

# In[1]:


import evo_etc as evo
import sa_etc as sa
import numpy as np
import GEMS
import os
import pandas as pd
import logging
import pickle
from random_sampler import RV

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

# In[2]:


def main():
    task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    outdir = '../results/sa'
    candidate_frame: pd.DataFrame = pickle.load(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='rb'))
    entry = candidate_frame.iloc[task_idx]
    simulation, outfile, random_seed, final_temp, move_type, step_size = entry[["simulation", "outfile","random_seed",
    "final_temp", "move_type", "step_size"]]
    maxiter = 1000
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
    population_size = 128
    n_children = 128
    

    logging.info('Initialize model')
    model = sa.SimulatedAnnealing(simulator= GEMS.simulate_at_two_conditions_2,
                             priors=priors,
                             Yobs=Yobs,
                             maxiter=maxiter,
                             generation_size=population_size,
                             min_epsilon=min_epsilon,
                             outfile=outfile,
                             rng=rng,
                             distance_function=GEMS.distance_2,
                             version = 2,
                             normalize=False,
                             final_temp=final_temp,
                             save_intermediate=True,
                             step_size=step_size,
                             move_type=move_type)
    
    
    logging.info(f"""Start evolutionary simulations with Simulated Annealing,
     final_temp: {final_temp}, normalize: False, simulation: {simulation}, step_size: {step_size}, move_type: {move_type}""")

    model.run_simulation()
    logging.info("DONE")


if __name__ == '__main__':
    main()
