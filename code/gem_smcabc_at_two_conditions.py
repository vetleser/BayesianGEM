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
import abc_etc as abc

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# In[2]:


def main():
    task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    outdir = '../results/reduced_smcabc_res'
    candidate_frame: pd.DataFrame = pickle.load(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='rb'))
    entry = candidate_frame.iloc[task_idx]
    simulation, outfile, random_seed = entry[["simulation", "outfile","random_seed"]]
    maxiter = 500
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
    population_size = 100
    

    logging.info('Initialize model')
    # Parameters are set to generate 128 children per generation in order to provide 
    # comparable results with the Bayesian fitting algorithm
    model = abc.SMCABC(GEMS.simulate_at_two_conditions_2,
                                priors=priors,
                                min_epsilon=min_epsilon,
                                population_size=population_size,
                                distance_function=GEMS.distance_2,
                                Yobs=Yobs,
                                outfile=outfile,
                                generation_size=128,
                                maxiter=500)

    
    np.random.seed(random_seed)
    logging.info(f'Start simulation {simulation+1} with at two conditions')
    model.run_simulation()

    logging.info(f'DONE')


if __name__ == '__main__':
    main()