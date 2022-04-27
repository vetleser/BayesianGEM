#!/usr/bin/env python
# coding: utf-8

# This is a companion script to gem_smcabc_at_three_conditions_prepare.py
# Its rationale is to devide the jobs to the aforementioned script on serveral nodes
# such that several nodes can work on the computations at once.
# This process is so computationally expensive that it is difficult for a single node
# to generate results in a reasonable amount of time
#In[]
import abc_etc as abc
import numpy as np
import GEMS
import os
import pandas as pd
import pickle
import logging

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])


# In[2]:


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

outdir = '../results/permuted_smcabc_res'
candidate_frame: pd.DataFrame = pickle.load(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='rb'))
entry = candidate_frame.iloc[task_idx]



# In[]

min_epsilon = -1.0 # equivalent to r2 score of 1
population_size = 100
logging.info('Initialize model')


priors, origin, status, outfile, random_seed = entry[["priors", "origin", "status", "outfile","random_seed"]]
model = abc.SMCABC(GEMS.simulate_at_three_conditions_2,
                                priors=priors,
                                min_epsilon=min_epsilon,
                                population_size=population_size,
                                distance_function=GEMS.distance_2,
                                Yobs=Yobs,
                                outfile=outfile,
                                generation_size=128,
                                 maxiter=500)




# In[ ]:
# #### Run simulation

np.random.seed(random_seed)
logging.info(f'Start simulations with prior parameter set {origin}, {status}')
model.run_simulation()

logging.info(f'DONE')
