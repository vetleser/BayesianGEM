#!/usr/bin/env python
# coding: utf-8

# In[1]:


from random_sampler import RV
import numpy as np
import GEMS
import os
import pandas as pd
import pickle
import logging

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

# In[ ]:

num_elites = [16,32,64,128]
n_replicates = 2
overall_random_seed = 158
rng = np.random.default_rng(overall_random_seed)
path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)

candidate_frame: pd.DataFrame = pd.DataFrame(index=pd.MultiIndex.from_product([num_elites, range(n_replicates)], names = ["num_elites", "simulation"])).reset_index()




# In[ ]:


min_epsilon = -1.0 # equivalent to r2 score of 1
population_size = 100
outdir = '../results/evo_truncation'
if not os.path.exists(outdir):
    os.makedirs(outdir)
candidate_frame['outfile'] = [f'{outdir}/smcevo_gem_truncation_{num_elites}_{simulation}_save_all_particles.pkl' for
 num_elites, simulation in zip(candidate_frame['num_elites'],candidate_frame['simulation'])]
candidate_frame['random_seed'] = rng.choice(range(0,100000), candidate_frame.shape[0])

pickle.dump(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='wb'),obj=candidate_frame)
