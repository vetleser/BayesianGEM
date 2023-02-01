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

scaling_factor = [0.5,1]
crossover_prob = [0.9,0.99,0.999]
n_replicates = 2
overall_random_seed = 158
rng = np.random.default_rng(overall_random_seed)
path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)

candidate_frame: pd.DataFrame = pd.DataFrame(index=pd.MultiIndex.from_product([scaling_factor, crossover_prob, range(n_replicates)],
 names = ["scaling_factor","crossover_prob", "simulation"])).reset_index()




# In[ ]:


min_epsilon = -1.0 # equivalent to r2 score of 1
population_size = 100
outdir = '../results/crowdingDE'
if not os.path.exists(outdir):
    os.makedirs(outdir)
candidate_frame['outfile'] = [f'{outdir}/smcevo_gem_{scaling_factor}_{crossover_prob}_{simulation}.pkl' for
 scaling_factor, crossover_prob, simulation in zip(candidate_frame['scaling_factor'],candidate_frame['crossover_prob'],candidate_frame['simulation'])]
candidate_frame['random_seed'] = rng.choice(range(0,100000), candidate_frame.shape[0])

pickle.dump(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='wb'),obj=candidate_frame)
