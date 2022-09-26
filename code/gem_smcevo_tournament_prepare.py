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

localities = [1,64,128,256]
n_replicates = 2
overall_random_seed = 158
rng = np.random.default_rng(overall_random_seed)
path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)


candidate_frame: pd.DataFrame = pd.DataFrame(index=pd.MultiIndex.from_product([localities, range(n_replicates)], names = ["locality", "simulation"])).reset_index()




# In[ ]:

outdir = '../results/evo_tournament'
if not os.path.exists(outdir):
    os.makedirs(outdir)
candidate_frame['outfile'] = [f'{outdir}/smcevo_gem_tournament_{locality}_{simulation}.pkl' for
 locality, simulation in zip(candidate_frame['locality'],candidate_frame['simulation'])]
candidate_frame['random_seed'] = rng.choice(range(0,100000), candidate_frame.shape[0])

pickle.dump(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='wb'),obj=candidate_frame)
