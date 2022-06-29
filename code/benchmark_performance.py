#!/usr/bin/env python
# coding: utf-8
# The purpose of this script is to benchmark etcFBA model evaluation across the COBRApy and
# reframed implementations. Note that due to history constrains, this script cannot be used to benchmark
# both frameworks at the same time. Therefore, some git trickery is needed, consult benchmark_performance.sh
# for this challenge

# In[]
import time
import GEMS
import os
import pandas as pd
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
N_ITERATIONS = 1

path = '../'
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
particle = dict()
for ind in params.index: 
    for col in ['Tm','Topt','dCpt']: 
        particle['{0}_{1}'.format(ind,col)] = params.loc[ind,col]

# In[]
tic = time.perf_counter()
for _ in range(N_ITERATIONS):
    total_res = GEMS.simulate_at_three_conditions_2(particle)
toc = time.perf_counter()
print(f"Computed results at three conditions with an average of {(toc-tic)/N_ITERATIONS} seconds.")
