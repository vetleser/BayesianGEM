#!/usr/bin/env python
# coding: utf-8

import dill
import os
import numpy as np


task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

rng = np.random.default_rng(task_idx)

dill.dump(rng,file=open(f"../results/debug_results_{task_idx}.pkl",'wb'))
