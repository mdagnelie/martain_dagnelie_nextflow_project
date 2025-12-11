#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy.stats import qmc
import sys
import json

# Settings

n_samples = int(sys.argv[1])  # number of initial guesses
n_dims = 9           # number of parameters
ym0 = 1990           # experimental value used for init_epo
np.random.seed(42)

# Latin Hypercube Sampling 
sampler = qmc.LatinHypercube(d=n_dims)
sample = sampler.random(n=n_samples)  # shape: (n_samples, n_dims)
ite, dim = sample.shape

param_list=[]
# Generate parameter sets lying in the bounds
for idx in range(ite):
    phi0 = np.empty(dim)
    
    # Parameters 0-5: [-7, 3] in log10 space
    phi0[:6] = sample[idx, :6]*10 - 7
    
    # Parameter 6: init_epo ~ 0.9–1.1 * ym0
    phi0[6] = np.log10(ym0 * (0.2*sample[idx, 6] + 0.9))
    
    # Parameters 7-8: sa, sb ~ [-3, 3]
    phi0[7:] = sample[idx, 7:]*6 - 3
    
    with open(f"{idx}_phi.json", "w") as f:
        json.dump(phi0.tolist(), f)

