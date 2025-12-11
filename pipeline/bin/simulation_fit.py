#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import sys
import json


#A. Raue et al., PLOS ONE 8 (2013) e74335: Lessons Learned from Quantitative Dynamical Modeling in Systems Biology - models used in the article

def epo_receptor(t,y,theta):
    """
ODE system for EPO-receptor dynamics.

Parameters
----------
t : float
    Current time (minutes)
y : list or np.ndarray
    Current state vector [epo, epor, epo_epor, epo_epor_i, depo_i, depo_e]
theta : list or np.ndarray
    6 kinetic parameters [k_on, k_t, k_e, k_ex, k_di, k_de] in linear space

Returns
-------
dydt : list
    Time derivatives of the state vector
"""
    y_clamped = np.maximum(y,0) # to avoid negative concentrations
    epo,epor,epo_epor,epo_epor_i,depo_i,depo_e = y_clamped
    k_on,k_t,k_e,k_ex,k_di,k_de=theta #logarithmized parameter space : phi goes from -3 to 3 --> theta goes from 10^-3 to 10^3  --> theta_i=10^{phi_i}

    kd=164
    k_off=kd*k_on
    b_max=129
    
    #rate equations
    v1=k_on*epo*epor
    v2=k_off*epo_epor
    v3=k_t*b_max #b_max=init_epor
    v4=k_t*epor
    v5=k_e*epo_epor
    v6=k_ex*epo_epor_i
    v7=k_di*epo_epor_i
    v8=k_de*epo_epor_i

    #ODEs
    ydot1=-v1+v2+v6
    ydot2=-v1+v2+v3-v4+v6
    ydot3=v1-v2-v5
    ydot4=v5-v6-v7-v8
    ydot5=v7
    ydot6=v8

    return [ydot1, ydot2, ydot3, ydot4, ydot5, ydot6]

def objective(phi, ym, tm):
    """
    Compute negative log-likelihood for given parameters.
    
    phi: list or array of 9 parameters [k_on,k_t,k_e,k_ex,k_di,k_de, init_epo, sa, sb]
    ym: experimental measurements (timepoints x observables)
    tm: time points array
    """
    try:
        theta = 10**np.array(phi)
        init_epo = theta[6]
        init_epor = 129  # b_max
        sa = theta[7]
        sb = theta[8]

        y0 = [init_epo, init_epor, 0, 0, 0, 0]
        t_span = (0, tm[-1])
        sol = solve_ivp(
            epo_receptor,
            t_span,
            y0,
            method='BDF',
            t_eval=tm,
            args=(theta[:6],)
        )
        states = sol.y.T
        epo_medium = states[:,0] + states[:,5]
        epo_surf = states[:,2]
        epo_cells = states[:,3] + states[:,4]
        obs = np.column_stack([epo_medium, epo_surf, epo_cells])

        eps = 1e-6
        like = np.empty(obs.shape[1])
        log_ym = np.log10(np.maximum(ym, eps))
        log_obs = np.log10(np.maximum(obs, eps))
        var = (sa + sb*log_obs)**2
        var = np.maximum(var, eps)

        for k in range(obs.shape[1]):
            like[k] = np.sum(
                np.log(2*np.pi*var[1:,k]) +
                ((log_ym[1:,k]-log_obs[1:,k])**2 / var[1:,k])
            )

        likelihood = np.sum(like) + np.log(2*np.pi*var[0,0]) + ((log_ym[0,0]-log_obs[0,0])**2 / var[0,0])
        return likelihood
    except Exception:
        return np.inf

    # Compute likelihood
    likelihood = objective(phi, ym, tm)

 #MAIN: Nextflow integration
if __name__ == "__main__":
    phi_json = sys.argv[1]     # From Nextflow: "[-2.278,-2.642,...]"
    dataset_file = sys.argv[2] # time_series.csv
    
    # Parse phi0 (your 9 parameters)
    with open(phi_json, 'r') as f:
        phi = json.load(f)
    
    # Load dataset CSV
    df = pd.read_csv(dataset_file, index_col=0).transpose()
    tm = df['time_points'].values
    ym = df[['epo_medium', 'epo_surf', 'epo_cells']].values
    
    # Run likelihood computation
    likelihood = objective(phi, ym, tm)
    
    # Save result
    pd.DataFrame({
        'phi0': [phi],
        'likelihood': [likelihood]
    }).to_csv('fit_result.csv', index=False, header=False)