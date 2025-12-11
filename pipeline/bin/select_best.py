#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sys, json

def epo_receptor_scipy(t,y,theta):
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

if __name__ == "__main__":
    all_fits = sys.argv[1]      # "all_fits.csv"
    dataset_file = sys.argv[2]  # "all_timeseries.csv"
    
    # Load your CSV format perfectly
    df = pd.read_csv(all_fits,header=None, names=['phi0','likelihood'])
    obj_opt = df['likelihood'].values
    
    # 1. Objective plot 
    plt.figure(figsize=(8,6))
    obj_order = np.sort(obj_opt)
    plt.plot(np.arange(0,len(obj_opt))+1,obj_order,'-o',ms=3,mfc='w')
    plt.yscale("log")
    plt.xlabel('iteration'); plt.ylabel('Neg Log-Likelihood')
    plt.title('Objective Distribution')
    plt.savefig('objective_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. best parameters
    best_idx = np.argmin(obj_opt)
    phi_str = df.iloc[best_idx]['phi0']  # JSON string
    phi_best = json.loads(phi_str)       # Parse to list
    
    with open('best_phi.json', 'w') as f:
        json.dump(phi_best, f)
    
    # 3. Best fit plot 
    df2 = pd.read_csv(dataset_file, index_col=0).transpose()
    tm = df2['time_points'].values                  
    ym = df2[['epo_medium', 'epo_surf', 'epo_cells']].values.transpose()
    
    theta = 10**np.array(phi_best)
    init_epo, init_epor, sa, sb = theta[6], 129, theta[7], theta[8]
    init_epo= theta[6]
    init_epor= 129
    sa=theta[7]
    sb=theta[8]
    t_end=tm[-1]
    y0=[init_epo,init_epor,0,0,0,0] #epo,epor,epo_epor,epo_epor_i,depo_i,depo_e
    t_span=(0,t_end) #minutes
    sol = solve_ivp(epo_receptor_scipy, t_span, y0, method='BDF', args=(theta[:6],))
    states=sol.y.T #epo,epor,epo_epor,epo_epor_i,depo_i,depo_e
    epo_medium=states[:,0]+states[:,5] #epo+depo_e
    epo_surf=states[:,2] #epo_epor
    epo_cells=states[:,3]+states[:,4] #epo_epor_i+depo_i
    obs=np.array([epo_medium,epo_surf,epo_cells]).T 
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12, 3))


    ax1.scatter(tm,ym[0],label="Epo in medium",color='k',s=10)
    ax1.plot(sol.t,obs[:,0],color='red')

    ax2.scatter(tm[1:],ym[1][1:],color='k',label='Epo on surface',s=10)
    ax2.plot(sol.t,obs[:,1],color='blue')

    ax3.scatter(tm[1:],ym[2][1:],color='k',label='Epo in cells',s=10)
    ax3.plot(sol.t,obs[:,2],color='green')

    ax1.set_ylim([-50,2500])
    ax2.set_ylim([-5,350])
    ax3.set_ylim([-50,800])

    ax1.set_xlabel('Time [min]')
    ax2.set_xlabel('Time [min]')
    ax3.set_xlabel('Time [min]')
    ax1.set_ylabel('[Epo] (pM)')

    plt.savefig('best_fit.png', dpi=300, bbox_inches='tight')
    plt.close()