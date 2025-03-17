import pickle
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# sim_loh_az = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)/0.0075.nc'
# sim_loh_sum = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)/0.0075.pkl'
# hier_cnloh_az = '/Volumes/files/az/hier/comp/tet_comp.nc'
# hier_cnloh_sum = '/Volumes/files/az/hier/comp/tet_comp.pkl'
sim_loh_az = 'k4/30.nc'
sim_loh_sum = 'k4/30.pkl'
hier_cnloh_az = 'k4/20.nc'
hier_cnloh_sum = 'k4/20.pkl'


patient_idx = np.arange(1, 31) 

with open(sim_loh_sum, 'rb') as f:
    sim_sum = pickle.load(f)
sim_az = az.from_netcdf(sim_loh_az)

with open(hier_cnloh_sum, 'rb') as f:
    hier_sum = pickle.load(f)
hier_az = az.from_netcdf(hier_cnloh_az)

event_times = np.array([39,8,48,5,29,54,9,8,20,5,5,29,5,57,5,64,12,9,41,4,37,23,15,37,27,56,10,5,28,47])

t_sim = np.array(sim_sum[sim_sum.index.str.contains(r'^t\[\d+\]$')]['Mean'].values)
t_hier = np.array(hier_sum[hier_sum.index.str.contains(r'^t\[\d+\]$')]['Mean'].values)

sim_res = t_sim - event_times
hier_res = t_hier - event_times

me_sim = np.mean(sim_res)  
me_hier = np.mean(hier_res)  

plt.figure(figsize=(7, 5))
plt.axhline(0, color='gray', linestyle='--', linewidth=1)  

plt.scatter(patient_idx, sim_res, color='blue', alpha=1, 
            label=f'Pooled (ME={me_sim:.2f})')
plt.scatter(patient_idx, hier_res, color='red', alpha=1, 
            label=f'Hierarchical (ME={me_hier:.2f})')

plt.vlines(patient_idx, 0, sim_res, color='blue', alpha=1, linewidth=1)
plt.vlines(patient_idx, 0, hier_res, color='red', alpha=1, linewidth=1)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("Patient index", fontsize=14.5)
plt.ylabel(r"Posterior mean - true event time", fontsize=14.5)

plt.legend(loc='upper right', fontsize=11.7)
plt.ylim(-18, 18)
plt.tight_layout()
plt.savefig('res_tet.pdf', format='pdf', dpi=300)
plt.show()