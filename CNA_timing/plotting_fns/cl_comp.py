import pickle
import arviz as az
import numpy as np
import pandas as pd
from scipy.special import softmax
import corner
import matplotlib.pyplot as plt
import sys
import os
from itertools import chain
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from plot import plot_mu_uncertainty, plot_event_time_uncertainty, plot_ppc, plot_pair, plot_distributions


sim_loh_az = '/Volumes/files/az/pool/loh/mu_regime/norm(0,0.01)/0.0075.nc'
sim_loh_sum = '/Volumes/files/az/pool/loh/mu_regime/norm(0,0.01)/0.0075.pkl'
hier_cnloh_az = '/Volumes/files/az/hier/comp/loh_comp.nc'
hier_cnloh_sum = '/Volumes/files/az/hier/comp/loh_comp.pkl'
patient_idx = np.arange(1, 31)  

with open(sim_loh_sum, 'rb') as f:
    sim_sum = pickle.load(f)
sim_az = az.from_netcdf(sim_loh_az)

with open(hier_cnloh_sum, 'rb') as f:
    hier_sum = pickle.load(f)
hier_az = az.from_netcdf(hier_cnloh_az)

event_times =  [18, 12, 42, 17, 5, 47, 30, 17, 19, 9, 9, 30, 21, 40, 15, 9, 10, 9, 24, 7, 47, 18, 13, 36, 15, 50, 57, 19, 7, 32]

t_sim = np.array(sim_sum[sim_sum.index.str.contains(r'^t\[\d+\]$')]['Mean'].values)
t_hier = np.array(hier_sum[hier_sum.index.str.contains(r'^t\[\d+\]$')]['Mean'].values)

sim_ci = az.hdi(sim_az, var_names=["t"], hdi_prob=0.95) 
hier_ci = az.hdi(hier_az, var_names=["t"], hdi_prob=0.95)  

sim_ci_size = sim_ci["t"].sel(hdi="higher") - sim_ci["t"].sel(hdi="lower")
hier_ci_size = hier_ci["t"].sel(hdi="higher") - hier_ci["t"].sel(hdi="lower")
mean_sim = np.mean(np.array(sim_ci_size))
mean_hier = np.mean(np.array(hier_ci_size))


patient_idx = np.arange(1, len(sim_ci_size) + 1)

plt.figure(figsize=(7, 5))
plt.scatter(patient_idx, sim_ci_size, color="blue", alpha=1, label=f"Pooled (mean={mean_sim:.2f})")
plt.scatter(patient_idx, hier_ci_size, color="red", alpha=1, label=f"Hierarchical (mean={mean_hier:.2f})")

plt.vlines(patient_idx, 0, sim_ci_size, color="blue", alpha=0.6, linewidth=1)
plt.vlines(patient_idx, 0, hier_ci_size, color="red", alpha=0.6, linewidth=1)

plt.xlabel("Patient Index", fontsize=14.5)
plt.ylabel("CI Size (95% HDI)", fontsize=14.5)

plt.xticks(np.arange(0, 31, 5), fontsize=12)
plt.xlim(0, 31)

plt.yticks(fontsize=12)
plt.ylim(0, 50)

plt.gca().tick_params(axis='y', which='both', left=True, right=False)
plt.gca().tick_params(axis='x', which='both', bottom=True, top=False)

plt.legend(loc="upper right", fontsize=11.5, bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("loh_ci_comp.pdf", format="pdf", dpi=300)
plt.show()
