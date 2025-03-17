import pickle
import pandas as pd
import arviz as az 
import matplotlib.pyplot as plt
import numpy as np  
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from plot import convergence_check, plot_event_time_uncertainty, plot_ppc, posterior_single_mu_gamma

type = 1
cnloh_az = '/Volumes/files/az/pool/loh/mu_regime/norm(0,0.01)/0.025.nc'
cnloh_sum = '/Volumes/files/az/pool/loh/mu_regime/norm(0,0.01)/0.025.pkl'
tri_az = 'az/pool/tri/0.005/0.5/500.nc'
tri_sum = 'az/pool/tri/0.005/0.5/500.pkl'
tet_az =  'az/pool/tet/0.01/0.5/300.nc'
tet_sum = 'az/pool/tet/0.01/0.5/300.pkl'
      

if type == 1:
    with open(cnloh_sum, 'rb') as f:
        summary_df = pickle.load(f)

    az_fit = az.from_netcdf(cnloh_az)

if type == 2:
    az_fit = az.from_netcdf(tri_az)

    with open(tri_sum, 'rb') as f:
        summary_df = pickle.load(f)
if type == 3:
    az_fit = az.from_netcdf(tet_az)

    with open(tet_sum, 'rb') as f:
        summary_df = pickle.load(f)
    prefix = 'tet'
    
patient_ages = [44, 63, 53, 60, 34, 59, 51, 36, 35, 50, 44, 53, 60, 62, 38, 69, 32, 26, 46, 25, 53, 28, 36, 56, 37, 61, 66, 24, 33, 62]
event_times =  [18, 12, 42, 17, 5, 47, 30, 17, 19, 9, 9, 30, 21, 40, 15, 9, 10, 9, 24, 7, 47, 18, 13, 36, 15, 50, 57, 19, 7, 32]


plot_event_time_uncertainty(az_fit, event_times, patient_ages)


age_minus_event = np.array(patient_ages) - np.array(event_times)

posterior_t = az_fit.posterior['t'].values 
lower_ci = np.percentile(posterior_t, 2.5, axis=(0, 1))
upper_ci = np.percentile(posterior_t, 97.5, axis=(0, 1))
ci_widths = upper_ci - lower_ci

plt.figure(figsize=(7, 5))
scatter = plt.scatter(
    age_minus_event, 
    ci_widths, 
    c=patient_ages, 
    cmap='viridis', 
    alpha=0.8, 
    edgecolors='k', 
    s=70
)
plt.axvline(x=20, linestyle='--', color='black', linewidth=1, label='Non-identifiability Threshold\n($T_{after}=20$)')
plt.xlabel('Time since cn-LOH (years)', fontsize=13)
plt.ylabel('Credible Interval Width (years)', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.7)

cbar = plt.colorbar(scatter)
cbar.set_label('Patient Age (years)', fontsize=12)
plt.legend(fontsize=11, loc='lower right', bbox_to_anchor=(1, 0.05))
plt.savefig('/Users/finnkane/Desktop/Plots/mu_regime_exp.png', dpi=300, bbox_inches="tight")
plt.show()

