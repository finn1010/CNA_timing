import pickle
import pandas as pd
from scipy.special import softmax
import arviz as az 
import matplotlib.pyplot as plt
import numpy as np  
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from plot import plot_pair_pool, plot_mu_distribution, convergence_check, plot_event_time_uncertainty, plot_ppc

    

type = 1
cnloh_az = 'num_sites_comp_real.nc'
cnloh_sum = 'num_sites_comp_real.pkl'
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
# plot_ppc(az_fit,f'/Users/finnkane/Desktop/Plots/sim/loh_ppc0.01.png')




# # convergence_check(summary_df)
def posterior_mean_difference(az_fit, true_event_times):
    posterior_means = np.mean(az_fit.posterior['t'].values, axis=(0, 1))
    return posterior_means - true_event_times

# mu_summary = summary_df[summary_df.index.str.contains('mu')]
# gamma_summary = summary_df[summary_df.index.str.contains('gamma')]
kappa_summary = summary_df[summary_df.index.str.contains(r'^kappa\[\d+\]$')]
eta_summary = summary_df[summary_df.index.str.contains(r'^eta\[\d+\]$')]
delta_summary = summary_df[summary_df.index.str.contains(r'^delta\[\d+\]$')]
t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]

# print(mu_summary)
# print(gamma_summary)
print(eta_summary)
print(delta_summary)
print(kappa_summary)
print(kappa_summary['Mean'].mean())
print(t_summary)

patient_ages = [44, 63, 53, 60, 34, 59, 51, 36, 35, 50, 44, 53, 60, 62, 38, 69, 32, 26, 46, 25, 53, 28, 36, 56, 37, 61, 66, 24, 33, 62]
event_times =  [18, 12, 42, 17, 5, 47, 30, 17, 19, 9, 9, 30, 21, 40, 15, 9, 10, 9, 24, 7, 47, 18, 13, 36, 15, 50, 57, 19, 7, 32]
d = posterior_mean_difference(az_fit, event_times)

print(d.mean())
convergence_check(summary_df, rhat_threshold=1.01, neff_threshold=400)

# plot_pair_pool('/Users/finnkane/Desktop/Plots/sim/pairs/tetpairs.png', az_fit, 3)
# plot_mu_distribution('/Users/finnkane/Desktop/Plots/sim/loh_mu0.015.png', az_fit, true_mu=0.015, J=4000)

# plot_ppc(az_fit, 's.png')



# # posterior_single_mu_gamma(az_fit)
# # patient_num = 1
plot_event_time_uncertainty('loh_cred.png', az_fit, event_times, patient_ages)
