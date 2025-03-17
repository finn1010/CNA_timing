import arviz as az
import matplotlib.pyplot as plt
from state_evolve.diploid_evolution import ss_initialisation, state_initialisation
from cnLOH_inf import diploid_to_cnLOH
from tri_inf import diploid_to_trisomy
from tet_inf import diploid_to_tetraploidy
from diploid_cnLOH_evoln import diploid_to_cnLOH_prob_dist1
from diploid_trisomy_evoln import diploid_to_trisomy_prob_dist1
from diploid_tetraploidy_evoln import diploid_to_tetraploidy_prob_dist1
import numpy as np
import pandas as pd
from plot import hist_plot
import pickle
import xarray as xr

type = 1
cnloh_az = 'az/compare_mu/0.0005/0.0005.nc'
cnloh_sum = 'az/compare_mu/0.0005/0.0005.pkl'
tri_az = 'az/pool/tri/0.005/0.5/500.nc'
tri_sum = 'az/pool/tri/0.005/0.5/500.pkl'
tet_az =  'az/pool/tet/0.01/0.5/300.nc'
tet_sum = 'az/pool/tet/0.01/0.5/300.pkl'
      
patient_ages = [44, 63, 53, 60, 34, 59, 51, 36, 35, 50, 44, 53, 60, 62, 38, 69, 32, 26, 46, 25, 53, 28, 36, 56, 37, 61, 66, 24, 33, 62]

if type == 1:
    with open(cnloh_sum, 'rb') as f:
        summary_df = pickle.load(f)

    az_fit = az.from_netcdf(cnloh_az)

mu_summary = summary_df[summary_df.index.str.contains('mu')]
gamma_summary = summary_df[summary_df.index.str.contains('gamma')]
t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]

mu = mu_summary['Mean'].values
mu = np.array(mu)
mu = mu.item()  

gamma = gamma_summary['Mean'].values
gamma = np.array(gamma)
gamma = gamma.item() 

t_vals = list(t_summary['Mean'])
num_sites_arr = [1600] * len(t_summary)
vals = []
type = 4

if type == 4:
    for i in range(len(patient_ages)):
        noisy_beta_after = diploid_to_cnLOH_prob_dist1([0.5, 0, 0.5], mu, gamma, t_vals[i], patient_ages[i], num_sites_arr[i])
        vals.append(noisy_beta_after)
vals = np.concatenate(vals)
posterior_pred = az_fit.posterior_predictive["y_rep"]


observed_data = {"y": np.array(vals)}  
posterior_predictive = {"y_rep": az_fit.posterior_predictive["y_rep"].values}
new_az_data = az.from_dict(observed_data=observed_data, posterior_predictive=posterior_predictive)
az.plot_ppc(new_az_data, data_pairs={"y": "y_rep"}, kind="kde", figsize=(8, 6))

plt.show()





# # plt.hist(pp_flat, bins=50, density=True, alpha=0.5, label="Posterior Predictive")

# # plt.hist(obs_values, bins=30, density=True, alpha=0.7, label="Observed Data", color="red")

# # plt.xlabel("y values")
# # plt.ylabel("Density")
# # plt.legend()
]
# # plt.show()