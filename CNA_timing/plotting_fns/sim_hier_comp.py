import pickle
import pandas as pd
import arviz as az 
import matplotlib.pyplot as plt
import numpy as np  
import sys
import os
import tqdm
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
# from plot import convergence_check, plot_event_time_uncertainty, plot_ppc, posterior_single_mu_gamma

# sim_file_az = "/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)/0.0075.nc"
# sim_file_sum = "/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)/0.0075.pkl"
# hier_file_az = "/Volumes/files/az/hier/comp/tet_comp.nc"
# hier_file_sum = "/Volumes/files/az/hier/comp/tet_comp.pkl"
sim_file_az = "real_pooled_single.nc"
sim_file_sum = "real_pooled_single.pkl"
hier_file_az = 'cll_tri.nc'
hier_file_sum = 'cll_tri.pkl'

with open(sim_file_sum, 'rb') as f:
    summary_sim = pickle.load(f)

az_fit_sim = az.from_netcdf(sim_file_az)

with open(hier_file_sum, 'rb') as f:
    summary_hier = pickle.load(f)

az_fit_hier = az.from_netcdf(hier_file_az)

models = {'sim': az_fit_sim, 'hier': az_fit_hier}
loo_results = {}

with tqdm.tqdm(total=len(models), desc="Computing LOO") as pbar:
    for model_name, az_fit in models.items():
        loo_results[model_name] = az.loo(az_fit, pointwise=True)
        pbar.update(1)


print("\nLOO Results:")
for model_name, loo in loo_results.items():
    print(f"{model_name}: {loo}")

plt.figure(figsize=(10, 6))
for model_name, loo in loo_results.items():
    pareto_k_values = loo.pareto_k.values
    plt.hist(
        pareto_k_values, bins=20, alpha=0.5, edgecolor="black", label=model_name
    )

plt.axvline(0.7, color="red", linestyle="dashed", label="Threshold (0.7)")
plt.axvline(1.0, color="black", linestyle="dashed", label="Threshold (1.0)")
plt.xlabel("Pareto k values")
plt.ylabel("Count")
plt.legend()
plt.show()
#############################################################

# mu_summary_sim = summary_sim[summary_sim.index.str.contains('mu')]

# cred_int_sim = az.hdi(az_fit_sim, hdi_prob=0.95)
# lower_sim = cred_int_sim['mu'].sel(hdi="lower").values
# upper_sim = cred_int_sim['mu'].sel(hdi="higher").values
# mu_sim = mu_summary_sim['Mean'].values
# print(mu_sim)
# print(lower_sim)
# print(upper_sim)



# mu_hier_post = az_fit_hier.posterior['mu']
# credible_intervals = az.hdi(mu_hier_post, hdi_prob=0.95)
# lower_bounds = credible_intervals['mu'].sel(hdi="lower").values
# upper_bounds = credible_intervals['mu'].sel(hdi="higher").values
# posterior_means = mu_hier_post.mean(dim=("chain", "draw")).values 

# errors1 = posterior_means - 0.0075
# ME1 = np.mean(errors1)
# RMSE1 = np.sqrt(np.mean(errors1**2))

# print(f"Mean Error (ME): {ME1:.5f}")
# print(f"Root Mean Squared Error (RMSE): {RMSE1:.5f}")

# errors2 = mu_sim - 0.0075
# ME2 = np.mean(errors2)
# RMSE2 = np.sqrt(np.mean(errors2**2))

# print(f"Mean Error (ME): {ME2:.5f}")
# print(f"Root Mean Squared Error (RMSE): {RMSE2:.5f}")
# import numpy as np
# import matplotlib.pyplot as plt

# plt.figure(figsize=(10, 7))

# num_patients = len(lower_bounds)  
# patient_indices = np.arange(1, num_patients + 1)  

# for i, (lower, upper, posterior_mean, patient_idx) in enumerate(zip(lower_bounds, upper_bounds, posterior_means, patient_indices)):
#     plt.plot([lower, upper], [patient_idx, patient_idx], 'b-', alpha=0.7, linewidth=2, label='95% CI (hierarchical)' if i == 0 else None)
#     plt.plot(posterior_mean, patient_idx, 'ko', markersize=6.5, label=f"$\mu$ posterior mean\n(hierarchical)" if i == 0 else None)

# plt.axvspan(lower_sim, upper_sim, color='gray', alpha=0.3, label="95% CI (pooled)")

# plt.axvline(x=mu_sim, color='red', linestyle='--', linewidth=2, label=f"$\mu$ (pooled) = {float(mu_sim):.5f}")
# plt.axvline(x=0.0075, color='black', linestyle='-', linewidth=2, label="True $\mu$ = 0.0075")

# plt.xticks(fontsize=12)
# plt.ylim(0, 31)
# plt.yticks(np.arange(0, 31, 5), fontsize=12)  
# plt.xlabel("Methylation rate $\mu$ (per year)", fontsize=14.5)
# plt.ylabel("Patient Index", fontsize=14.5)
# plt.legend(fontsize=12, loc='upper right')
# # plt.xlim(0.004, 0.0095)
# plt.grid(alpha=0.3)

# plt.savefig('p_v_h_mu_tet.png', dpi=300, bbox_inches="tight")
# plt.show()

 