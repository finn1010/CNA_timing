posterior_pred = az_fit.posterior_predictive["y_rep"]

patient_id = 12
n_sites_per_patient = 300
start_idx = patient_id * n_sites_per_patient  
end_idx = (patient_id + 1) * n_sites_per_patient 

patient_slice = slice(start_idx, end_idx)

pp_subset = posterior_pred.isel(y_rep_dim_0=patient_slice)

obs_subset = az_fit.observed_data["y"].isel(y_dim_0=patient_slice)
# print(np.shape(obs_subset))

pp_flat = pp_subset.values.flatten() 

obs_values = obs_subset.values  

plt.hist(pp_flat, bins=50, density=True, alpha=0.5, label="Posterior Predictive")
plt.hist(obs_values, bins=30, density=True, alpha=0.7, label="Observed Data", color="red")

# Labels and legend
plt.xlabel("y values")
plt.ylabel("Density")
plt.legend()
plt.title("Posterior Predictive Check for Patient 1")

plt.show()