import matplotlib.pyplot as plt
import pickle
import arviz as az

patient_ages=[44,63,53,60,34,59,51,36,35,50,44,53,60,62,38,69,32,26,46,25,53,28,36,56,37,61,66,24,33,62]

event_times=[18,12,42,17,5,47,30,17,19,9,9,30,21,40,15,9,10,9,24,7,47,18,13,36,15,50,57,19,7,32]

cnloh_sum = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)/0.0075.pkl'
cnloh_az = '/Volumes/files/az/pool/tet/mu_regime/norm(0,0.01)/0.0075.nc'
az_fit = az.from_netcdf(cnloh_az)

posterior_median = az_fit.posterior.median(dim=["chain", "draw"])


with open(cnloh_sum, 'rb') as f:
    summary_df = pickle.load(f)

t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]

t_vals = list(t_summary['Mean'].values)

time_since_event = [age - event for age, event in zip(patient_ages, event_times)]

event_time_difference = [t_val - event for event, t_val in zip(event_times, t_vals)]

plt.figure(figsize=(8, 6))
scatter = plt.scatter(time_since_event, event_time_difference, 
                      c=patient_ages, cmap='viridis', alpha=0.7, edgecolor='k')

cbar = plt.colorbar(scatter)
cbar.set_label('Patient Age', fontsize=13)

plt.xlabel('Time since event (years)', fontsize=13)
plt.ylabel('Posterior mean - true event time (years)', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tick_params(axis='both', labelsize=12)
plt.axhline(y=0, color='r', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(f'/Users/finnkane/Desktop/Plots/sim/tet_acc.pdf', format='pdf', dpi=300)

plt.show()
