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
from plot import plot_mu_uncertainty, plot_event_time_uncertainty_hier, plot_ppc, plot_pair, plot_distributions
import seaborn as sns
from scipy.stats import linregress
from itertools import chain
type = 1

hier_cnloh_az = '/Volumes/files/az/hier/tet/tet_hier5.nc'
hier_cnloh_sum = '/Volumes/files/az/hier/tet/tet_hier5.pkl'

if type == 1:
    with open(hier_cnloh_sum, 'rb') as f:
        summary_df = pickle.load(f)
    az_fit = az.from_netcdf(hier_cnloh_az)

data = pd.read_csv('hier_sim_data5.csv')
patient_ages = data['age']
data = data.drop(columns=['age'])

num_cna = []
for index, row in data.iterrows():
    num_cna.append(len(row))
patient_idx = []
patient_idx = list(chain.from_iterable([[i] * x for i, x in enumerate(num_cna)]))
patient_idx = [i + 1 for i in patient_idx]


mu_values =[]
gamma_values=[]
num_sites = len(patient_idx) * [500]
J = len(patient_ages)
num_samples = len(patient_ages)
np.random.seed(4)

patient_ages = np.repeat(patient_ages, 5).tolist()
event_times = list(chain.from_iterable(data.values.tolist()))
for i in range(max(patient_idx)):
    mu_mean = np.random.normal(-4.402, 0.0100)
    mu_sdd = np.random.normal(-0.0007, 0.4073)
    mu_raw = np.random.normal(0, 1)
    
    mu = max(np.exp(mu_mean + mu_sdd * mu_raw), 0.001)
    gamma = max(mu, 0.001) 
    mu_values.append(mu)
    gamma_values.append(gamma)
mu_summary = summary_df[summary_df.index.str.contains('mu')]

t_summary = list(summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]['Mean'].values)
mu_summary = list(summary_df[summary_df.index.str.contains(r'^mu\[\d+\]$')]['Mean'].values)

mu_summary = np.repeat(mu_summary, 5).tolist()
mu_values = np.repeat(mu_values, 5).tolist()



mu_diff = np.array(mu_summary) - np.array(mu_values)
t_diff = np.array(t_summary) - np.array(event_times)

slope, intercept, r_value, p_value, std_err = linregress(mu_diff, t_diff)

t_diff_pred = slope * mu_diff + intercept

plt.figure(figsize=(8, 6))
sc = plt.scatter(mu_diff, t_diff, c=event_times, cmap='viridis', alpha=0.7, edgecolors='k')

plt.plot(mu_diff, t_diff_pred, color='red', linewidth=2, label=f'$R^2$={r_value**2:.3f}, p={p_value:.3g}')

cbar = plt.colorbar(sc)
cbar.set_label("Event time (years)", fontsize=14.5)

plt.xlabel("Error on methylation rate posterior mean (per year)",fontsize=14.5)
plt.ylabel("Error on event time posterior mean (years)",fontsize=14.5)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('tri_err_comp_mu_t5.png', dpi=300, bbox_inches="tight")

plt.show()

print(f"R²: {r_value**2:.3f}")
print(f"P-value: {p_value:.3g}")

