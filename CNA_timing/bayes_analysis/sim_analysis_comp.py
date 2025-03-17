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

type = 1

hier_cnloh_az = '/Volumes/files/az/hier/comp/loh_comp.nc'
hier_cnloh_sum = '/Volumes/files/az/hier/comp/loh_comp.pkl'
sim_loh_az = '/Volumes/files/az/pool/loh/mu_regime/norm(0,0.01)/0.075.nc'
sim_loh_sum = '/Volumes/files/az/pool/loh/mu_regime/norm(0,0.01)/0.075.pkl'


with open(hier_cnloh_sum, 'rb') as f:
    hier_sum = pickle.load(f)
hier_az = az.from_netcdf(hier_cnloh_az)

with open(sim_loh_sum, 'rb') as f:
    sim_sum = pickle.load(f)
sim_az = az.from_netcdf(sim_loh_az)



data = pd.read_csv('hier_data_comp.csv')
patient_ages = data['age']
data = data.drop(columns=['age'])

num_cna = []
for index, row in data.iterrows():
    num_cna.append(len(row))
patient_idx = []
patient_idx = list(chain.from_iterable([[i] * x for i, x in enumerate(num_cna)]))
patient_idx = [i + 1 for i in patient_idx]

data = pd.read_csv('hier_data_comp.csv')
patient_ages = list(np.array(data['age'].values))
event_times = list(np.array(data['event_times'].values))
mu_values = [0.0075] * len(event_times)



def plot_error_on_mean(fig_name, az_fit1, az_fit2):


    means1 = az_fit1.posterior['t'].mean(dim=("chain", "draw")).values
    stds1 = az_fit1.posterior['t'].std(dim=("chain", "draw")).values / np.sqrt(az_fit1.posterior['t'].size)
    
    means2 = az_fit2.posterior['t'].mean(dim=("chain", "draw")).values
    stds2 = az_fit2.posterior['t'].std(dim=("chain", "draw")).values / np.sqrt(az_fit2.posterior['t'].size)
    
    plt.figure(figsize=(7, 5))
    plt.errorbar(means1, np.arange(len(means1)), xerr=stds1, fmt='bo', label='Fit 1')
    plt.errorbar(means2, np.arange(len(means2)), xerr=stds2, fmt='go', label='Fit 2')
    
    plt.xlabel("Mean Event Time (Years)", fontsize=14)
    plt.ylabel("Index", fontsize=14)
    plt.title("Comparison of Individual Means and Errors", fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(alpha=0.3)
    plt.savefig(fig_name, dpi=300, bbox_inches="tight")
    plt.show()


plot_error_on_mean('t.png', sim_az, hier_az)
plot_mu_uncertainty('hier_gamma_comp.png',hier_az, mu_values, patient_ages)

# plot_distributions(az_fit, mu_values, gamma_values, 4000)

patient_num = 0
# plot_pair(az_fit,patient_num)



