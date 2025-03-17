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

type = 1

hier_cnloh_az = '/Volumes/files/az/hier/loh/3.nc'
hier_cnloh_sum = '/Volumes/files/az/hier/loh/3.pkl'
hier_tri_az = '/Volumes/files/az/hier/tri.nc'
hier_tri_sum = '/Volumes/files/az/hier/tri.pkl'
hier_tet_az = 'az/hier/tet/tet_azn=100.pkl'
hier_tet_sum = 'az/hier/tet/tet_sumn=100.pkl'
if type == 1:
    with open(hier_cnloh_sum, 'rb') as f:
        summary_df = pickle.load(f)
    az_fit = az.from_netcdf(hier_cnloh_az)

    # az_fit = az.from_netcdf(hier_cnloh_az)
if type == 2:
    az_fit = az.from_netcdf(hier_tri_az)

    with open(hier_tri_sum, 'rb') as f:
        summary_df = pickle.load(f)
if type == 3:
    with open(hier_tet_az, 'rb') as f:
        az_fit = pickle.load(f)

    with open(hier_tet_sum, 'rb') as f:
        summary_df = pickle.load(f)
# print(summary_df)
# print(az_fit)
# print(max(az_fit.posterior["mu"].values.flatten()))

mu_values =[]
gamma_values=[]
np.random.seed(4)

# for i in range(10):
#     mu_mean = np.random.normal(-4.8502, 0.0100)
#     mu_sdd = np.random.normal(-0.0007, 0.9973)
#     mu_raw = np.random.normal(0, 1)
    
#     mu = max(np.exp(mu_mean + mu_sdd * mu_raw), 0.001)
#     gamma = max(mu, 0.001)  
    
#     mu_values.append(mu)
#     gamma_values.append(gamma)
# print(mu_values)


mu_summary = summary_df[summary_df.index.str.contains('mu')]
gamma_summary = summary_df[summary_df.index.str.contains('gamma')]
kappa_summary = summary_df[summary_df.index.str.contains(r'^kappa\[\d+\]$')]
eta_summary = summary_df[summary_df.index.str.contains(r'^eta\[\d+\]$')]
delta_summary = summary_df[summary_df.index.str.contains(r'^delta\[\d+\]$')]
t_summary = summary_df[summary_df.index.str.contains(r'^t\[\d+\]$')]

# print(mu_summary)
# print(gamma_summary)
# print(eta_summary)
# print(delta_summary)
# print(kappa_summary)
# print(kappa_summary['Mean'].mean())
# print(t_summary)
# plot_ppc(az_fit, 't.png')


data = pd.read_csv('hier_sim_data3.csv')
patient_ages = data['age']
data = data.drop(columns=['age'])

num_cna = []
for index, row in data.iterrows():
    num_cna.append(len(row))
patient_idx = []
patient_idx = list(chain.from_iterable([[i] * x for i, x in enumerate(num_cna)]))
patient_idx = [i + 1 for i in patient_idx]

# data = pd.read_csv('hier_sim_data5.csv')
# patient_ages = list(np.array(data['age'].values))
# event_times = list(np.array(data['event_times'].values))
# mu_values = [0.0075] *len(event_times)
# num_cna = []
# for index, row in data.iterrows():
#     num_cna.append(len(row))
# patient_idx = []
# patient_idx = list(chain.from_iterable([[i] * x for i, x in enumerate(num_cna)]))
# patient_idx = [i + 1 for i in patient_idx]

mu_values =[]
gamma_values=[]
num_sites = len(patient_idx) * [500]
J = len(patient_ages)
num_samples = len(patient_ages)
np.random.seed(4)

patient_ages = np.repeat(patient_ages, 1).tolist()
event_times = list(chain.from_iterable(data.values.tolist()))
for i in range(max(patient_idx)):
    mu_mean = np.random.normal(-4.402, 0.0100)
    mu_sdd = np.random.normal(-0.0007, 0.4073)
    mu_raw = np.random.normal(0, 1)
    
    mu = max(np.exp(mu_mean + mu_sdd * mu_raw), 0.001)
    gamma = max(mu, 0.001) 
    mu_values.append(mu)
    gamma_values.append(gamma)
# mu_values = [0.0075]*len(patient_ages)

print(len(patient_ages))
# plot_event_time_uncertainty_hier('hier_cred_loh_2.png',az_fit, event_times, patient_ages,2)
plot_mu_uncertainty('hier_mu_loh_3.png',az_fit, mu_values, patient_ages)

# plot_distributions(az_fit, mu_values, gamma_values, 4000)

patient_num = 0
# plot_pair(az_fit,patient_num)



