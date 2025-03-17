import pandas as pd
from itertools import chain
import cmdstanpy
import arviz as az
import matplotlib.pyplot as plt
from state_evolve.diploid_evolution import ss_initialisation, state_initialisation
from state_evolve.cnLOH_inf import diploid_to_cnLOH
from tri_inf import diploid_to_trisomy
from tet_inf import diploid_to_tetraploidy
from diploid_cnLOH_evoln import diploid_to_cnLOH_prob_dist1
from diploid_trisomy_evoln import diploid_to_trisomy_prob_dist1
from diploid_tetraploidy_evoln import diploid_to_tetraploidy_prob_dist1
import numpy as np
import pandas as pd
from plot import hist_plot
from scipy.special import softmax
import random
import pickle
from redundant.trial_dip import cnloh_sim

data = pd.read_csv('hier_sim_data1.csv')
patient_ages = data['age']
data = data.drop(columns=['age'])

num_cna = []
for index, row in data.iterrows():
    num_cna.append(len(row))
patient_idx = []
patient_idx = list(chain.from_iterable([[i] * x for i, x in enumerate(num_cna)]))
patient_idx = [i + 1 for i in patient_idx]

num_sites = len(patient_idx) * [500]
J = len(patient_ages)
num_samples = len(patient_ages)
patient_ages = np.repeat(patient_ages, 1).tolist()
event_times = list(chain.from_iterable(data.values.tolist()))
print(event_times)


mu_values =[]
gamma_values=[]
np.random.seed(4)
mu_values = [0.0075]*len(patient_ages)
gamma_values = [0.0075]*len(patient_ages)

print(patient_ages)
# for i in range(max(patient_idx)):
#     mu_mean = np.random.normal(-4.402, 0.0100)
#     mu_sdd = np.random.normal(-0.0007, 0.4073)
#     mu_raw = np.random.normal(0, 1)
    
#     mu = max(np.exp(mu_mean + mu_sdd * mu_raw), 0.001)
#     gamma = max(mu, 0.001) 
#     mu_values.append(mu)
#     gamma_values.append(gamma)

# print(mu_values)
# print(np.array(mu_values).mean())

ground_states = []
for type in types:
    vals = []
    # if type == 1:
    #     for i in range(len(patient_ages)):
    #         noisy_beta_before, noisy_beta_after = diploid_to_cnLOH(mu, gamma, state_initialisation, num_sites[i], event_times[i], patient_ages[i])
    #         vals.append(noisy_beta_after)
    #         K=2
    #     # prefix = f'/Users/finnkane/Desktop/ICR/inf_plots/ss_cnloh/t={event_time}/'
    # elif type == 2:
    #     for i in range(len(patient_ages)):
    #         noisy_beta_before, noisy_beta_after = diploid_to_trisomy(mu, gamma, state_initialisation, num_sites[i], event_times[i], patient_ages[i])
    #         K = 3
    #     # prefix = f'/Users/finnkane/Desktop/ICR/inf_plots/ss_tri/t={event_time}/'
    # elif type == 3:
    #     for i in range(len(patient_ages)):
    #         noisy_beta_before, noisy_beta_after = diploid_to_tetraploidy(mu, gamma, state_initialisation, num_sites[i], event_times[i], patient_ages[i])
    #         K = 4

    #     # prefix = f'/Users/finnkane/Desktop/ICR/inf_plots/ss_tet/t={event_time}/'
    #         vals.append(noisy_beta_after)

    # if type == 4:
    #     for i, x in enumerate(patient_idx):
    #         noisy_beta_after = diploid_to_cnLOH(mu_values[x-1], gamma_values[x-1],state_initialisation, num_sites[i], event_times[i], patient_ages[i])
    #         vals.append(noisy_beta_after)
    #         prefix1 = 'loh_1cna'  
    #         prefix2 = 'loh_1cna'
    #         K = 2
    #         print(type)
    #     # prefix = f'/Users/finnkane/Desktop/ICR/inf_plots/dip_cnloh/t={event_time}/'
    # elif type == 5:
    #     for i, x in enumerate(patient_idx):
    #         noisy_beta_after = diploid_to_trisomy(mu_values[x-1], gamma_values[x-1],state_initialisation, num_sites[i], event_times[i], patient_ages[i])
    #         vals.append(noisy_beta_after)
    #         prefix1 = 'tri_1cna'
    #         prefix2 = 'tri_1cna'
    #         print(type)
    #         K = 3
    # elif type == 6:
    #     for i, x in enumerate(patient_idx):
    #         noisy_beta_after = diploid_to_tetraploidy(mu_values[x-1], gamma_values[x-1],state_initialisation, num_sites[i], event_times[i], patient_ages[i])
    #         vals.append(noisy_beta_after)
    #         K = 4
    #         prefix1 = 'tet_1cna'
    #         prefix2 = 'tet_1cna'
    #         print(type)
    vals = np.array(vals)  
    vals = np.concatenate(vals)
    model = cmdstanpy.CmdStanModel(stan_file='hier_sim.stan')
    data = {
        'K': K, 
        'J': num_cna,
        'P': max(patient_idx),
        'patient_idx': patient_idx,
        'n': num_sites, 
        'y': vals,  
        'age': patient_ages,
        'type': type,
    }
    init_values = lambda: {
    "mu_mean": np.random.normal(-4.3, 0.2),
    # "mu_sdd": np.random.uniform(0.001, 1),
    # "mu_raw": np.random.uniform(-1, 1, P).tolist(),
    # "gamma_raw": np.random.uniform(0.5, 1.5, P).tolist(),
    "t_raw": np.random.uniform(0.01, 0.99, sum(num_cna)).tolist(),
    "kappa": np.random.uniform(10, 40, sum(num_cna)).tolist(),
    "eta": np.random.uniform(0.8, 1.0, sum(num_cna)).tolist(),
    "delta": np.random.uniform(0.0, 0.2, sum(num_cna)).tolist(),
    }

    # init_values = [
    #     {
    #         't': [
    #             np.random.uniform(0.01 * patient_ages[j], 0.99 * patient_ages[j])
    #             for j in range(J)
    #         ]
    #     }
    #     for _ in range(4) 
    # ]

    fit = model.sample(
        data=data,
        iter_sampling=1000,  
        iter_warmup=1000,  
        adapt_delta=0.9,
        max_treedepth=12,
        inits=init_values(),
        show_console=False
    )
    summary_df = fit.summary()

    # print(filtered_summary)
    # print(event_times)

    az_fit = az.from_cmdstanpy(
        fit,
        posterior_predictive="y_rep",  
        observed_data={"y": vals}      
    )

    az_fit.to_netcdf(f"{prefix1}.nc")

    with open(f'{prefix2}.pkl', 'wb') as f:
        pickle.dump(summary_df, f)
# filtered_summary = summary_df[summary_df.index.str.contains(r"^(t(\[\d+\]|_raw\[\d+\])?|mu(\[\d+\])?|gamma(\[\d+\])?)$", regex=True)]
