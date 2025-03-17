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
import pickle

mu_vals = [0.0075]
kappa_vals = [35]
# gamma = 0.01
# mu=0.01
# mu_vals = [0.02]
# event_times_dict = {
#     # 5: [40, 8, 20, 30, 31],
#     # 10: [23, 36, 51, 42, 25, 30, 15, 23, 6, 17],
#     # 15: [49, 66, 13, 22, 12, 26, 72, 24, 13, 22, 15, 8, 24, 40, 42],
#     # 20: [33, 39, 14, 29, 45, 18, 22, 17, 43, 46, 29, 22, 7, 9, 16, 20, 45, 33, 35, 12],
#     # 25: [25, 18, 10, 21, 57, 31, 44, 60, 19, 25, 14, 9, 29, 46, 5, 37, 55, 44, 19, 16, 8, 44, 22, 44, 55],
#     30: [18, 12, 42, 17, 5, 47, 30, 17, 19, 9, 9, 30, 21, 40, 15, 9, 10, 9, 24, 7, 47, 18, 13, 36, 15, 50, 57, 19, 7, 32]
# }

event_times_dict = {
    # 1: [22,30,27,33,16,27,28,18,16,25,20,24,29,24,15,31,14,13,20,10,29,13,17,23,17,29,28,12,15,29],
    # 2: [15,58,25,13,23,11,27,5,7,26,29,27,27,25,9,21,12,19,26,8,29,11,11,35,27,43,18,10,18,44],
    # 3: [11,24,11,5,29,54,23,31,23,8,29,48,28,57,11,60,17,9,24,7,20,18,31,13,5,14,60,15,10,45],
    # 4: [24,58,5,15,7,10,35,24,17,15,5,9,14,8,14,53,27,16,30,11,5,13,19,51,13,42,30,8,28,57],
    5: [39,8,48,5,29,54,9,8,20,5,5,29,5,57,5,64,12,9,41,4,37,23,15,37,27,56,10,5,28,47]
}

patient_ages_dict = {
    # 1: [44,63,53,60,34,59,51,36,35,50,44,53,60,62,38,69,32,26,46,25,53,28,36,56,37,61,66,24,33,62],
    # 2: [44,63,53,60,34,59,51,36,35,50,44,53,60,62,38,69,32,26,46,25,53,28,36,56,37,61,66,24,33,62],
    # 3: [44,63,53,60,34,59,51,36,35,50,44,53,60,62,38,69,32,26,46,25,53,28,36,56,37,61,66,24,33,62],
    # 4: [44,63,53,60,34,59,51,36,35,50,44,53,60,62,38,69,32,26,46,25,53,28,36,56,37,61,66,24,33,62],
    5: [44,63,53,60,34,59,51,36,35,50,44,53,60,62,38,69,32,26,46,25,53,28,36,56,37,61,66,24,33,62]
}


# patient_ages_dict = {
#     # 5: [75, 30, 52, 37, 62],
#     # 10: [53, 65, 75, 70, 64, 34, 44, 52, 33, 46],
#     # 15: [80, 76, 46, 54, 36, 55, 80, 57, 47, 54, 48, 45, 57, 72, 74],
#     # 20: [58, 65, 41, 38, 70, 45, 48, 70, 69, 57, 54, 49, 22, 33, 43, 46, 70, 59, 60, 38],
#     # 25: [56, 36, 33, 55, 60, 67, 71, 73, 54, 52, 65, 35, 65, 58, 39, 55, 66, 78, 54, 52, 55, 71, 61, 54, 69],
#     30: [44, 63, 53, 60, 34, 59, 51, 36, 35, 50, 44, 53, 60, 62, 38, 69, 32, 26, 46, 25, 53, 28, 36, 56, 37, 61, 66, 24, 33, 62]
# }
mu = 0.0075
# event_times_30 = np.array(event_times_dict[25])
# ages_30 = np.array(patient_ages_dict[25])

# mean_ratio_30 = np.mean(event_times_30)/np.mean(ages_30) 
# var_event_30 = np.var(event_times_30, ddof=1)
# var_age_30 = np.var(ages_30, ddof=1)
# print(mean_ratio_30)

# data = pd.read_csv('patient_data.csv')
# event_times = data['CNA Time']        
# patient_ages = data['Patient Age']
num_sites = [40]

for type in types:
    for kappa in kappa_vals:
        for num in num_sites:
            for key in event_times_dict.keys():
                event_times = event_times_dict[key]
                patient_ages = patient_ages_dict[key]
                num_sites_arr = [num] * len(patient_ages)

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
                #     for i in range(len(patient_ages)):
                #         noisy_beta_after = diploid_to_cnLOH_prob_dist1([0.5, 0, 0.5],kappa, mu, mu, event_times[i], patient_ages[i], num_sites_arr[i])
                #         vals.append(noisy_beta_after)
                #         prefix1 = f'kap{kappa}_4'
                #         prefix2 = f'kap{kappa}_4'
                #         # m, k, w = states
                #         # ground_states.extend([2] * m + [1] * k + [0] * w)
                #         K = 2
                #     # prefix = f'/Users/finnkane/Desktop/ICR/inf_plots/dip_cnloh/t={event_time}/'
                # if type == 5:
                #     for i in range(len(patient_ages)):
                #         noisy_beta_after = diploid_to_trisomy_prob_dist1([0.5, 0, 0.5], mu, mu, event_times[i], patient_ages[i], num_sites_arr[i])
                #         vals.append(noisy_beta_after)
                #         # vals = np.array(vals)  
                #         # print('sss', vals.mean())
                #         prefix1 = f'num_sites_comp_real'
                #         prefix2 = f'num_sites_comp_real'
                #         K = 3
                #         # prefix = f'/Users/finnkane/Desktop/ICR/inf_plots/dip_tri/t={event_time}/'
                # if type == 6:
                #     for i in range(len(patient_ages)):
                #         noisy_beta_after = diploid_to_tetraploidy_prob_dist1([0.5,0,0.5],kappa, mu, mu, event_times[i], patient_ages[i], num_sites_arr[i])
                #         K = 4
                #         prefix1 = f'kap{kappa}_6'
                #         prefix2 = f'kap{kappa}_6'
                # #     # prefix = f'/Users/finnkane/Desktop/ICR/inf_plots/dip_tet/t={event_time}/'
                #         vals.append(noisy_beta_after)
                # print(vals)
                # vals = np.array(vals)  
                vals = np.concatenate(vals)

                model = cmdstanpy.CmdStanModel(stan_file='multi_inf.stan')
                data = {
                    'K': K, 
                    'J': len(patient_ages),
                    'n': num_sites_arr, 
                    'y': vals,  
                    'age': patient_ages,
                    'type': type,
                }
                df = pd.DataFrame({'y': data['y']})
                print(df.describe())  
                print(df[df.isna().any(axis=1)]) 
                print(df[df['y'] < 0])  
                print((np.array(event_times)/np.array(patient_ages)).mean())

                init_values = [
                    {
                        't': [
                            np.random.uniform(0.01 * patient_ages[j], 0.99 * patient_ages[j])
                            for j in range(len(patient_ages))
                        ]
                    }
                    for _ in range(4) 
                ]
                fit = model.sample(
                    data=data,
                    iter_sampling=1000,  
                    iter_warmup=1000,  
                    adapt_delta=0.8,
                    max_treedepth=12,
                    inits=init_values,
                    show_console=False
                )
                summary_df = fit.summary()

                filtered_summary = summary_df[summary_df.index.str.contains(r"^(t(\[\d+\]|_raw\[\d+\])?|mu(\[\d+\])?|gamma(\[\d+\])?)$", regex=True)]
                print(filtered_summary)
                print(event_times)

                az_fit = az.from_cmdstanpy(
                    fit,
                    posterior_predictive="y_rep",  
                    observed_data={"y": vals}      
                )
                with open(f'{prefix2}.pkl', 'wb') as f:
                    pickle.dump(summary_df, f)

                az_fit.to_netcdf(f"{prefix1}.nc")






